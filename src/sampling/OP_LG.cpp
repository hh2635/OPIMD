#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <cstdint>
#include <filesystem>
#include <sys/stat.h>

#include "../inputs/input.h"
#include "OP_init.h"
#include "OP_staging.h"
#include "OP_force.h"
// #include "energy.h"

void update_stg_pos(
    const SimulationConfig& cfg,
    std::vector<std::vector<double>>& stg_pos,
    std::vector<std::vector<double>>& stg_mom,
    const std::vector<double>& stg_mass
) {
    size_t N = cfg.N;
    size_t P_plus_1 = cfg.P_plus_1;
    double dt = cfg.dt;

    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < P_plus_1; ++k) {
            stg_pos[i][k] += 0.5 * dt * stg_mom[i][k] / stg_mass[k];
        }
    }
}

void update_stg_mom(
    const SimulationConfig& cfg,
    std::vector<std::vector<double>>& stg_mom,
    std::vector<std::vector<double>>& stg_force_array
) {
    size_t N = cfg.N;
    size_t P_plus_1 = cfg.P_plus_1;
    double dt = cfg.dt;
    
    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < P_plus_1; ++k) {
            stg_mom[i][k] += 0.5 * dt * stg_force_array[i][k];
        }
    }
}

void rescale_stg_mom_LG(
    const SimulationConfig& cfg,
    std::vector<std::vector<double>>& stg_mom, 
    const std::vector<double>& stg_mass
) {
    size_t N = cfg.N;
    size_t P = cfg.P;
    size_t P_plus_1 = cfg.P_plus_1;
    double dt = cfg.dt;
    double gamma = cfg.gamma;
    double kB = cfg.kB;
    double T = cfg.T;

    std::vector<double> sigma_vector(P_plus_1, 0.0);
    for (size_t j = 0; j <= P; ++j) {
        sigma_vector[j] = std::sqrt(2.0 * kB * T * gamma * stg_mass[j]);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> gaussian(0.0, 1.0);

    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < P_plus_1; ++k) {
            double r = gaussian(gen);
            double decay = std::exp(-gamma * dt);
            double fluct = std::sqrt((1.0 - std::exp(-2.0 * gamma * dt)) / (2.0 * gamma));
            stg_mom[i][k] = stg_mom[i][k] * decay + r * sigma_vector[k] * fluct;
        }
    }
}

// void append_stg_bin(const std::string& path,
//                            const std::vector<std::vector<double>>& stg_pos,
//                            size_t col_index)
// {
//     if (stg_pos.empty()) return;

//     const size_t N = stg_pos.size();
//     std::ofstream ofs(path, std::ios::binary | std::ios::app);
//     if (!ofs) throw std::runtime_error("Failed to open file for appending: " + path);

//     // 收集这一列到缓冲区，一次性写入
//     std::vector<double> buffer; buffer.reserve(N);
//     for (size_t i = 0; i < N; ++i) buffer.push_back(stg_pos[i][col_index]);

//     ofs.write(reinterpret_cast<const char*>(buffer.data()),
//               static_cast<std::streamsize>(N * sizeof(double)));
// }

void append_tuple_bin(const std::string& path,
                   const std::vector<std::vector<double>>& stg_pos,
                   size_t q_index,
                   size_t s_index)
{
    if (stg_pos.empty()) return;

    const size_t N = stg_pos.size();
    std::ofstream ofs(path, std::ios::binary | std::ios::app);
    if (!ofs)
        throw std::runtime_error("Failed to open file for appending: " + path);

    // 一次性写入 N 行，每行包含两个 double (q, s)
    std::vector<double> buffer;
    buffer.reserve(N * 2);

    for (size_t i = 0; i < N; ++i) {
        buffer.push_back(stg_pos[i][q_index]);  // q
        buffer.push_back(stg_pos[i][s_index]);  // s
    }

    ofs.write(reinterpret_cast<const char*>(buffer.data()),
              static_cast<std::streamsize>(buffer.size() * sizeof(double)));
}

void run_NVT_LG(
    const SimulationConfig& cfg,
    const std::vector<double>& stg_mass,
    std::vector<std::vector<double>>& car_pos,
    std::vector<std::vector<double>>& stg_pos,
    std::vector<std::vector<double>>& stg_mom,
    const std::string& QS_bin_file
) {
    size_t steps = cfg.steps;
    // double gamma_LG = cfg.gamma;
    // double time = -dt; 

    std::vector<std::vector<double>> car_force_array, stg_force_array;

    for (size_t step = 0; step <= steps; step++) {
        car_force(cfg, car_pos, car_force_array);
        stg_force(cfg, car_force_array, stg_force_array);
        update_stg_mom(cfg, stg_mom, stg_force_array);
        update_stg_pos(cfg, stg_pos, stg_mom, stg_mass);
        rescale_stg_mom_LG(cfg, stg_mom, stg_mass); // Langevin thermostat
        update_stg_pos(cfg, stg_pos, stg_mom, stg_mass);
        car_pos = staging_backward(cfg, stg_pos);
        
        car_force(cfg, car_pos, car_force_array);
        stg_force(cfg, car_force_array, stg_force_array);
        update_stg_mom(cfg, stg_mom, stg_force_array);

        // 每10步保存一次 staging 坐标
        if (step > steps / 2 && step % cfg.save_interval == 0) {
            append_tuple_bin(QS_bin_file, stg_pos, 0, cfg.P);
        }
    }
    
    std::cout << "q-s tuple (stg_pos[:,0], (stg_pos[:,P])) appended to: " << std::filesystem::absolute(QS_bin_file) << std::endl;
}

// int main(){
//     SimulationConfig cfg = read_config("../../inputs/parameters.json");
//     std::string result_dir = cfg.result_dir;

//     namespace fs = std::filesystem;
//     if (fs::exists(result_dir)) {
//         fs::remove_all(result_dir);
//     }
//     fs::create_directories(result_dir);
//     std::cout << result_dir << " is created." << std::endl;

//     std::filesystem::path q_path = result_dir + "q.bin";
//     std::filesystem::path s_path = result_dir + "s.bin";

//     std::ofstream(q_path, std::ios::binary | std::ios::trunc).close();
//     std::ofstream(s_path, std::ios::binary | std::ios::trunc).close();

//     std::vector<std::vector<double>> car_pos = pos_sampling(cfg);
//     std::vector<std::vector<double>> stg_pos = staging_forward(cfg, car_pos);
//     std::vector<double> stg_mass = staging_mass(cfg);
//     std::vector<std::vector<double>> stg_mom = mom_sampling(cfg);

//     run_NVT_LG(cfg, stg_mass, car_pos, stg_pos, stg_mom, q_path, s_path);

//     std::cout << "Done. q -> " << std::filesystem::absolute(q_path)
//                   << ", s -> " << std::filesystem::absolute(s_path) << "\n";

//     return 0;
// }