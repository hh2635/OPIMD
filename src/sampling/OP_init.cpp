// sys_init.cpp
#include "OP_init.h"
#include <random>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
// #include <random>
// #include <vector>
// #include <cmath>
// #include "input.h"

// Generate the beads' position according to uniform distribution
std::vector<std::vector<double>> pos_sampling(const SimulationConfig& cfg) {
    int P = cfg.P;
    int N = cfg.N;
    double x_min = cfg.x_min; 
    double x_max = cfg.x_max; 

    std::random_device rd;
    // std::mt19937 gen(42); 
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(x_min, x_max);
    std::vector<std::vector<double>> car_pos(N, std::vector<double>(P+1, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < P+1; ++j) {
            car_pos[i][j] = dis(gen);  
        }
    }
    return car_pos;
}

// Calculate the staging mass for momentum sampling. 
// Each bead has different staging mass, therefore, they have different momentum distribuion. 
std::vector<double> staging_mass(const SimulationConfig& cfg) {
    int P = cfg.P;
    double mass = cfg.mass;
    std::vector<double> result(P + 1, mass);

    for (int k = 1; k < P; ++k) {
        result[k] = mass * (k+1) / k; 
    }

    result[P] = mass / P; 

    return result;
}

// Generate the correct momentum for each bead, according to its staging mass. 
std::vector<std::vector<double>> mom_sampling(const SimulationConfig& cfg) {
    int P = cfg.P;
    int N = cfg.N;
    double T = cfg.T;
    double kB = cfg.kB;

    std::vector<double> stg_mass = staging_mass(cfg);
    std::vector<double> stg_sigma(P + 1, 0.0);
    for (int k = 0; k <= P; ++k) {
        stg_sigma[k] = std::sqrt(stg_mass[k] * kB * T);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<std::vector<double>> mom_matrix(N, std::vector<double>(P + 1, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= P; ++j) {
            std::normal_distribution<> gaussian(0.0, stg_sigma[j]);
            mom_matrix[i][j] = gaussian(gen);
        }
    }
    return mom_matrix;
}

// save sampling results to csv file
void save_matrix_to_csv(const std::string& filename, const std::vector<std::vector<double>>& matrix) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j];
            if (j != row.size() - 1)
                file << ",";
        }
        file << "\n";
    }

    file.close();
    std::cout << "CSV file saved to: " << filename << std::endl;
}

// load sampling results from csv file
std::vector<std::vector<double>> load_matrix_from_csv(const std::string& filename) {
    std::vector<std::vector<double>> matrix;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            try {
                row.push_back(std::stod(cell));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Warning: unable to convert cell to double: " << cell << std::endl;
                row.push_back(0.0); // fallback value
            }
        }

        if (!row.empty()) {
            matrix.push_back(row);
        }
    }

    file.close();
    std::cout << "CSV file loaded from: " << filename << std::endl;
    return matrix;
}