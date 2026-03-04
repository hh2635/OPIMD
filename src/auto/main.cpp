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
#include "../sampling/OP_init.h"
#include "../sampling/OP_staging.h"
#include "../sampling/OP_force.h"
#include "../sampling/OP_LG.h"

int main(){
    // for sampling only as an intermediate test, not complied in the project main file. 
    SimulationConfig cfg = read_config("../inputs/parameters.json");
    std::string result_dir = cfg.result_dir;

    // namespace fs = std::filesystem;
    // if (fs::exists(result_dir)) {
    //     fs::remove_all(result_dir);
    // }
    // fs::create_directories(result_dir);
    // std::cout << result_dir << " is created." << std::endl;

    std::filesystem::path qs_path = result_dir + "qs_tuple.bin";
    std::ofstream(qs_path, std::ios::binary | std::ios::trunc).close();

    std::vector<std::vector<double>> car_pos = pos_sampling(cfg);
    std::vector<std::vector<double>> stg_pos = staging_forward(cfg, car_pos);
    std::vector<double> stg_mass = staging_mass(cfg);
    std::vector<std::vector<double>> stg_mom = mom_sampling(cfg);

    run_NVT_LG(cfg, stg_mass, car_pos, stg_pos, stg_mom, qs_path);

    std::cout << "Open Polymer Tuple -> " << std::filesystem::absolute(qs_path) << std::endl;

    return 0;
}