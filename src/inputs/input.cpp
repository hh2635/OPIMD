// input.cpp
#include <iostream>
#include <sstream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <cmath>
#include "input.h"

using json = nlohmann::json;

SimulationConfig read_config(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Unable to open " << filename << "\n";
        std::exit(1);
    }

    json config_json;
    file >> config_json;

    SimulationConfig cfg;

    // read parameters
    cfg.hbar = config_json["physical_constants"]["hbar"];
    cfg.kB   = config_json["physical_constants"]["kB"];
    cfg.mass = config_json["open_chain_parameters"]["mass"];
    cfg.P = config_json["open_chain_parameters"]["P"];
    cfg.P_plus_1 = cfg.P + 1; 
    cfg.ExDoF = config_json["NHC_parameters"]["ExDoF"];
    cfg.tau = config_json["NHC_parameters"]["tau"];
    cfg.dN = config_json["NHC_parameters"]["dN"];
    cfg.gamma = config_json["LG_parameters"]["gamma"];
    cfg.T = config_json["simulation_parameters"]["T"];
    cfg.beta = 1.0 / (cfg.kB * cfg.T);
    cfg.omegaP = std::sqrt(static_cast<double>(cfg.P)) / (cfg.beta * cfg.hbar);
    cfg.N = config_json["simulation_parameters"]["N"];
    cfg.x_min = config_json["simulation_parameters"]["x_min"];
    cfg.x_max = config_json["simulation_parameters"]["x_max"];
    cfg.dt = config_json["simulation_parameters"]["dt"];
    cfg.steps = config_json["simulation_parameters"]["steps"];
    cfg.save_interval = config_json["simulation_parameters"]["save_interval"];
    cfg.potential_type = config_json["simulation_parameters"]["potential_type"];
    std::ostringstream ospath;
    ospath << "../../results/" << cfg.potential_type << "_beta=" << std::fixed << std::setprecision(2) << cfg.beta << "/";
    cfg.result_dir = ospath.str();

    cfg.s_min = config_json["measurement_parameters"]["s_min"];
    cfg.s_max = config_json["measurement_parameters"]["s_max"];
    cfg.p_min = config_json["measurement_parameters"]["p_min"];
    cfg.p_max = config_json["measurement_parameters"]["p_max"];
    cfg.q_min = config_json["measurement_parameters"]["q_min"];
    cfg.q_max = config_json["measurement_parameters"]["q_max"];
    cfg.nbins = config_json["measurement_parameters"]["nbins"];
    
    return cfg;
}