// input.h
#pragma once
#include <string>

struct SimulationConfig {
    double hbar, kB;
    double mass;
    int P, P_plus_1;
    int ExDoF, dN;
    double tau, gamma;
    double T, beta, omegaP, x_min, x_max, dt;
    double s_min, s_max, p_min, p_max, q_min, q_max;
    int N, steps, save_interval, nbins;
    std::string potential_type, result_dir;
};

SimulationConfig read_config(const std::string& filename);