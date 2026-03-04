// force.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include "../inputs/input.h"

// Calculate the force in the Cartesian Coordiante
void car_force(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& car_pos,
    std::vector<std::vector<double>>& car_force_array
) {
    double mass = cfg.mass;
    double omegaP = cfg.omegaP;
    std::string type = cfg.potential_type;
    double force_HO=0.0, force_PES=0.0, force_ext=0.0; 
    size_t N = cfg.N;
    size_t P_plus_1 = cfg.P_plus_1;

    // Resize output arrays
    car_force_array.resize(N, std::vector<double>(P_plus_1));

    for (size_t i = 0; i < N; i++){
        for (size_t k = 0; k < P_plus_1; k++){
            if (type == "Harmonic") { // U(x) = \frac{1}{2} x^2
                force_PES = -car_pos[i][k];
            }
            else if (type=="Quartic") { // U(x) = \frac{1}{4} x^4
                force_PES = -std::pow(car_pos[i][k],3);
            }
            else if (type=="Morse") { // U(x) = 4 (1 - e^{-x/4})^2
                force_PES = -2 * std::exp(-car_pos[i][k]/4) * (1 - std::exp(-car_pos[i][k]/4));
            }
            else if (type=="DoubleWell") {
                force_PES = -2 * car_pos[i][k] * (car_pos[i][k] * car_pos[i][k] - 1); // U(x) = \frac{1}{2} (x^2 - 1)^2
                // force_PES = -4 * car_pos[i][k] * (car_pos[i][k] * car_pos[i][k] - 1); // U(x) = (x^2 - 1)^2
            }
            else {
                throw std::runtime_error("Undefined PES type: " + type);
            }

            if (k==0) {
                force_HO = -mass * omegaP * omegaP * (car_pos[i][0]-car_pos[i][1]);
                force_ext = 0.5 * force_PES / (P_plus_1-1);
            }
            else if (k < P_plus_1 - 1) {
                force_HO = -mass * omegaP * omegaP * (-car_pos[i][k-1]+2*car_pos[i][k]-car_pos[i][k+1]);
                force_ext = force_PES / (P_plus_1-1);
            }
            else if (k == P_plus_1 - 1){
                force_HO = -mass * omegaP * omegaP * (car_pos[i][P_plus_1-1]-car_pos[i][P_plus_1-2]);
                force_ext = 0.5 * force_PES / (P_plus_1-1);
            }
            car_force_array[i][k] = force_HO + force_ext;
        }
    }
}

// Transform the Cartesian force to Staging force for simulation
void stg_force(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& car_force_array,
    std::vector<std::vector<double>>& stg_force_array
) {
    size_t N = cfg.N;
    size_t P_plus_1 = cfg.P_plus_1;
    stg_force_array.resize(N, std::vector<double>(P_plus_1));

    for (size_t i = 0; i < N; ++i){
        for (size_t k = 0; k < P_plus_1; ++k){
            double stg_force = 0.0;

            if (k == 0){
                for (size_t j = 0; j < P_plus_1; ++j){
                    stg_force += car_force_array[i][j];
                }
                stg_force_array[i][k] = stg_force;
            }

            else if (k < P_plus_1-1){
                stg_force = car_force_array[i][k];
                stg_force += (k - 1) * stg_force_array[i][k-1] / k; 
                stg_force_array[i][k] = stg_force;
            }

            else if (k == P_plus_1-1){
                stg_force = 0.5 * (car_force_array[i][0] - car_force_array[i][P_plus_1 - 1]);
                for (size_t j = 1; j < P_plus_1 - 1; ++j){
                    stg_force += ((0.5*(P_plus_1 - 1) - j) / (P_plus_1 - 1)) * car_force_array[i][j];
                }
                stg_force_array[i][k] = stg_force;
            }
        }
    }
}

