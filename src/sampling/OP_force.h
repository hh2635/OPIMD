// force.h
#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <string>
#include "../inputs/input.h"  

// Calcualte force in Cartesian Coordinate
void car_force(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& car_pos,
    std::vector<std::vector<double>>& car_force_array
);

// Calculate the force in staging coordinate by staging transform (chain rule) 
void stg_force(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& car_force_array,
    std::vector<std::vector<double>>& stg_force_array
);

#endif  // FORCE_H