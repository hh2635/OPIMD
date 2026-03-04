// staging.h

#ifndef STAGING_H
#define STAGING_H

#include <vector>


// staging forward 
std::vector<std::vector<double>> staging_forward(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& car_pos);

// staging backward 
std::vector<std::vector<double>> staging_backward(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& stg_pos);

#endif  // STAGING_H