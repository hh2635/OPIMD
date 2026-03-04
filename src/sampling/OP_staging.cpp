// staging.cpp
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include "OP_init.h"

std::vector<std::vector<double>> staging_forward(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& car_pos)
{
    int P = cfg.P;
    int N = cfg.N;

    std::vector<std::vector<double>> stg_pos(N, std::vector<double>(P + 1, 0.0));

    for (int i = 0; i < N; ++i) {
        for (int k = 0; k <= P; ++k) {
            double y = 0.0;
            if (k == 0) {
                y = 0.5 * (car_pos[i][0] + car_pos[i][P]);
            } else if (k != P) {
                y = car_pos[i][k] - (k * car_pos[i][k + 1] + car_pos[i][0]) / (k + 1);
            } else {
                y = car_pos[i][0] - car_pos[i][P];
            }
            stg_pos[i][k] = y;
        }
    }
    return stg_pos;
}

std::vector<std::vector<double>> staging_backward(
    const SimulationConfig& cfg,
    const std::vector<std::vector<double>>& stg_pos)
{
    int P = cfg.P;
    int N = cfg.N;
    std::vector<std::vector<double>> car_pos(N, std::vector<double>(P + 1, 0.0));

    for (int i = 0; i < N; ++i) {
        for (int k = 0; k <= P; ++k) {
            double x = 0.0;
            if (k == 0) {
                x = stg_pos[i][0] + 0.5 * stg_pos[i][P];
            } else if (k != P) {
                x = stg_pos[i][0] + (0.5 * P - k) * stg_pos[i][P] / P;
                for (int l = k; l < P; ++l) {
                    x += (k) * stg_pos[i][l] / (l);
                }
            } else {
                x = stg_pos[i][0] - 0.5 * stg_pos[i][P];
            }
            car_pos[i][k] = x;
        }
    }
    return car_pos;
}