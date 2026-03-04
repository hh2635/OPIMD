#ifndef NVT_LG_H
#define NVT_LG_H

#include <vector>
#include <string>
#include "../inputs/input.h"

// 更新 staging 坐标
void update_stg_pos(
    const SimulationConfig& cfg,
    std::vector<std::vector<double>>& stg_pos,
    std::vector<std::vector<double>>& stg_mom,
    const std::vector<double>& stg_mass
);

// 更新 staging 动量
void update_stg_mom(
    const SimulationConfig& cfg,
    std::vector<std::vector<double>>& stg_mom,
    std::vector<std::vector<double>>& stg_force_array
);

// Langevin 热化
void rescale_stg_mom_LG(
    const SimulationConfig& cfg,
    std::vector<std::vector<double>>& stg_mom, 
    const std::vector<double>& stg_mass
);

// // 保存一维数组为 bin
// void append_stg_bin(const std::string& path,
//                            const std::vector<std::vector<double>>& stg_pos,
//                            size_t col_index);

// 保存二维数组为 bin
void append_tuple_bin(const std::string& path,
                   const std::vector<std::vector<double>>& stg_pos,
                   size_t q_index,
                   size_t s_index);
// 运行主程序
void run_NVT_LG(
    const SimulationConfig& cfg,
    const std::vector<double>& stg_mass,
    std::vector<std::vector<double>>& car_pos,
    std::vector<std::vector<double>>& stg_pos,
    std::vector<std::vector<double>>& stg_mom,
    const std::string& QS_bin_file
);

#endif // NVT_LG_H