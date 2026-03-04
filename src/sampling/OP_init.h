// sys_init.h

#ifndef SYS_INIT_H
#define SYS_INIT_H

#include <vector>
#include "../inputs/input.h"  // 需要用到 SimulationConfig

// 返回 staging mass（P+1 个）
std::vector<double> staging_mass(const SimulationConfig& cfg);

// 位置采样：返回 N × (P+1) 的矩阵
std::vector<std::vector<double>> pos_sampling(const SimulationConfig& cfg);

// 动量采样：返回 N × (P+1) 的矩阵
std::vector<std::vector<double>> mom_sampling(const SimulationConfig& cfg);

// // 保存采样结果到csv文件
// void save_matrix_to_csv(const std::string& filename, const std::vector<std::vector<double>>& matrix);

// // 读取csv文件并加载到内存
// std::vector<std::vector<double>> load_matrix_from_csv(const std::string& filename);

#endif // SYS_INIT_H