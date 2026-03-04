import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.signal import savgol_filter
from scipy.integrate import trapezoid

# 导入输入参数 (保持原逻辑)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../OPIMD_LG/src/')))
from inputs.input import InputParameters

def FK_1D(cfg, qs_path, position):
    """
    计算在特定位置 q=position 时，对应所有 p 轴网格点的量子有效力 F(q, p)。
    """
    rho_qs = np.loadtxt(qs_path, delimiter=",")

    q_min, q_max = cfg.q_min, cfg.q_max
    s_min, s_max = cfg.s_min, cfg.s_max
    p_min, p_max = cfg.p_min, cfg.p_max
    nbins = cfg.nbins
    hbar = cfg.hbar

    q_grid = np.linspace(q_min, q_max, nbins)
    p_grid = np.linspace(p_min, p_max, nbins)
    s_grid = np.linspace(s_min, s_max, nbins)
    idx_q = np.argmin(np.abs(q_grid - position))
    q0 = q_grid[idx_q]
    print(q0)

    # obtain the local quantum fluctuation
    s_prob = rho_qs[:, idx_q]
    s_prob = 0.5 * (s_prob + s_prob[::-1])
    if len(s_prob) > 11: 
        s_prob = savgol_filter(s_prob, window_length=11, polyorder=3)
        

    # Shifted potential
    eps = 1e-10
    q_plus = q0 + s_grid / 2
    q_minus = q0 - s_grid / 2

    # Potential Type
    if cfg.potential_type == "Harmonic":
        DU = np.where(np.abs(s_grid) < 1e-5, q0, (0.5 * q_plus**2 - 0.5 * q_minus**2) / (s_grid + eps))
    elif cfg.potential_type == "Quartic":
        DU = np.where(np.abs(s_grid) < 1e-5, q0**3, (0.25 * q_plus**4 - 0.25 * q_minus**4) / (s_grid + eps))
    elif cfg.potential_type == "DoubleWell":
        DU = np.where(np.abs(s_grid) < 1e-5, 4 * q0 * (q0**2 - 1), ((q_plus**2-1)**2 - (q_minus**2-1)**2) / (s_grid + eps))
    else:
        raise ValueError(f"Unknown potential type: {cfg.potential_type}")

    cos_kernel = np.cos(p_grid[:, np.newaxis] * s_grid / hbar)
    numerator_FT = np.trapezoid(DU * s_prob * cos_kernel, s_grid, axis=1)
    denominator_FT = np.trapezoid(s_prob * cos_kernel, s_grid, axis=1)
    quantum_force = - numerator_FT / (denominator_FT + 1e-12)
    
    return s_grid, s_prob, DU, p_grid, quantum_force

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(script_dir, "../OPIMD_LG/src/inputs/parameters.json")
    results_dir = os.path.join(script_dir, "../OPIMD_LG/results/Quartic_beta=8.00/") 

    cfg = InputParameters(json_path)
    qs_path = os.path.join(results_dir, "rho_qs.csv")

    test_pos = 1.0
    s_grid, s_prob, DU, p_axis, forces = FK_1D(cfg, qs_path, test_pos)
    # print(s_grid)
    # print(DU)
    
    fig, (ax_top, ax_bot) = plt.subplots(1, 2, figsize=(9.6, 4), layout='constrained')
    for ax in [ax_top, ax_bot]:
        ax.tick_params(direction='in', right=True, which='both')
    #     for spine in ax.spines.values():
    #         spine.set_linewidth(1.6)

    ax_top.plot(s_grid[32:97], DU[32:97], color='black', lw=1.5)
    ax_top.set_xlabel(r"$s$", fontweight='bold', fontsize=12)
    ax_top.set_ylabel(r"$D(\bar{q_0})$", fontweight='bold', fontsize=12)
    ax_top.set_xlim(-3, 3)

    ax_inset = inset_axes(ax_top, width="45%", height="45%", loc='upper center', borderpad=2)
    # for spine in ax_inset.spines.values():
    #     spine.set_linewidth(1.0)
    ax_inset.plot(s_grid[32:97], s_prob[32:97], color='skyblue', lw=1.2)
    ax_inset.fill_between(s_grid[32:97], s_prob[32:97], alpha=0.2, color='skyblue')
    ax_inset.set_xlim(-3, 3)
    ax_inset.set_xlabel(r"$s$", fontweight='bold', fontsize=12)
    ax_inset.set_ylabel(r"$\mathbb{P}(q_0, s)$", fontweight='bold', fontsize=12)
    ax_inset.tick_params(direction='in', right=True, which='both')

    # ax_bot.plot(p_axis, forces, color='red', lw=1.5)
    ax_bot.plot(p_axis[24:105], forces[24:105], color='red', lw=1.5)
    ax_bot.set_xlim(-3, 3)
    ax_bot.set_xlabel(r'$p$', fontweight='bold', fontsize=12)
    ax_bot.set_ylabel(r'$\dot{p}$', fontweight='bold', fontsize=12)

    plt.savefig("FK_QF.png", dpi=600, bbox_inches='tight')

    # fig, ax1 = plt.subplots(figsize=(8, 6))
    # ax1.plot(s_grid[32:97], DU[32:97], color='darkorange', lw=2, label="Derivative of Potential with quantum fluctuation (DU)")
    # ax1.set_xlabel(r"Relative Coordinate $s$")
    # ax1.set_ylabel(r"Effective Potential Gradient")
    # ax1.set_xlim(-3, 3)
    # ax_inset = inset_axes(ax1, width="40%", height="40%", loc='upper center', borderpad=2)
    # ax_inset.plot(s_grid[32:97], s_prob[32:97], color='royalblue', lw=1.5, label="local quantum fluctuation")
    # ax_inset.tick_params(axis='both', which='major', labelsize=8)
    # ax_inset.fill_between(s_grid[32:97], s_prob[32:97], alpha=0.2, color='royalblue')
    # ax_inset.set_xlim(-3,3)
    # plt.savefig("quantum_fluctuation_inset.png", dpi=600, bbox_inches='tight')

    # plt.figure(figsize=(8,6))
    # plt.plot(p_axis[24:105], forces[24:105], label=f'q={test_pos}')
    # plt.xlim(-3, 3)
    # plt.xlabel('Momentum p')
    # plt.ylabel('Effective Force F(p)')
    # plt.title(f'Quantum Effective Force in {cfg.potential_type} Potential')
    # plt.legend()
    # plt.savefig("quantum_effective_force.png")