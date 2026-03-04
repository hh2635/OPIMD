import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import RegularGridInterpolator

# import input parameters
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../OPIMD_LG/src/')))
from inputs.input import InputParameters


def FK_FF(cfg, qs_path):
    q_min, q_max = cfg.q_min, cfg.q_max
    s_min, s_max = cfg.s_min, cfg.s_max
    p_min, p_max = cfg.p_min, cfg.p_max
    nbins = cfg.nbins

    hbar = cfg.hbar
    q_grid = np.linspace(q_min, q_max, nbins)
    p_grid = np.linspace(p_min, p_max, nbins)
    s_grid = np.linspace(s_min, s_max, nbins)
    ds = s_grid[1] - s_grid[0]

    potential_type = cfg.potential_type

    rho_qs = np.loadtxt(qs_path, delimiter=",")
    force_field = np.zeros((nbins, nbins), dtype=float)

    for i in range(nbins):
        q0 = q_grid[i]
        s_prob = rho_qs[:, i]
        s_prob = savgol_filter(s_prob, window_length=11, polyorder=3)

        eps = 1e-10
        q_plus = q0 + s_grid / 2
        q_minus = q0 - s_grid / 2

        # Potential Type
        if potential_type == "Harmonic":
            DU = np.where(np.abs(s_grid) < 1e-4, q0, (0.5 * q_plus**2 - 0.5 * q_minus**2) / (s_grid + eps))
        elif potential_type == "Quartic":
            DU = np.where(np.abs(s_grid) < 1e-4, q0**3, (0.25 * q_plus**4 - 0.25 * q_minus**4) / (s_grid + eps))
        elif potential_type == "Morse":
            DU = np.where(np.abs(s_grid) < 1e-4, 2 * np.exp(-q0/4) * (1 - np.exp(-q0/4)), (4 * (1-np.exp(-q_plus/4))**2 - 4 * (1-np.exp(-q_minus/4))**2) / (s_grid + eps))
        elif potential_type == "DoubleWell":
            DU = np.where(np.abs(s_grid) < 1e-4, 4 * q0 * (q0 * q0 - 1), ((q_plus**2-1)**2 - (q_minus**2-1)**2) / (s_grid + eps))

        shifted = DU * s_prob
        phase = np.cos(p_grid[:, np.newaxis] * s_grid / hbar)
        enumerator_FT = np.trapz(shifted * phase, s_grid)
        denominator_FT = np.trapz(s_prob * phase, s_grid) 

        force_field[:, i] = -1 * enumerator_FT / (denominator_FT + 1e-10)

    return force_field 

def force_interpolator(force_field, q_min, q_max, p_min, p_max, nbins):
    q_grid = np.linspace(q_min, q_max, nbins)
    p_grid = np.linspace(p_min, p_max, nbins)
    FF_interpolator = RegularGridInterpolator((p_grid, q_grid), force_field, bounds_error=False, fill_value=None)
    return FF_interpolator

def get_quantum_force(FF_interpolator, position, momentum):
    # point = np.array([momentum, position])
    points = np.column_stack((momentum, position))
    force = FF_interpolator(points)
    return force

def PP_quantum(cfg, FF_interpolator, position, momentum):
    mass = cfg.mass
    dt = cfg.dt
    steps = cfg.steps
    momentums, positions = [momentum], [position]
    for step in range(steps):
        force = get_quantum_force(FF_interpolator, position, momentum)
        momentum += 0.5 * force[0] * dt 
        position += momentum * dt / mass
        force = get_quantum_force(FF_interpolator, position, momentum)
        momentum += 0.5 * force[0] * dt 
        momentums.append(momentum)
        positions.append(position)
    return momentums, positions

def PP_classical(cfg, position, momentum):
    mass = cfg.mass
    dt = cfg.dt
    steps = cfg.steps
    momentums, positions = [momentum], [position]
    for step in range(steps):
        if cfg.potential_type == "Harmonic":
            force = -1 * position
        elif cfg.potential_type == "Quartic":
            force = -1 * position**3
        elif cfg.potential_type == "Morse":
            force = -2 * np.exp(-position/4) * (1 - np.exp(-position/4))
        elif cfg.potential_type == "DoubleWell":
            force = -4 * position * (position**2 - 1)
        momentum += 0.5 * force * dt 
        position += momentum * dt / mass

        if cfg.potential_type == "Harmonic":
            force = -1 * position
        elif cfg.potential_type == "Quartic":
            force = -1 * position**3
        elif cfg.potential_type == "Morse":
            force = -2 * np.exp(-position/4) * (1 - np.exp(-position/4))
        elif cfg.potential_type == "DoubleWell":
            force = -4 * position * (position**2 - 1)
        momentum += 0.5 * force * dt 

        momentums.append(momentum)
        positions.append(position)
    return momentums, positions

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(script_dir, "../OPIMD_LG/src/inputs/parameters.json")
    results_dir = os.path.join(script_dir, "../OPIMD_LG/results/DoubleWell_beta=8.00/") 

    cfg = InputParameters(json_path)
    qs_path = os.path.join(results_dir, "rho_qs.csv")
    force_field = FK_FF(cfg, qs_path)
    FF_interpolator = force_interpolator(force_field, cfg.q_min, cfg.q_max, cfg.p_min, cfg.p_max, cfg.nbins)
    position, momentum = -1.375, 0.0
    test_force = get_quantum_force(FF_interpolator, position, momentum)
    print(test_force[0])

    quantum_momentums, quantum_positions = PP_quantum(cfg, FF_interpolator, position, momentum) 
    classical_momentums, classical_positions = PP_classical(cfg, position, momentum)

    plt.figure(figsize=(4.8, 4))
    # ax = plt.gca()
    # for spine in ax.spines.values():
    #     spine.set_linewidth(1.0)

    x_grid = np.linspace(-1.6, 1.6, 321)
    p_grid = np.linspace(-1.8, 1.8, 361)
    X, P = np.meshgrid(x_grid, p_grid)
    Energy = (P**2) / 2.0 + (X**2 - 1)**2
    levels = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    contour = plt.contour(X, P, Energy, levels=levels, colors='gray', alpha=0.2, linestyles='--')
    # plt.clabel(contour, inline=True, fontsize=8, fmt='E=%.1f')

    plt.plot(classical_positions[:800], classical_momentums[:800], lw=1.5, color="black")
    plt.plot(quantum_positions[:800], quantum_momentums[:800], lw=1.5, color="red")
    plt.xlabel(r"$q$")
    plt.ylabel(r"$p$")
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.savefig("QC_PP.png", dpi=600)