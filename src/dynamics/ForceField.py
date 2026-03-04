import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import RegularGridInterpolator

# import input parameters
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
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

    # edge_mask = np.abs(q_grid) > 3.5
    # force_field[:, edge_mask] = -(q_grid[edge_mask]**3)

    # plt.pcolormesh(q_array, p_array, force_field, shading="auto", cmap="RdBu_r")
    plt.plot(p_grid, force_field[:,77])
    # plt.colorbar(label="Effective Force")
    plt.xlabel("q")
    plt.ylabel("p") 
    plt.savefig("FK_FF.png", dpi=600)

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

if __name__ == "__main__":
    cfg = InputParameters("../inputs/parameters.json")
    qs_path = os.path.join(cfg.result_dir, "rho_qs.csv")
    force_2d = FK_FF(cfg, qs_path)
    FF_interpolator = force_interpolator(force_2d, cfg.q_min, cfg.q_max, cfg.p_min, cfg.p_max, cfg.nbins)
    test_force = get_quantum_force(FF_interpolator, position=np.array([1.0, 1.0, 1.0, 1.0, 1.0]), momentum=np.array([-2.0, -1.0, 0.0, 1.0, 2.0]))
    print("Test Quantum Force:", test_force)