import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import RegularGridInterpolator

# import input parameters
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from inputs.input import InputParameters
from ForceField import FK_FF, force_interpolator, get_quantum_force

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
    cfg = InputParameters("../inputs/parameters.json")
    qs_path = os.path.join(cfg.result_dir, "rho_qs.csv")
    force_2d = FK_FF(cfg, qs_path)
    FF_interpolator = force_interpolator(force_2d, cfg.q_min, cfg.q_max, cfg.p_min, cfg.p_max, cfg.nbins)
    position, momentum = -1.375, 0.0
    test_force = get_quantum_force(FF_interpolator, position, momentum)
    print("Test Quantum Force:", test_force[0])

    quantum_momentums, quantum_positions = PP_quantum(cfg, FF_interpolator, position, momentum) 
    classical_momentums, classical_positions = PP_classical(cfg, position, momentum)
    plt.figure(figsize=(6,6))
    plt.plot(quantum_positions[:800], quantum_momentums[:800], label="Quantum PP_traj")
    plt.plot(classical_positions[:800], classical_momentums[:800], label="Classical PP_traj")
    plt.title("PP Trajectory Comparison")
    plt.xlabel("Position")
    plt.ylabel("Momentum")
    plt.legend()
    plt.savefig("PP_Comparison.png", dpi=600)
    print("Phase point dynamics trajector is saved as PP_Comparison.png")

