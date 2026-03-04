#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# ---------- (I) Solve TISE for Weights and Wavefunctions ----------
def Weight_Wavefunction(x, dx, potential, beta, n, mass, hbar):
    """
    Solve TDSE for Weights and Wavefunctions using Finite Difference Method (FDM)
    """
    if potential == "Harmonic":
        V = 0.5 * (x**2)
    elif potential == "Quartic":
        V = 0.25 * (x**4)
    elif potential == "Morse":
        V = 4 * (1 - np.exp(-x/4))**2
    elif potential == "DoubleWell":
        # V = 0.5 * (x**2 - 1.0)**2
        V = (x**2 - 1.0)**2
    else:
        raise ValueError("Potential type is not defined.")

    # 1. Construct the Finite-Difference Hamiltonian Matrix
    T_term = hbar**2 / (mass * dx**2)
    diag = T_term + V
    off  = -0.5 * T_term * np.ones(x.size - 1)

    # 2. Solve for the weight and wavefunctions for n states
    eigvals, eigvecs = eigh_tridiagonal(diag, off, select='i', select_range=(0, n-1))
    weights = np.exp(-beta * eigvals)
    probs = weights / np.sum(weights)

    for j in range(eigvecs.shape[1]):
        norm_factor = np.sqrt(np.trapz(np.abs(eigvecs[:, j])**2, x))
        eigvecs[:, j] /= norm_factor

    return probs, eigvecs

# ---------- (II) Interpolator for wavefunction evaluation ----------
def interpolation(q, psi):
    return np.interp(q, x, psi, left=0.0, right=0.0)

# ---------- (III) Single State Wigner Function (q,p) Computation ----------
def Single_WignerFunction(qW, pW, s, hbar, psi):
    """
    Computes the Wigner function for a given wavefunction psi
    """
    # 1. Initialize Wigner function array
    single_Wigner = np.zeros((qW.size, pW.size))

    # 2. Pre-calculate the Fourier kernel exp(-i p s / ħ)
    phase = np.exp(1j * (pW[:, None] * s[None, :]) / hbar)

    # 3. main loop over position grid qW
    for idx, q in enumerate(qW):
        psi_minus = interpolation(q - s / 2, psi)
        psi_plus  = interpolation(q + s / 2, psi)
        base = psi_minus * np.conj(psi_plus)
        I = np.trapz(phase * base[None, :], s, axis=1) 
        single_Wigner[idx, :] = (I.real) / (2.0 * np.pi * hbar) 

    return single_Wigner

# ---------- (IV) Thermal Wigner Function (q,p) Computation ----------
def Thermal_WignerFunction(qW, pW, s, hbar, probs, eigvecs):
    """
    Computes the Thermal Wigner function as a weighted sum of single state Wigner functions
    """
    # 1. Initialize Thermal Wigner function array
    thermal_Wigner = np.zeros((qW.size, pW.size))

    # 2. Loop over each state to accumulate the weighted Wigner functions
    for j in range(len(probs)):
        single_Wigner = Single_WignerFunction(qW, pW, s, hbar, eigvecs[:, j])
        thermal_Wigner += probs[j] * single_Wigner

    np.savetxt(f"ThermalState_Wigner_{potential}.csv", thermal_Wigner, delimiter=",")
    print(f"Thermal Wigner function saved to ThermalState_Wigner_{potential}.csv")
    return thermal_Wigner

# ---------- (V) Visualize Wigner Function ----------
def plot_Wigner(qW, pW, Wigner, potential):
    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [qW.min(), qW.max(), pW.min(), pW.max()]
    im=ax.imshow(Wigner.T, origin='lower', aspect='auto', extent=extent, interpolation='bicubic')
    fig.colorbar(im, ax=ax)
    ax.set_xlabel('q')
    ax.set_ylabel('p')
    plt.tight_layout()
    plt.savefig(f"ThermalState_Wigner_{potential}.png", dpi=600)
    print(f"Thermal Wigner function plot saved to ThermalState_Wigner_{potential}.png")

# ---------- (VI) Main Execution ----------
if __name__ == "__main__":
    # 1. Parameters for wavefunction calculation
    mass, hbar, beta, n = 1.0, 1.0, 1.0, 10
    # x_min, x_max, gird = -5.0, 5.0, 1024
    x_min, x_max, gird = -4.0, 6.0, 1024 # For Morse Potential
    x = np.linspace(x_min, x_max, gird)
    dx = x[1] - x[0]
    # potential = "Harmonic"  
    # potential = "Quartic"
    potential = "DoubleWell"
    # potential = "DoubleWell"
    probs, eigvecs = Weight_Wavefunction(x, dx, potential, beta, n, mass, hbar)
    print("highest probs", probs[-1])

    # 2. Parameters for Wigner function calculation
    q_min, q_max = -5.0, 5.0
    # q_min, q_max = -4.0, 6.0 # For Morse Potential
    p_min, p_max = -5.0, 5.0
    s_min, s_max = -6.0, 6.0
    nbins = 129
    qW = np.linspace(q_min, q_max, nbins)
    pW = np.linspace(p_min, p_max, nbins)
    s = np.linspace(s_min, s_max, nbins)

    Wigner = Thermal_WignerFunction(qW, pW, s, hbar, probs, eigvecs)
    plot_Wigner(qW, pW, Wigner, potential)
