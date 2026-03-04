#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# ---------- (I) Solve Time-Independent Schrodinger Equation in Hilbert Space (x-grid) by Finite Difference Method ----------
def Wavefunction(x, dx, potential, mass, hbar):
    """
    Solve time-independent Schrodinger equation using the Finite Difference method (FDM)
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

    # 2. Solve for the ground state (select_range=(0, 0))
    eigvals, eigvecs = eigh_tridiagonal(diag, off, select='i', select_range=(0, 0))
    E0  = eigvals[0]
    psi0 = eigvecs[:, 0]

    # 3. Normalize the wavefunction
    psi0 /= np.sqrt(np.trapz(np.abs(psi0)**2, x)) 

    return E0, psi0

# ---------- (II) Interpolator for wavefunction evaluation ----------
def interpolation(q, psi):
    return np.interp(q, x, psi, left=0.0, right=0.0)

# ---------- (III) Ground State Wigner Function W(p, q) Computation ----------
def WignerFunction(qW, pW, s, hbar, psi):
    """
    Computes the Wigner function for a given wavefunction psi
    """
    # 1. Initialize Wigner function array
    Wigner = np.zeros((qW.size, pW.size))

    # 2. Pre-calculate the Fourier kernel exp(-i p s / ħ)
    phase = np.exp(1j * (pW[:, None] * s[None, :]) / hbar)

    # 3. main loop over position grid qW
    for idx, q in enumerate(qW):
        psi_minus = interpolation(q - s / 2, psi)
        psi_plus  = interpolation(q + s / 2, psi)
        base = psi_minus * np.conj(psi_plus)
        I = np.trapz(phase * base[None, :], s, axis=1) 
        Wigner[idx, :] = (I.real) / (2.0 * np.pi * hbar) 

    np.savetxt(f"GroundState_Wigner_{potential}.csv", Wigner, delimiter=",")
    print(f"Ground Wigner function data saved to GroundState_Wigner_{potential}.csv")
    return Wigner

# ---------- (IV) Visualize Wigner Function ----------
def plot_Wigner(qW, pW, Wigner, potential):
    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [qW.min(), qW.max(), pW.min(), pW.max()]
    im=ax.imshow(Wigner.T, origin='lower', aspect='auto', extent=extent, interpolation='bicubic')
    fig.colorbar(im, ax=ax)
    ax.set_xlabel('q')
    ax.set_ylabel('p')
    plt.tight_layout()
    plt.savefig(f"GroundState_Wigner_{potential}.png", dpi=600)
    print(f"Ground Wigner function plot saved to GroundState_Wigner_{potential}.png")

# ---------- (V) Main Execution ----------
if __name__ == "__main__": 
    # 1. Parameters for wavefunction calculation
    mass, hbar = 1.0, 1.0
    x_min, x_max, gird = -5.0, 5.0, 1024
    # x_min, x_max, gird = -4.0, 6.0, 1024 # For Morse Potential
    x = np.linspace(x_min, x_max, gird)
    dx = x[1] - x[0]
    # potential = "Harmonic"  
    # potential = "Quartic"
    # potential = "Morse"
    potential = "DoubleWell"
    E0, psi0 = Wavefunction(x, dx, potential, mass, hbar)
    print(f"Ground State Energy E0: {E0}")

    # 2. Parameters for Wigner function calculation
    q_min, q_max = -5.0, 5.0
    # q_min, q_max = -4.0, 6.0 # For Morse Potential
    p_min, p_max = -5.0, 5.0
    s_min, s_max = -6.0, 6.0
    nbins = 129
    qW = np.linspace(q_min, q_max, nbins)
    pW = np.linspace(p_min, p_max, nbins)
    s = np.linspace(s_min, s_max, nbins)

    Wigner = WignerFunction(qW, pW, s, hbar, psi0)
    plot_Wigner(qW, pW, Wigner, potential)
