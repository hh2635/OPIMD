#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# import input parameters
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from inputs.input import InputParameters

def Wigner_marginal(cfg):
    hbar = cfg.hbar
    q_min, q_max = cfg.q_min, cfg.q_max
    s_min, s_max = cfg.s_min, cfg.s_max
    p_min, p_max = cfg.p_min, cfg.p_max
    nbins = cfg.nbins
    result_dir = cfg.result_dir

    bin_path = os.path.join(result_dir, "qs_tuple.bin")
    data = np.fromfile(bin_path, dtype=np.float64)
    tuples = data.reshape((-1, 2))
    tuples = np.vstack([tuples, tuples * np.array([1.0, -1.0])])
    q = tuples[:, 0]
    s = tuples[:, 1]

    # compute marginal distribution of q, assuming the independence 
    dist_q, q_edges = np.histogram(q, bins=nbins, range=(q_min, q_max), density=True)
    q_centers = 0.5 * (q_edges[:-1] + q_edges[1:])
    marginal_q = np.column_stack((q_centers, dist_q))
    q_path = os.path.join(result_dir, "q_distribution.csv")
    np.savetxt(q_path, marginal_q, delimiter=",", header="q,dist_q")
    print(f"Marginal q distribution saved to {q_path}")

    # compute marginal distribution of s, assuming the independence
    dist_s, s_edges = np.histogram(s, bins=nbins, range=(s_min, s_max), density=True)
    s_centers = 0.5 * (s_edges[:-1] + s_edges[1:])
    marginal_s = np.column_stack((s_centers, dist_s))
    s_path = os.path.join(result_dir, "s_distribution.csv")
    np.savetxt(s_path, marginal_s, delimiter=",", header="s,dist_s")
    print(f"Marginal s distribution saved to {s_path}")

    # compute the Gaussian distribution by the standard deviation of s
    sigma_s = np.std(s)
    s_grid = np.linspace(s_min, s_max, nbins)
    Gaussian = np.exp(-0.5 * (s_grid / sigma_s)**2) / (np.sqrt(2 * np.pi) * sigma_s)
    Gaussian_s = np.column_stack((s_grid, Gaussian))
    Gaussian_path = os.path.join(result_dir, "s_Gaussian.csv")
    np.savetxt(Gaussian_path, Gaussian_s, delimiter=",", header="s,dist_s")
    print(f"Gaussian with sampled standard deviation saved to {Gaussian_path}")

    # compute marginal momentum distribution via Fourier transform over s
    p_grid = np.linspace(p_min, p_max, nbins)
    dist_p = np.zeros(nbins, dtype=float)
    kernel = np.exp(1j * (p_grid[:, None] * s_centers[None, :]) / hbar)  # (p_n, s_bins)
    dist_p = np.trapz(kernel * dist_s[None, :], x=s_centers, axis=1).real 
    marginal_p = np.column_stack((p_grid, dist_p))
    p_path = os.path.join(result_dir, "p_distribution.csv")
    np.savetxt(p_path, marginal_p, delimiter=",", header="p,dist_p")
    print(f"Marginal p distribution saved to {p_path}")

    # compute Wigner function with the independence assumption
    Q_OP = dist_s[len(s_centers)//2]
    Wigner = (1 / (2 * np.pi * hbar)) * (1 / Q_OP) * np.outer(dist_q, dist_p) 
    Wigner_path = os.path.join(result_dir, "Wigner_marginal.csv")
    np.savetxt(Wigner_path, Wigner, delimiter=",")
    print(f"Wigner function with marginal distributions saved to {Wigner_path}")

    # plot marginal distribution of end-to-end distance s, and compare with Gaussian
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(s_centers, dist_s, width=(s_edges[1]-s_edges[0]), color="skyblue", label='OPIMD')
    ax.plot(s_grid, Gaussian, color="red", label='Gaussian', linewidth="2", linestyle='-')
    ax.set_xlabel('s')
    ax.set_ylabel('$g(s)$')
    ax.set_xlim(s_min, s_max)
    ax.tick_params(direction="in", which="both", right=True)
    ax.legend(frameon=False)
    plt.tight_layout()
    fluct_path = os.path.join(result_dir, "s_distribution.png")
    plt.savefig(fluct_path, dpi=600)
    print(f"End-to-end distance distribution plot saved to {fluct_path}")

    # plot Wigner function
    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [q_min, q_max, p_min, p_max]
    im = ax.imshow(Wigner.T, origin='lower', aspect='auto', extent=extent, interpolation='bicubic')
    ax.set_xlabel('q')
    ax.set_ylabel('p')
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    Wigner_path = os.path.join(result_dir, "Wigner_marginal.png")
    plt.savefig(Wigner_path, dpi=600)
    print(f"Wigner function plot saved to {Wigner_path}")


if __name__ == "__main__":
    # load input parameters
    cfg = InputParameters("../inputs/parameters.json")
    Wigner_marginal(cfg)