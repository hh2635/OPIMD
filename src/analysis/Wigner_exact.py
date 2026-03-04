#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# import input parameters
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from inputs.input import InputParameters

def rho_Wigner(cfg):
    """Compute Wigner function from binary or sampled data."""
    # read parameters from cfg
    hbar = cfg.hbar
    q_min, q_max = cfg.q_min, cfg.q_max
    s_min, s_max = cfg.s_min, cfg.s_max
    p_min, p_max = cfg.p_min, cfg.p_max
    nbins = cfg.nbins
    result_dir = cfg.result_dir

    # read (q,s) tuples from binary file
    bin_path = os.path.join(result_dir, "qs_tuple.bin")
    data = np.fromfile(bin_path, dtype=np.float64)
    tuples = data.reshape((-1, 2))
    tuples = np.vstack([tuples, tuples * np.array([1.0, -1.0])])
    q = tuples[:, 0]
    s = tuples[:, 1]

    # do statistics to rho(q,s), normalize by trace condition with q length, save and plot rho_qs
    H_counts, q_edges, s_edges = np.histogram2d(q, s, bins=[nbins, nbins], range=[[q_min, q_max], [s_min, s_max]])
    dq = np.diff(q_edges)
    s_idx = nbins // 2
    denom_counts = H_counts[:, s_idx].sum()
    rho_qs = H_counts / (denom_counts * dq[:, None])
    rho_path = os.path.join(result_dir, "rho_qs.csv")
    np.savetxt(rho_path, rho_qs, delimiter=",")
    print(f"Density matrix rho(q,s) saved to {rho_path}")

    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [q_min, q_max, s_min, s_max]
    im = ax.imshow(rho_qs.T, origin='lower', aspect='auto', extent=extent, interpolation='bicubic')
    ax.set_xlabel('q')
    ax.set_ylabel('s')
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    rho_path = os.path.join(result_dir, "rho_qs.png")
    plt.savefig(rho_path, dpi=600)
    print(f"Wigner function plot saved to {rho_path}")

    # calculate Wigner function W(q,p) via Fourier transform over s, save and plot Wigner function
    p_grid = np.linspace(p_min, p_max, nbins)
    s_centers = 0.5 * (s_edges[:-1] + s_edges[1:])
    W = np.zeros((nbins, nbins), dtype=float)

    for iq in range(nbins):
        counts_s = rho_qs[iq, :]                  # length s_bins
        # integrand shape (p_n, s_bins): counts_s * exp(+i * p * s)
        phase = np.exp(1j * (p_grid[:, None] * s_centers[None, :]) / hbar)  # (p_n, s_bins)
        integrand = counts_s[None, :] * phase
        Fp = np.trapz(integrand, x=s_centers, axis=1).real   # complex, length p_n
        # normalization: divide by (2*pi*hbar) to match Wigner convention used earlier
        W[iq, :] = Fp / (2.0 * np.pi * hbar)

    Wigner_path = os.path.join(result_dir, "Wigner_exact.csv")
    np.savetxt(Wigner_path, W, delimiter=",")
    print(f"Wigner function saved to {Wigner_path}")

    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [q_min, q_max, p_min, p_max]
    im = ax.imshow(W.T, origin='lower', aspect='auto', extent=extent, interpolation='bicubic')
    ax.set_xlabel('q')
    ax.set_ylabel('p')
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    Wigner_path = os.path.join(result_dir, "Wigner_exact.png")
    plt.savefig(Wigner_path, dpi=600)
    print(f"Wigner function plot saved to {Wigner_path}")
    
# def plot_Wigner(cfg):
#     """Plot Wigner function from saved CSV."""
#     q_min, q_max = cfg.q_min, cfg.q_max
#     p_min, p_max = cfg.p_min, cfg.p_max
#     result_dir = cfg.result_dir

#     # === 加载数据 ===
#     input_path = os.path.join(result_dir, "Wigner.csv")
#     W = np.loadtxt(input_path, delimiter=',')

#     fig, ax = plt.subplots(figsize=(6, 5))
#     extent = [q_min, q_max, p_min, p_max]
#     im = ax.imshow(W.T, origin='lower', aspect='auto', extent=extent, interpolation='bicubic')
#     ax.set_xlabel('q')
#     ax.set_ylabel('p')
#     fig.colorbar(im, ax=ax)
#     plt.tight_layout()

#     plt_path = os.path.join(result_dir, "Wigner.png")
#     plt.savefig(plt_path, dpi=600)
#     print(f"Wigner function plot saved to {plt_path}")


if __name__ == "__main__":
    # === 初始化参数 ===
    cfg = InputParameters("../inputs/parameters.json")  # 路径按需修改

    # === 计算与绘图 ===
    rho_Wigner(cfg)
    # plot_Wigner(cfg)
