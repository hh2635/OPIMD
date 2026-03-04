import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mass, hbar = 1.0, 1.0
beta = 8.00
# beta = 1.00
# beta = 0.25
s_min, s_max, nbins = -6.0, 6.0, 129
# s_min, s_max, nbins = -3.0, 3.0, 129 # beta = 0.25
s_grid = np.linspace(s_min, s_max, nbins)
sigma = np.sqrt(beta * hbar * hbar / mass)
Classical_Gaussian = (1 / (np.sqrt(2 * np.pi * sigma * sigma))) * np.exp(-0.5 * s_grid * s_grid / (sigma * sigma))

OPIMD_Fluctuation_Harmonic = pd.read_csv(f"../OPIMD_LG/results/Harmonic_beta={beta:.2f}/s_distribution.csv")
OPIMD_Fluctuation_Quartic = pd.read_csv(f"../OPIMD_LG/results/Quartic_beta={beta:.2f}/s_distribution.csv")
OPIMD_Fluctuation_Morse = pd.read_csv(f"../OPIMD_LG/results/Morse_beta={beta:.2f}/s_distribution.csv")
OPIMD_Fluctuation_DoubleWell = pd.read_csv(f"../OPIMD_LG/results/DoubleWell_beta={beta:.2f}/s_distribution.csv")

OPIMD_Gaussian_Harmonic = pd.read_csv(f"../OPIMD_LG/results/Harmonic_beta={beta:.2f}/s_Gaussian.csv")
OPIMD_Gaussian_Quartic = pd.read_csv(f"../OPIMD_LG/results/Quartic_beta={beta:.2f}/s_Gaussian.csv")
OPIMD_Gaussian_Morse = pd.read_csv(f"../OPIMD_LG/results/Morse_beta={beta:.2f}/s_Gaussian.csv")
OPIMD_Gaussian_DoubleWell = pd.read_csv(f"../OPIMD_LG/results/DoubleWell_beta={beta:.2f}/s_Gaussian.csv")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.ravel()

# for ax in axes:
#     ax.tick_params(direction='in', right=True, which='both', width=1.2)
#     for spine in ax.spines.values():
#         spine.set_linewidth(1.6)

axes[0].bar(OPIMD_Fluctuation_Harmonic["# s"].to_numpy(), OPIMD_Fluctuation_Harmonic["dist_s"].to_numpy(), width=(s_grid[1]-s_grid[0]), color="skyblue", label='OPIMD')
# axes[0].plot(OPIMD_Gaussian_Harmonic["# s"].to_numpy(), OPIMD_Gaussian_Harmonic["dist_s"].to_numpy(), color = "red", label="Gaussian")
axes[0].plot(s_grid ,Classical_Gaussian, color = "black", linestyle="--", label="Classical", lw=1.5)
axes[0].tick_params(direction="in", which="both", right=True)
# axes[0].legend(frameon=False)
axes[0].set_xlabel(r'$y_{P+1}$')
axes[0].set_ylabel(r'$g(y_{P+1})$')
axes[0].set_xlim(s_min, s_max)
axes[0].set_ylim(0, 0.34) # beta = 8.00
# axes[0].set_ylim(0, 0.46) # beta = 1.00
# axes[0].set_ylim(0, 0.84) # beta = 0.25

axes[1].bar(OPIMD_Fluctuation_Quartic["# s"].to_numpy(), OPIMD_Fluctuation_Quartic["dist_s"].to_numpy(), width=(s_grid[1]-s_grid[0]), color="skyblue", label='OPIMD')
# axes[1].plot(OPIMD_Gaussian_Quartic["# s"].to_numpy(), OPIMD_Gaussian_Quartic["dist_s"].to_numpy(), color = "red", label="Gaussian")
axes[1].plot(s_grid ,Classical_Gaussian, color = "black", linestyle="--", label="Classical", lw=1.5)
axes[1].tick_params(direction="in", which="both", right=True)
# axes[1].legend(frameon=False)
axes[1].set_xlabel(r'$y_{P+1}$')
axes[1].set_ylabel(r'$g(y_{P+1})$')
axes[1].set_xlim(s_min, s_max)
axes[1].set_ylim(0, 0.34) # beta = 8.00
# axes[1].set_ylim(0, 0.46) # beta = 1.00
# axes[1].set_ylim(0, 0.84) # beta = 0.25

axes[2].bar(OPIMD_Fluctuation_Morse["# s"].to_numpy(), OPIMD_Fluctuation_Morse["dist_s"].to_numpy(), width=(s_grid[1]-s_grid[0]), color="skyblue", label='OPIMD')
# axes[2].plot(OPIMD_Gaussian_Morse["# s"].to_numpy(), OPIMD_Gaussian_Morse["dist_s"].to_numpy(), color = "red", label="Gaussian")
axes[2].plot(s_grid ,Classical_Gaussian, color = "black", linestyle="--", label="Classical", lw=1.5)
axes[2].tick_params(direction="in", which="both", right=True)
# axes[2].legend(frameon=False)
axes[2].set_xlabel(r'$y_{P+1}$')
axes[2].set_ylabel(r'$g(y_{P+1})$')
axes[2].set_xlim(s_min, s_max)
axes[2].set_ylim(0, 0.34) # beta = 8.00
# axes[2].set_ylim(0, 0.46) # beta = 1.00
# axes[2].set_ylim(0, 0.84) # beta = 0.25

axes[3].bar(OPIMD_Fluctuation_DoubleWell["# s"].to_numpy(), OPIMD_Fluctuation_DoubleWell["dist_s"].to_numpy(), width=(s_grid[1]-s_grid[0]), color="skyblue", label='OPIMD')
# axes[3].plot(OPIMD_Gaussian_DoubleWell["# s"].to_numpy(), OPIMD_Gaussian_DoubleWell["dist_s"].to_numpy(), color = "red", label="Gaussian")
axes[3].plot(s_grid ,Classical_Gaussian, color = "black", linestyle="--", label="Classical", lw=1.5)
axes[3].tick_params(direction="in", which="both", right=True)
# axes[3].legend(frameon=False)
axes[3].set_xlabel(r'$y_{P+1}$')
axes[3].set_ylabel(r'$g(y_{P+1})$')
axes[3].set_xlim(s_min, s_max)
axes[3].set_ylim(0, 0.34) # beta = 8.00
# axes[3].set_ylim(0, 0.46) # beta = 1.00
# axes[3].set_ylim(0, 0.84) # beta = 0.25

plt.tight_layout()
plt.savefig(f"OPIMD_Fluctuation_Marginal_beta={beta:.2f}.png", dpi=600, bbox_inches='tight')
