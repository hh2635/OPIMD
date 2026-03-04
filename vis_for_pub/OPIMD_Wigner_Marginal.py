import os
import numpy as np
import matplotlib.pyplot as plt

# beta = 8.00
beta = 1.00
Wigner_Harmonic = np.loadtxt(f"../OPIMD_LG/results/Harmonic_beta={beta:.2f}/Wigner_marginal.csv", delimiter=",")
Wigner_Quartic = np.loadtxt(f"../OPIMD_LG/results/Quartic_beta={beta:.2f}/Wigner_marginal.csv", delimiter=",")
Wigner_Morse = np.loadtxt(f"../OPIMD_LG/results/Morse_beta={beta:.2f}/Wigner_marginal.csv", delimiter=",")
Wigner_DoubleWell = np.loadtxt(f"../OPIMD_LG/results/DoubleWell_beta={beta:.2f}/Wigner_marginal.csv", delimiter=",")

q_min, q_max = -5.0, 5.0
Morse_min, Morse_max = -4.0, 6.0
p_min, p_max = -5.0, 5.0

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.ravel()

im0 = axes[0].imshow(Wigner_Harmonic.T, origin='lower', aspect='auto', extent=[q_min, q_max, p_min, p_max], cmap="viridis", interpolation='bicubic')
axes[0].tick_params(axis='both', which='major', direction='out', bottom=True, top=False, left=True, right=False)
axes[0].set_xlabel('q')
axes[0].set_ylabel('p')

im1 = axes[1].imshow(Wigner_Quartic.T, origin='lower', aspect='auto', extent=[q_min, q_max, p_min, p_max], cmap="viridis", interpolation='bicubic')
axes[1].tick_params(axis='both', which='major', direction='out', bottom=True, top=False, left=True, right=False)
axes[1].set_xlabel('q')
axes[1].set_ylabel('p')

im2 = axes[2].imshow(Wigner_Morse.T, origin='lower', aspect='auto', extent=[Morse_min, Morse_max, p_min, p_max], cmap="viridis", interpolation='bicubic')
axes[2].tick_params(axis='both', which='major', direction='out', bottom=True, top=False, left=True, right=False)
axes[2].set_xlabel('q')
axes[2].set_ylabel('p')

im3 = axes[3].imshow(Wigner_DoubleWell.T, origin='lower', aspect='auto', extent=[q_min, q_max, p_min, p_max], cmap="viridis", interpolation='bicubic')
axes[3].tick_params(axis='both', which='major', direction='out', bottom=True, top=False, left=True, right=False)
axes[3].set_xlabel('q')
axes[3].set_ylabel('p')

plt.colorbar(im0, ax=axes[0])
plt.colorbar(im1, ax=axes[1])
plt.colorbar(im2, ax=axes[2])
plt.colorbar(im3, ax=axes[3])

plt.tight_layout()
plt.savefig(f"OPIMD_Wigner_Marginal_beta={beta:.2f}.png", dpi=600, bbox_inches='tight')
