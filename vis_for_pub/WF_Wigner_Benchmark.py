import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

# For GroundState Wigner Functions
Wigner_Harmonic = np.loadtxt(f"../Benchmark/GroundState_Wigner_Harmonic.csv", delimiter=",")
Wigner_Quartic = np.loadtxt(f"../Benchmark/GroundState_Wigner_Quartic.csv", delimiter=",")
Wigner_Morse = np.loadtxt(f"../Benchmark/GroundState_Wigner_Morse.csv", delimiter=",")
Wigner_DoubleWell = np.loadtxt(f"../Benchmark/GroundState_Wigner_DoubleWell.csv", delimiter=",")
temp = "GroundState" 

# # For ThermalState Wigner Functions
# Wigner_Harmonic = np.loadtxt(f"../Benchmark/ThermalState_Wigner_Harmonic.csv", delimiter=",")
# Wigner_Quartic = np.loadtxt(f"../Benchmark/ThermalState_Wigner_Quartic.csv", delimiter=",")
# Wigner_Morse = np.loadtxt(f"../Benchmark/ThermalState_Wigner_Morse.csv", delimiter=",")
# Wigner_DoubleWell = np.loadtxt(f"../Benchmark/ThermalState_Wigner_DoubleWell.csv", delimiter=",")
# temp = "ThermalState" 

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

if temp == "GroundState":
    vmin, vmax = -0.05, 0.35
    frac_zero = (0.0 - vmin) / (vmax - vmin)
    cdict = {
    'red': [
        (0.0,       0.0, 0.0),   # v = vmin  -> 蓝 (R=0)
        (frac_zero, 1.0, 1.0),   # v = 0     -> 白 (R=1)
        (1.0,       1.0, 1.0)    # v = vmax  -> 红 (R=1)
    ],
    'green': [
        (0.0,       0.0, 0.0),
        (frac_zero, 1.0, 1.0),
        (1.0,       0.0, 0.0)
    ],
    'blue': [
        (0.0,       1.0, 1.0),   # v = vmin  -> 蓝 (B=1)
        (frac_zero, 1.0, 1.0),   # v = 0     -> 白 (B=1)
        (1.0,       0.0, 0.0)    # v = vmax  -> 红 (B=0)
    ]
    }
    cmap = LinearSegmentedColormap('blue_white_red', segmentdata=cdict, N=256)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    Wigner_DoubleWell[np.abs(Wigner_DoubleWell) < 0.001] = 0.0    
    im3 = axes[3].imshow(Wigner_DoubleWell.T, origin='lower', aspect='auto', extent=[q_min, q_max, p_min, p_max], interpolation='bicubic', cmap=cmap, norm=norm)
    axes[3].tick_params(axis='both', which='major', direction='out', bottom=True, top=False, left=True, right=False)
    axes[3].set_xlabel('q')
    axes[3].set_ylabel('p')
else:
    im3 = axes[3].imshow(Wigner_DoubleWell.T, origin='lower', aspect='auto', extent=[q_min, q_max, p_min, p_max], interpolation='bicubic', cmap="viridis")
    axes[3].tick_params(axis='both', which='major', direction='out', bottom=True, top=False, left=True, right=False)
    axes[3].set_xlabel('q')
    axes[3].set_ylabel('p')

plt.colorbar(im0, ax=axes[0])
plt.colorbar(im1, ax=axes[1])
plt.colorbar(im2, ax=axes[2])
plt.colorbar(im3, ax=axes[3])

plt.tight_layout()
plt.savefig(f"WF_Wigner_GrondState.png", dpi=600, bbox_inches='tight') # GroundState
# plt.savefig(f"WF_Wigner_ThermalState.png", dpi=600, bbox_inches='tight') # ThermalState
