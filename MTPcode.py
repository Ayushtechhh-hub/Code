import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# =================================================
# Fixed physical constants (from thesis)
# =================================================
Mw = 18.01528          # g/mol
Ci = 1.1               # mol/m^3
Co = 0.55              # mol/m^3
Hp = 62.5e-6           # paper thickness (m)
D_PLA_imperm = 1e-15   # impermeable PLA (m^2/s)

# =================================================
# Experimental data
# =================================================

# 1× Biopolymer-1
H_1x = np.array([57, 60, 62]) * 1e-6
WVTR_1x = np.array([3363, 2232, 1922])

# 2× Biopolymer-1
H_2x = np.array([72, 74, 79]) * 1e-6
WVTR_2x = np.array([2395, 1673, 1513])

# 3× Biopolymer-1
H_3x = np.array([73, 76, 79]) * 1e-6
WVTR_3x = np.array([2177, 1739, 1470])

# 4× Biopolymer-1 + PLA
H_PLA = np.array([64, 79, 81, 83]) * 1e-6
WVTR_PLA = np.array([1880, 452, 414, 212])

# =================================================
# Linear regression + diffusion extraction
# =================================================
def fit_linear(H, WVTR):
    slope, intercept, r, p, se = linregress(H, WVTR)
    return slope, intercept, r**2

def extract_D(slope):
    return Mw * (Ci - Co) / abs(slope)

s1, i1, r1 = fit_linear(H_1x, WVTR_1x)
s2, i2, r2 = fit_linear(H_2x, WVTR_2x)
s3, i3, r3 = fit_linear(H_3x, WVTR_3x)
sP, iP, rP = fit_linear(H_PLA, WVTR_PLA)

Ds_1x = extract_D(s1)
Ds_2x = extract_D(s2)
Ds_3x = extract_D(s3)
D_PLA_perm = extract_D(sP)

# =================================================
# Print results (Results section ready)
# =================================================
print("Extracted diffusion coefficients (linear regression):")
print(f"Biopolymer-1 Ds (1×) = {Ds_1x:.2e} m²/s")
print(f"Biopolymer-1 Ds (2×) = {Ds_2x:.2e} m²/s")
print(f"Biopolymer-1 Ds (3×) = {Ds_3x:.2e} m²/s")
print(f"PLA Ds (permeable) = {D_PLA_perm:.2e} m²/s")
print(f"PLA Ds (impermeable) = {D_PLA_imperm:.2e} m²/s")

# =================================================
# Fit lines
# =================================================
Hfit_1 = np.linspace(H_1x.min(), H_1x.max(), 100)
Hfit_2 = np.linspace(H_2x.min(), H_2x.max(), 100)
Hfit_3 = np.linspace(H_3x.min(), H_3x.max(), 100)
Hfit_P = np.linspace(H_PLA.min(), H_PLA.max(), 200)

WVTR_fit_1 = s1 * Hfit_1 + i1
WVTR_fit_2 = s2 * Hfit_2 + i2
WVTR_fit_3 = s3 * Hfit_3 + i3
WVTR_PLA_perm = sP * Hfit_P + iP
WVTR_PLA_imperm = Mw * (Ci - Co) / (Hp / D_PLA_imperm + Hfit_P / D_PLA_imperm)

# =================================================
# FIGURE 1: RESULTS SECTION (4 SUBPLOTS)
# =================================================
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

axs[0,0].scatter(H_1x*1e6, WVTR_1x, color="orange")
axs[0,0].plot(Hfit_1*1e6, WVTR_fit_1, '--', color="orange")
axs[0,0].set_title(f"1× Biopolymer-1 (R² = {r1:.2f})")
axs[0,0].set_ylabel("WVTR (g/m²·day)")
axs[0,0].grid(True)

axs[0,1].scatter(H_2x*1e6, WVTR_2x, color="blue")
axs[0,1].plot(Hfit_2*1e6, WVTR_fit_2, '--', color="blue")
axs[0,1].set_title(f"2× Biopolymer-1 (R² = {r2:.2f})")
axs[0,1].grid(True)

axs[1,0].scatter(H_3x*1e6, WVTR_3x, color="green")
axs[1,0].plot(Hfit_3*1e6, WVTR_fit_3, '--', color="green")
axs[1,0].set_title(f"3× Biopolymer-1 (R² = {r3:.2f})")
axs[1,0].set_xlabel("Thickness (µm)")
axs[1,0].set_ylabel("WVTR (g/m²·day)")
axs[1,0].grid(True)

axs[1,1].scatter(H_PLA*1e6, WVTR_PLA, color="red", label="Measured")
axs[1,1].plot(Hfit_P*1e6, WVTR_PLA_perm, '--', color="red", label="PLA permeable")
axs[1,1].plot(Hfit_P*1e6, WVTR_PLA_imperm, ':', color="black", label="PLA impermeable")
axs[1,1].set_title(f"4× Biopolymer-1 + PLA (R² = {rP:.2f})")
axs[1,1].set_xlabel("Thickness (µm)")
axs[1,1].legend()
axs[1,1].grid(True)

plt.suptitle("WVTR vs Coating Thickness — Individual Coating Regimes", fontsize=14)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# =================================================
# FIGURE 2: ALL SCENARIOS IN ONE GRAPH
# =================================================
plt.figure(figsize=(10,7))

plt.scatter(H_1x*1e6, WVTR_1x, color="orange", label="1× Biopolymer-1")
plt.plot(Hfit_1*1e6, WVTR_fit_1, '--', color="orange")

plt.scatter(H_2x*1e6, WVTR_2x, color="blue", label="2× Biopolymer-1")
plt.plot(Hfit_2*1e6, WVTR_fit_2, '--', color="blue")

plt.scatter(H_3x*1e6, WVTR_3x, color="green", label="3× Biopolymer-1")
plt.plot(Hfit_3*1e6, WVTR_fit_3, '--', color="green")

plt.scatter(H_PLA*1e6, WVTR_PLA, color="red", label="4× Biopolymer-1 + PLA")
plt.plot(Hfit_P*1e6, WVTR_PLA_perm, '--', color="red", label="PLA permeable")
plt.plot(Hfit_P*1e6, WVTR_PLA_imperm, ':', color="black", label="PLA impermeable")

plt.xlabel("Coating Thickness (µm)")
plt.ylabel("WVTR (g/m²·day)")
plt.title("WVTR vs Coating Thickness — All Coating Scenarios")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

