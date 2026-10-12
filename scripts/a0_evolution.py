"""Predictions of the two hypotheses on the origin of a0 (docs/A0_ORIGIN.md).

H_A : a0(z) = c H(z) / (2 pi)          a0 tied to the expansion rate, grows with redshift
H_B : a0     = c sqrt(Lambda/3) / (2 pi) a0 tied to the cosmological constant, constant in time
H_C : a0 = 1.2e-10 fixed empirical, for reference

Observable consequences at fixed baryonic distribution (deep regime g = sqrt(g_bar a0)):
  BTFR  V_flat^4 = G M_b a0(z)   ->  Delta log V | M_b = 1/4 log(a0(z)/a0(0)),  Delta log M_b | V = -log(a0(z)/a0(0))
  RAR   Delta log g_obs | g_bar  = 1/2 log(a0(z)/a0(0))
  dark fraction inside a radius grows with a0.
Writes results/a0_evolution.md
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
c = 2.99792458e8
H0 = 70e3 / 3.0857e22
OM, OL = 0.3, 0.7
A0_OBS = 1.2e-10

hbar, kB = 1.0546e-34, 1.3806e-23
a0_A0 = c * H0 / (2 * np.pi)
a0_B = c * np.sqrt(OL) * H0 / (2 * np.pi)          # sqrt(Lambda/3) = H0 sqrt(Omega_Lambda)
lines = []
lines.append(f"a0 observed = {A0_OBS:.2e} m/s^2 (SPARC)")
lines.append(f"H_A today : c H0 / 2 pi = {a0_A0:.2e} (H0 = 70) -> ratio to observed {a0_A0 / A0_OBS:.2f}; with H0 = 73: {c * 73e3 / 3.0857e22 / (2 * np.pi):.2e}")
lines.append(f"H_B       : c sqrt(Lambda/3) / 2 pi = {a0_B:.2e} -> ratio {a0_B / A0_OBS:.2f}")
TU = hbar * A0_OBS / (2 * np.pi * c * kB)
TdS = hbar * np.sqrt(OL) * H0 / (2 * np.pi * kB)
lines.append(f"Unruh temperature at a0: {TU:.1e} K ; de Sitter temperature: {TdS:.1e} K ; ratio {TU / TdS:.2f} (Milgrom 1999: the MOND regime is where a body's Unruh temperature falls below the de Sitter temperature)")
lines.append("\n| z | H(z)/H0 | a0_A(z)/a0(0) | BTFR: Delta log V at fixed M_b | Delta log M_b at fixed V | RAR deep: Delta log g_obs | H_B |")
lines.append("|---|---|---|---|---|---|---|")
for z in (0.0, 0.5, 1.0, 1.5, 2.0, 2.3, 3.0):
    E = np.sqrt(OM * (1 + z) ** 3 + OL)
    r = E
    lines.append(f"| {z} | {E:.2f} | {r:.2f} | {0.25 * np.log10(r):+.3f} | {-np.log10(r):+.2f} | {0.5 * np.log10(r):+.3f} | 0 |")
lines.append("\nH_A predicts, at z = 2.3, baryonic masses 0.5 dex LOWER at fixed flat velocity than today, and accelerations 0.25 dex higher at fixed g_bar.")
lines.append("H_B predicts no evolution of either relation. The two are separated by the measured zero point of the baryonic Tully-Fisher relation at z ~ 1-2.3,")
lines.append("provided the rotation velocities are corrected for pressure support and the baryonic masses include the (large) gas fractions of those discs.")
with open(ROOT / "results/a0_evolution.md", "w", encoding="utf-8") as fh:
    fh.write("# a0: expansion-tied or vacuum-tied? Predictions (2026-10-10)\n\n" + "\n".join(lines) + "\n")
print("\n".join(lines))
