"""Temporary mass from energy in transit, per SPARC galaxy and per cosmic epoch.

Checks the claim that energy passing through a galaxy (internal and external
sources), held temporarily as mass m = E/c^2, plus the galaxy's own mass
variation, builds up a drastic mass excess over the galaxy's life.

Resident energy of a flow is  E_res = P * tau  (power x residence time); it is a
steady-state quantity, not a time integral. Every factor below is chosen in
favour of the mechanism:
- starlight: L_bol = L[3.6] (twice the docs' convention), residence 3 R/c
  (free streaming x3 for dust reprocessing and scattering)
- supernova neutrinos: 3e46 J per SN, one SN per 100 Msun formed, residence R/c
- CMB, cosmic neutrino background (0.7 CMB), extragalactic background light
  (0.05 CMB): energy density x the full sphere of radius R_last
- cosmic rays, magnetic field, warm-gas thermal energy: 1 eV/cm^3 each over the
  full sphere (they are confined to the disc, so this overcounts)
Epoch scaling: CMB and CnuB x (1+z)^4; starlight, SNe, cosmic rays, thermal
x SFR(z)/SFR(0) from Madau & Dickinson 2014; stellar mass M(z) = fraction
formed by z; radius R(z) = R/(1+z).

Writes results/temporary_mass_budget.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
C = 2.99792458e8
LSUN = 3.828e26
EV_CM3 = 1.602e-19 * 1e6                     # 1 eV/cm^3 in J/m^3
U_CMB = 4.17e-14                              # J/m^3 today
U_CNUB = 0.7 * U_CMB
U_EBL = 0.05 * U_CMB
YR = 3.156e7
H0 = 70e3 / 3.0857e22
OM = 0.3

cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")


def budget(R_m, L_bol_W, SFR_Msun_yr):
    V = 4 / 3 * np.pi * R_m**3
    t_cross = R_m / C
    e = dict(
        starlight=L_bol_W * 3 * t_cross,
        sn_neutrinos=(SFR_Msun_yr / 100.0) * 3e46 / YR * t_cross,
        CMB=U_CMB * V, CnuB=U_CNUB * V, EBL=U_EBL * V,
        cosmic_rays=EV_CM3 * V, magnetic=EV_CM3 * V, gas_thermal=EV_CM3 * V,
    )
    return e


rows = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    n = f.stem.replace("_rotmod", "")
    c = cat.loc[n]
    R = df.R.values[-1] * KPC
    L36 = float(c.sparc_L36_1e9_Lsun) * 1e9
    L_bol = L36 * LSUN
    sfr = float(c.SFR_Msun_yr) if pd.notna(c.SFR_Msun_yr) else 1e-10 * L36
    vb2 = fb.baryonic2(df, 0.5, 1.0, 0.7)[-1]
    M_bar = max(vb2, 1.0) * R / G
    M_dyn = (df.Vobs.values[-1] * 1e3) ** 2 * R / G
    e = budget(R, L_bol, sfr)
    E_tot = sum(e.values())
    rows.append(dict(Galaxy=n, R_kpc=R / KPC, L36=L36, SFR=sfr, M_bar=M_bar / MSUN, M_missing=max(M_dyn - M_bar, 0) / MSUN,
                     m_temp=E_tot / C**2 / MSUN, **{k: v / C**2 / MSUN for k, v in e.items()}))
T = pd.DataFrame(rows)
T["ratio_bar"] = T.m_temp / T.M_bar
T["ratio_missing"] = T.m_temp / T.M_missing.replace(0, np.nan)

# ---- epoch dependence for the median galaxy -------------------------------------
med = T.sort_values("M_bar").iloc[len(T) // 2]


def sfr_ratio(z):  # Madau & Dickinson 2014 shape, normalised to z = 0
    f = lambda zz: (1 + zz) ** 2.7 / (1 + ((1 + zz) / 2.9) ** 5.6)
    return f(z) / f(0.0)


def mass_fraction_formed(z):  # cumulative cosmic SFR density up to z, normalised to today
    zz = np.linspace(z, 20, 2000)
    E = np.sqrt(OM * (1 + zz) ** 3 + 1 - OM)
    num = np.trapezoid(sfr_ratio(zz) / ((1 + zz) * E), zz)
    zz0 = np.linspace(0, 20, 4000)
    E0 = np.sqrt(OM * (1 + zz0) ** 3 + 1 - OM)
    den = np.trapezoid(sfr_ratio(zz0) / ((1 + zz0) * E0), zz0)
    return num / den


epochs = []
for z in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
    fm = mass_fraction_formed(z)
    R_z = med.R_kpc * KPC / (1 + z)
    L_z = med.L36 * LSUN * sfr_ratio(z) * max(fm, 1e-3)  # luminosity tracks star formation and the mass already formed
    sfr_z = med.SFR * sfr_ratio(z) * max(fm, 1e-3)
    e = budget(R_z, L_z, sfr_z)
    e["CMB"] *= (1 + z) ** 4
    e["CnuB"] *= (1 + z) ** 4
    e["EBL"] *= (1 + z) ** 2
    for k in ("cosmic_rays", "magnetic", "gas_thermal"):
        e[k] *= sfr_ratio(z)
    m_temp = sum(e.values()) / C**2 / MSUN
    M_bar_z = med.M_bar * max(fm, 1e-3)
    epochs.append(dict(z=z, age_Gyr=None, M_bar_z=M_bar_z, m_temp=m_temp, ratio=m_temp / M_bar_z,
                       dominant=max(e, key=e.get)))

# ---- 'integrated since birth' totals, i.e. what a 100 % retention would give ------
zz = np.linspace(0, 10, 20001)
E = np.sqrt(OM * (1 + zz) ** 3 + 1 - OM)
# CMB energy intercepted by the (shrinking) galaxy, per unit present-day area: int c u0 (1+z)^4 (1+z)^-2 dt
I_cmb = np.trapezoid((1 + zz) ** 2 / ((1 + zz) * E), zz) / H0
E_cmb_int = C * U_CMB * I_cmb * np.pi * (med.R_kpc * KPC) ** 2
# starlight: total energy ever radiated ~ 7e-4 of the rest energy of all stars ever formed (~1.6 x present stellar mass)
E_star_int = 7e-4 * 1.6 * (0.5 * med.L36 * MSUN) * C**2

with open(ROOT / "results/temporary_mass_budget.md", "w", encoding="utf-8") as fh:
    fh.write("# Temporary mass from energy in transit (2026-10-10)\n\n")
    fh.write(f"{len(T)} SPARC galaxies, inside R_last, every factor generous to the mechanism (see script docstring).\n\n")
    fh.write("## Today, per galaxy (medians over the sample)\n\n| source | m_temp / M_bary (median) | m_temp / M_missing (median) |\n|---|---|---|\n")
    for k in ["starlight", "sn_neutrinos", "CMB", "CnuB", "EBL", "cosmic_rays", "magnetic", "gas_thermal", "m_temp"]:
        fh.write(f"| {k} | {(T[k] / T.M_bar).median():.1e} | {(T[k] / T.M_missing.replace(0, np.nan)).median():.1e} |\n")
    fh.write(f"\nLargest m_temp / M_bary in the sample: {T.ratio_bar.max():.1e} ({T.loc[T.ratio_bar.idxmax(), 'Galaxy']}).\n")
    fh.write(f"\n## Through cosmic time, median galaxy ({med.Galaxy}, M_bary = {med.M_bar:.1e} Msun, R_last = {med.R_kpc:.1f} kpc)\n\n")
    fh.write("| z | M_bary(z) / Msun | m_temp(z) / Msun | m_temp / M_bary | dominant source |\n|---|---|---|---|---|\n")
    for r in epochs:
        fh.write(f"| {r['z']} | {r['M_bar_z']:.1e} | {r['m_temp']:.1e} | {r['ratio']:.1e} | {r['dominant']} |\n")
    fh.write("\n## If energy were retained for ever (100 % absorption, no re-emission)\n\n")
    fh.write(f"- all starlight ever emitted: {E_star_int / C**2 / MSUN:.2e} Msun = {E_star_int / (med.M_bar * MSUN * C**2):.1e} M_bary\n")
    fh.write(f"- all CMB intercepted since z = 10, galaxy shrinking as 1/(1+z), opacity 1: {E_cmb_int / C**2 / MSUN:.2e} Msun = {E_cmb_int / (med.M_bar * MSUN * C**2):.1e} M_bary; "
             f"with the measured millimetre opacity 1e-5: {1e-5 * E_cmb_int / (med.M_bar * MSUN * C**2):.1e} M_bary\n")
    fh.write(f"- required: M_missing / M_bary = {(T.M_missing / T.M_bar).median():.1f} (median)\n")

print(f"N={len(T)}  median m_temp/M_bary = {T.ratio_bar.median():.1e}  max = {T.ratio_bar.max():.1e} ({T.loc[T.ratio_bar.idxmax(), 'Galaxy']})")
print(f"median m_temp/M_missing = {T.ratio_missing.median():.1e}")
print("by source, median m_temp/M_bary:", {k: f"{(T[k] / T.M_bar).median():.1e}" for k in ["starlight", "sn_neutrinos", "CMB", "CnuB", "EBL", "cosmic_rays", "magnetic", "gas_thermal"]})
print(f"\nmedian galaxy {med.Galaxy}: M_bar={med.M_bar:.2e} Msun, R_last={med.R_kpc:.1f} kpc")
for r in epochs:
    print(f"  z={r['z']:3.1f}  M_bar(z)={r['M_bar_z']:.2e}  m_temp={r['m_temp']:.2e} Msun  ratio={r['ratio']:.1e}  dominant={r['dominant']}")
print(f"\nretained for ever: starlight {E_star_int / (med.M_bar * MSUN * C**2):.1e} M_bary ; CMB opacity 1: {E_cmb_int / (med.M_bar * MSUN * C**2):.1e} M_bary ; CMB opacity 1e-5: {1e-5 * E_cmb_int / (med.M_bar * MSUN * C**2):.1e} M_bary")
print(f"required: {(T.M_missing / T.M_bar).median():.1f} M_bary")
print("wrote results/temporary_mass_budget.md")
