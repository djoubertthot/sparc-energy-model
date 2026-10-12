"""a0 from 3D-modelled rotation curves (Iorio et al. 2017, LITTLE THINGS in 3D) vs the SPARC curves.

For DDO 154, DDO 168 and NGC 2366, in both SPARC and LITTLE THINGS: the SPARC baryonic
components are rescaled to the Iorio distance (R ~ D, V ~ sqrt D) and the Iorio circular
velocity V_c (asymmetric-drift corrected, 3D inclination) is fitted with the acceleration law,
a0 free, Y_disk marginalised with a 0.1 dex prior, inclination marginalised within the quoted
3D uncertainty, distance within 5 %. The same fit is done on the SPARC curve with SPARC errors.
Writes results/littlethings3d.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KPC = fb.KPC
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")
A0S = np.logspace(-11, -9, 81)
LT = {"DDO154": ("ddo154", 3.7, 67.9, 1.1), "DDO168": ("ddo168", 4.3, 62.0, 2.6), "NGC2366": ("ngc2366", 3.4, 65.1, 4.2)}


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def fit(r_kpc, vobs, evobs, df, D_sparc, D, eD, inc, einc, ysed):
    """a0 profile with Y, D, i nuisances; baryons interpolated from the SPARC components rescaled to D."""
    sD = D / D_sparc
    rs = df.R.values * sD
    gas, disk, bul = df.Vgas.values * np.sqrt(sD), df.Vdisk.values * np.sqrt(sD), df.Vbul.values * np.sqrt(sD)
    prof = np.full(len(A0S), np.inf)
    for yd in ysed * 10 ** np.array([-0.2, -0.1, 0.0, 0.1, 0.2]):
        pen_y = (np.log10(yd / ysed) / 0.1) ** 2
        for kd in (-2, -1, 0, 1, 2):
            fd = 1 + kd * eD / D
            r = r_kpc * KPC * fd
            vg = np.interp(r_kpc * fd, rs * fd, gas) * 1e3 * np.sqrt(fd)
            vd = np.interp(r_kpc * fd, rs * fd, disk) * 1e3 * np.sqrt(fd)
            vb = np.interp(r_kpc * fd, rs * fd, bul) * 1e3 * np.sqrt(fd)
            vb2 = vg * np.abs(vg) + yd * vd * np.abs(vd) + 1.4 * yd * vb * np.abs(vb)
            gN = np.maximum(vb2, 1.0) / r
            for ki in (-2, -1, 0, 1, 2):
                i1 = np.clip(inc + ki * einc, 10, 89)
                s = np.sin(np.radians(inc)) / np.sin(np.radians(i1))
                v, e = vobs * 1e3 * s, np.maximum(evobs, 0.5) * 1e3 * s
                pen = pen_y + kd**2 + ki**2
                for j, a0 in enumerate(A0S):
                    chi2 = np.sum(((v - np.sqrt(gN * nu(gN / a0) * r)) / e) ** 2) + pen
                    if chi2 < prof[j]:
                        prof[j] = chi2
    j = int(np.argmin(prof))
    sc = max(prof[j] / len(r_kpc), 1.0)
    ok = prof <= prof[j] + sc
    return A0S[j], A0S[ok].min(), A0S[ok].max(), prof[j] / len(r_kpc)


L = ["| galaxy | source | D (Mpc) | i (deg) | N | R_last (kpc) | V_last (km/s) | a0 (1e-10) | 1 sigma | chi2/N |", "|---|---|---|---|---|---|---|---|---|---|"]
for sp, (fn, D_lt, i_lt, ei_lt) in LT.items():
    c = cat.loc[sp]
    df = fb.load_rotmod(ROOT / f"data/rotmod/{sp}_rotmod.dat")
    ysed = float(c.Y36) if pd.notna(c.Y36) else 0.5
    t = pd.read_csv(ROOT / f"data/external/littlethings3d/finalrot/{fn}_onlinetab.txt", sep=r"\s+", comment="#", header=None,
                    names=["Ras", "R", "Vr", "eVr", "Va", "eVa", "Vc", "eVc", "Vd", "eVd", "S", "eS"])
    a_s = fit(df.R.values, df.Vobs.values, df.eVobs.values, df, float(c.D_Mpc), float(c.D_Mpc), float(c.e_D_Mpc), float(c.Inc_deg), float(c.e_Inc_deg), ysed)
    a_l = fit(t.R.values, t.Vc.values, t.eVc.values, df, float(c.D_Mpc), D_lt, 0.05 * D_lt, i_lt, ei_lt, ysed)
    a_r = fit(t.R.values, t.Vr.values, t.eVr.values, df, float(c.D_Mpc), D_lt, 0.05 * D_lt, i_lt, ei_lt, ysed)
    L.append(f"| {sp} | SPARC | {c.D_Mpc} | {c.Inc_deg} +/- {c.e_Inc_deg} | {len(df)} | {df.R.values[-1]:.1f} | {df.Vobs.values[-1]:.1f} +/- {df.eVobs.values[-1]:.1f} | {a_s[0] * 1e10:.2f} | {a_s[1] * 1e10:.2f}-{a_s[2] * 1e10:.2f} | {a_s[3]:.1f} |")
    L.append(f"| {sp} | LT 3D, V_c (drift-corrected) | {D_lt} | {i_lt} +/- {ei_lt} | {len(t)} | {t.R.values[-1]:.1f} | {t.Vc.values[-1]:.1f} +/- {t.eVc.values[-1]:.1f} | {a_l[0] * 1e10:.2f} | {a_l[1] * 1e10:.2f}-{a_l[2] * 1e10:.2f} | {a_l[3]:.1f} |")
    L.append(f"| {sp} | LT 3D, V_rot (no drift correction) | {D_lt} | {i_lt} | {len(t)} | {t.R.values[-1]:.1f} | {t.Vr.values[-1]:.1f} | {a_r[0] * 1e10:.2f} | {a_r[1] * 1e10:.2f}-{a_r[2] * 1e10:.2f} | {a_r[3]:.1f} |")
    # direct comparison of the curves at common radii (SPARC rescaled to the LT distance)
    sD = D_lt / float(c.D_Mpc)
    vs = np.interp(t.R.values, df.R.values * sD, df.Vobs.values)
    m = (t.R.values >= df.R.values[0] * sD) & (t.R.values <= df.R.values[-1] * sD)
    L.append(f"|  | ratio V_c(LT) / V_obs(SPARC) at common radii: median {np.median(t.Vc.values[m] / vs[m]):.3f}, asymmetric drift at R_last: {t.Va.values[-1]:.1f} km/s | | | | | | | | |")
with open(ROOT / "results/littlethings3d.md", "w", encoding="utf-8") as fh:
    fh.write("# a0 from the LITTLE THINGS 3D curves (Iorio et al. 2017) vs SPARC (2026-10-10)\n\n" + "\n".join(L) + "\n")
print("\n".join(L))
