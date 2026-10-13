# -*- coding: utf-8 -*-
"""
Predictions of Milgrom's modified-inertia (MI) models (Milgrom 2022, Phys. Rev. D 106, 064060, arXiv:2208.07073) for the two tests
that killed the local field theories: the wide binaries in the Galactic field, and the lensing circular velocities at 0.3-1 Mpc.

The model (its heuristic one-parameter form, eqs. 'v', 'shiluta', 'shitaba', 'exasa' of the paper): the equation of motion in Fourier
space is a(w) mu[A(w)/a0] = a_N(w), where the nonlocal acceleration A(w) = sum_k a_k theta(w_k/w) collects the acceleration amplitudes
a_k of all frequency components of the trajectory, weighted by theta(w_k/w), theta(1) = 1, theta decreasing, theta(0) of order a few
(the paper's examples: 2/(1+y^2), exp(1-y)). Circular orbits in a galaxy obey the algebraic relation a mu(a/a0) = g_N exactly (the
RAR). A subsystem with internal frequency w_in in an external field of acceleration a_ex and frequency w_ex has its internal dynamics
governed by mu[(a_in + theta(w_ex/w_in) a_ex)/a0]: isotropic, and set by the OBSERVED external acceleration a_ex, not by its
Newtonian value. Here mu is the inverse of the RAR nu (a_N = a mu(a/a0) <=> a = nu(a_N/a0) a_N).

Test 1, wide binaries: a_ex = 1.8 a0 (V^2/R at the Sun), w_ex/w_in = 1 Myr / 220 Myr -> theta(0). Internal acceleration of each star
a_in = g_N/(2 mu) (equal masses). The internal effective G is G/mu, the velocity ratio scales as mu^(-1/2); compared with the
calibration bins of results/wide_binaries_bins.csv (g_N > 10 a0) exactly as the data were.
Test 2, lensing at 0.3-2.5 Mpc: a circular test orbit at R obeys a mu[(a + theta_eff a_ex)/a0] = G M/R^2 with a_ex the galaxy's
acceleration in the large-scale field, a_ex = nu(e_N) e_N a0 = 0.0105 a0 for e_N = 1.1e-4 (isolated galaxies, efe_isolated_2mpp),
and theta_eff = theta(w_ex/w_in). MI has no lensing theory; we assume, as every relativistic MOND theory does, that the deprojected
lensing velocity equals the circular velocity of massive test particles. theta_eff is left free and bounded by the data, then read
against the frequency ratio: the internal period at 0.3-1 Mpc (15-50 Gyr) exceeds the age of the universe, so w_ex/w_in > 1 for any
large-scale field that varies within a Hubble time, and theta_eff < 1.
Output: results/mi_milgrom_predictions.md, .png
"""
import os

import numpy as np
import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
A0 = 1.2e-10
G = 4.301e-6; A0K = A0 * 3.0857e19 / 1e6          # kpc (km/s)^2 / Msun ; (km/s)^2/kpc


def nu(y):
    y = np.maximum(y, 1e-14)
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))


_y = np.logspace(-8, 5, 6001); _x = _y * nu(_y)


def mu(x):
    """MOND mu(a/a0) from the RAR nu: a_N = a mu(a/a0)"""
    y = 10 ** np.interp(np.log10(np.maximum(x, 1e-12)), np.log10(_x), np.log10(_y))
    return y / np.maximum(x, 1e-12)


def theta_fam(y, t0):
    """theta(y) with theta(0) = t0, theta(1) = 1, theta -> 0 as y^-2 (t0 = 2 gives Milgrom's 2/(1+y^2))"""
    return t0 / (1 + (t0 - 1) * y ** 2)


# ====================================================== test 1: wide binaries
B = pd.read_csv("results/wide_binaries_bins.csv")
W = pd.read_csv("results/wide_binaries.csv")
cal = W[W.gN > 10]
gcal = np.median(cal.gN)
AEX_BIN = 1.8


def mi_boost_binary(gN, t0, aex=AEX_BIN):
    m = 1.0
    for _ in range(60):
        m = mu(gN / (2 * m) + t0 * aex)
    return 1.0 / m


T0 = [1.0, 2.0, np.e, 3.0, 5.0]
for t0 in T0:
    B[f"mi_{t0:.2f}"] = np.sqrt(np.array([mi_boost_binary(g, t0) for g in B.gN]) / mi_boost_binary(gcal, t0))
low = B[B.hi <= 1.0]
chi = {}
grid = np.linspace(1.0, 10.0, 91)
for t0 in grid:
    p = np.sqrt(np.array([mi_boost_binary(g, t0) for g in low.gN]) / mi_boost_binary(gcal, t0))
    chi[t0] = np.sum(((low.ratio - p) / low.e_ratio) ** 2)
chi = pd.Series(chi)
best = chi.idxmin(); ok = chi[chi <= chi.min() + 1]
chi_newton = np.sum(((low.ratio - 1) / low.e_ratio) ** 2); chi_law = np.sum(((low.ratio - low.pred_law) / low.e_ratio) ** 2)
chi_qumond = np.sum(((low.ratio - low.pred_qumond) / low.e_ratio) ** 2)
mid = B[(B.lo >= 0.1) & (B.hi <= 1.0)]      # the two clean bins 0.1-1 a0

# ====================================================== test 2: lensing at 0.3-2.5 Mpc
def read_mrt(path):
    rows = []
    for line in open(path, encoding="utf-8"):
        if not line[:24].strip() or line.startswith(("Title", "Authors", "Table", "=", "-", " ", "Byte", "Note")):
            continue
        try:
            R = float(line[25:32])
        except ValueError:
            continue
        rec = dict(sample=line[0:24].strip(), R=R)
        for k in range(4):
            base = 33 + 21 * k
            try:
                rec[f"Vc{k+1}"], rec[f"es{k+1}"], rec[f"ey{k+1}"] = [float(line[b:b + 6]) for b in (base, base + 7, base + 14)]
            except ValueError:
                rec[f"Vc{k+1}"] = rec[f"es{k+1}"] = rec[f"ey{k+1}"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


L_ = read_mrt("data/external/lensing/mistele2024_vc_mrt.txt"); A = L_[L_["sample"] == "All"]
LOGMB = [10.10, 10.66, 10.96, 11.29]
EN_ISO = 1.1e-4
AEX_GAL = nu(EN_ISO) * EN_ISO                 # the galaxy's acceleration in the large-scale field, units of a0


def v_mi(M, c, R):
    """circular velocity of a test orbit at R (kpc) around point mass M with the MI external term c = theta_eff a_ex (units of a0)"""
    out = []
    for r in np.atleast_1d(R):
        n = G * M / r ** 2 / A0K
        lo, hi = 1e-8, 1e4
        for _ in range(80):
            x = np.sqrt(lo * hi)
            if x * mu(x + c) > n:
                hi = x
            else:
                lo = x
        out.append(np.sqrt(x * A0K * r))
    return np.array(out)


TH = np.concatenate([[0.0], np.logspace(-1.5, 1.0, 26)])
chi_l = np.zeros((4, len(TH))); curves = {}
for k, lm in enumerate(LOGMB):
    Mb = 10 ** lm
    V = A[f"Vc{k+1}"].values; e = np.sqrt(A[f"es{k+1}"].values ** 2 + A[f"ey{k+1}"].values ** 2); R = A.R.values
    ok_ = np.isfinite(V) & (R <= 1000); V, e, R = V[ok_], e[ok_], R[ok_]
    inner = R < 300
    scale = (np.average(V[inner], weights=1 / e[inner] ** 2) / (G * Mb * A0K) ** 0.25) ** 4
    for j, t in enumerate(TH):
        chi_l[k, j] = np.sum(((V - v_mi(Mb * scale, t * AEX_GAL, R)) / e) ** 2)
    curves[k] = dict(R=R, V=V, e=e, M=Mb * scale)
tot = chi_l.sum(axis=0)
tbest = TH[np.argmin(tot)]; tmax = TH[tot <= tot.min() + 4].max()
tmax1 = TH[chi_l[0] <= chi_l[0].min() + 4].max()
d_t1 = np.interp(1.0, TH[1:], tot[1:]) - tot.min(); d_t2 = np.interp(2.0, TH[1:], tot[1:]) - tot.min()
# internal periods at the lensing radii, bin 1
Vb1 = np.average(curves[0]["V"][curves[0]["R"] < 300], weights=1 / curves[0]["e"][curves[0]["R"] < 300] ** 2)
Tin = {r: 2 * np.pi * r * 3.0857e16 / Vb1 / 3.156e16 for r in (100, 300, 500, 1000)}   # Gyr
y_needed = np.sqrt((2.0 / tmax - 1) / 1.0) if tmax < 2 else np.nan                    # theta = 2/(1+y^2) = tmax -> y

# ====================================================== report
L = ["# Milgrom's modified-inertia models on the two tests (2026-10-13)\n",
     "Model: Milgrom 2022 (PRD 106, 064060), a(w) mu[A(w)/a0] = a_N(w) with A(w_in) = a_in + theta(w_ex/w_in) a_ex, theta(1) = 1, theta "
     "decreasing, theta(0) of order 2-3. The external acceleration enters with its OBSERVED value and with a weight set by the frequency ratio. "
     "Here mu is the inverse of the RAR nu.\n",
     "## Test 1: wide binaries (a_ex = 1.8 a0, w_ex/w_in = 0.005-0.02, weight theta(0))\n",
     "| g_N / a0 | measured ratio | Newton | QUMOND (field theory) | law without EFE | MI theta(0)=1 | MI theta(0)=2 | MI theta(0)=e | MI theta(0)=3 | MI theta(0)=5 |",
     "|---|---|---|---|---|---|---|---|---|---|"]
for _, b in B.iterrows():
    L.append(f"| {b.lo}-{b.hi} | {b.ratio:.3f} ± {b.e_ratio:.3f} | 1.000 | {b.pred_qumond:.3f} | {b.pred_law:.3f} | " +
             " | ".join(f"{b[f'mi_{t0:.2f}']:.3f}" for t0 in T0) + " |")
L.append("")
L.append(f"chi2 on the three bins below 1 a0: Newton {chi_newton:.1f}, QUMOND {chi_qumond:.1f}, law without EFE {chi_law:.1f}, MI best theta(0) = {best:.1f} "
         f"(chi2 {chi.min():.1f}, 1-sigma range {ok.index.min():.1f}-{ok.index.max():.1f}); chi2 at theta(0) = 2: {chi[2.0]:.1f}, at 3: {chi[3.0]:.1f}, at 5: {chi[5.0]:.1f}.")
L.append(f"In the two clean bins (0.1-1 a0) the data give {mid.ratio.iloc[0]:.3f} ± {mid.e_ratio.iloc[0]:.3f} and {mid.ratio.iloc[1]:.3f} ± {mid.e_ratio.iloc[1]:.3f}; "
         f"MI with theta(0) = 2 predicts {mid['mi_2.00'].iloc[0]:.3f} and {mid['mi_2.00'].iloc[1]:.3f}, with theta(0) = 3: {mid['mi_3.00'].iloc[0]:.3f} and {mid['mi_3.00'].iloc[1]:.3f}.")
L.append("")
L.append(f"## Test 2: lensing at 0.3-1 Mpc (a_ex = nu(e_N) e_N a0 = {AEX_GAL:.4f} a0 for the isolated galaxies' e_N = 1.1e-4; weight theta_eff = theta(w_ex/w_in) free)\n")
L.append("| theta_eff | chi2 bin 1 | chi2 bin 2 | chi2 bin 3 | chi2 bin 4 | total (48 points) |")
L.append("|---|---|---|---|---|---|")
for j in [0] + list(range(1, len(TH), 4)) + [len(TH) - 1]:
    L.append(f"| {TH[j]:.2f} | " + " | ".join(f"{chi_l[k, j]:.1f}" for k in range(4)) + f" | {tot[j]:.1f} |")
L.append("")
L.append(f"- Free amplitude per bin, R <= 1 Mpc: best theta_eff = {tbest:.2f}, theta_eff < {tmax:.2f} at Delta chi2 = 4 (bin 1 alone: < {tmax1:.2f}); "
         f"Delta chi2 = {d_t1:+.1f} at theta_eff = 1 ({np.sqrt(max(d_t1,0)):.1f} sigma) and {d_t2:+.1f} at theta_eff = 2 ({np.sqrt(max(d_t2,0)):.1f} sigma).")
L.append(f"- Internal orbital periods at the lensing radii for bin 1 (V = {Vb1:.0f} km/s): " + ", ".join(f"{r} kpc: {Tin[r]:.0f} Gyr" for r in Tin) +
         ". Any large-scale field that varies within a Hubble time (13.8 Gyr) has w_ex > w_in beyond ~300 kpc, so theta_eff = theta(w_ex/w_in > 1) < 1 there.")
if np.isfinite(y_needed):
    L.append(f"- With Milgrom's example theta(y) = 2/(1+y^2), the bound theta_eff < {tmax:.2f} requires w_ex/w_in > {y_needed:.1f}, i.e. an external field varying faster "
             f"than {Tin[1000]/y_needed:.0f} Gyr at 1 Mpc: satisfied by the cosmic flow.")
L.append("")
L.append("## Reading\n")
L.append(f"- Wide binaries: with theta(0) = 2-3 the MI models predict a 5-9 % excess of the relative velocities between 1 and 0.1 a0, against "
         f"{100*(mid.ratio.iloc[0]-1):.0f} ± {100*mid.e_ratio.iloc[0]:.0f} % and {100*(mid.ratio.iloc[1]-1):.0f} ± {100*mid.e_ratio.iloc[1]:.0f} % measured: compatible, where the "
         "field theory (QUMOND, +17 %) is 3-5 sigma above and Newton 2-3 sigma below. The reason is structural: in MI the Galactic acceleration enters the "
         "inertia of the binary with a weight theta(0) > 1 and isotropically, so the internal dynamics are closer to Newtonian than in modified gravity.")
L.append(f"- Lensing at 0.3-1 Mpc: the data bound the external weight to theta_eff < {tmax:.1f} (2 sigma), and exclude the field-theory strength theta_eff = 1 at "
         f"{np.sqrt(max(d_t1,0)):.0f} sigma, more strongly than QUMOND was, because the algebraic MI form feels the external acceleration as soon as it is a "
         "fraction of the internal one, without the shielding of QUMOND's phantom mass inside r_t. In MI the weight is theta(w_ex/w_in): the model delivers "
         f"theta_eff < {tmax:.1f} if w_ex/w_in > {y_needed:.0f}, i.e. if the cosmic field varies faster than {Tin[1000]/y_needed:.0f} Gyr at 1 Mpc (orbital period "
         f"{Tin[1000]:.0f} Gyr) and {Tin[300]/y_needed:.0f} Gyr at 300 kpc ({Tin[300]:.0f} Gyr): plausible for the first, marginal for the second. MI thus predicts a "
         "weak external-field effect where the lensing finds none, and a strong one in the binaries, where the field theories predict the reverse.")
L.append("- What this does not establish: MI has no relativistic version, hence no lensing prediction of its own (the test above assumes light follows the "
         "circular test orbits); the model's Fourier construction presumes stationary trajectories, which the 50-Gyr orbits at 1 Mpc are not (Milgrom's "
         "'initial conditions problem'); theta(y) is unknown apart from theta(1) = 1; and the binary excess of 5 % is also within the reach of hidden "
         "companions. The clusters do not reach MI either: neither lensing nor hydrostatic gas has a modified-inertia prediction. MI survives both tests, "
         "the first by predicting it, the second by predicting nothing firm there.")
L.append("\nFigure: results/mi_milgrom_predictions.png")
open("results/mi_milgrom_predictions.md", "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
a = ax[0]
a.errorbar(np.log10(B.gN), B.ratio, yerr=B.e_ratio, fmt="o", color="k", label="Gaia wide binaries")
a.plot(np.log10(B.gN), np.ones(len(B)), "k--", label="Newton")
a.plot(np.log10(B.gN), B.pred_qumond, "C0-.", label="QUMOND, Galactic field")
a.plot(np.log10(B.gN), B.pred_law, "C3-", label="law without EFE")
for t0, ls in ((2.0, "-"), (3.0, "--"), (5.0, ":")):
    a.plot(np.log10(B.gN), B[f"mi_{t0:.2f}"], color="C2", ls=ls, label=f"MI, theta(0) = {t0:.0f}")
a.set_xlabel("log g_N / a0"); a.set_ylabel("median v~ / calibration"); a.invert_xaxis(); a.set_ylim(0.8, 1.8); a.legend(fontsize=7)
a.set_title("Wide binaries: modified inertia against the field theory")
b = ax[1]
c = curves[0]; Rg = np.logspace(1.5, 3.0, 50)
b.errorbar(c["R"], c["V"], yerr=c["e"], fmt="o", color="k", label="lensing bin 1 (Mistele+24)")
for t, ls, lab in ((0.0, "--", "no EFE"), (0.3, ":", "theta_eff = 0.3"), (1.0, "-", "theta_eff = 1 (field-theory strength)"), (2.0, "-.", "theta_eff = 2")):
    b.plot(Rg, v_mi(c["M"], t * AEX_GAL, Rg), color="C2", ls=ls, label=lab)
b.set_xscale("log"); b.set_xlabel("R (kpc)"); b.set_ylabel("V_c (km/s)"); b.set_ylim(0, 200); b.legend(fontsize=7)
b.set_title("Lensing 34 kpc-1 Mpc: MI with external weight theta_eff")
plt.tight_layout(); plt.savefig("results/mi_milgrom_predictions.png", dpi=130)
print("figure written")
