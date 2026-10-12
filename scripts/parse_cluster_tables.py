# -*- coding: utf-8 -*-
"""Build the cluster tables used by scripts/cluster_rar.py:
  data/external/clusters/tian2020_clash_rar.csv : Tian et al. 2020 (ApJ 896, 70; VizieR J/ApJ/896/70/fig2) log g_bar, log g_tot per
      cluster and radius (14-600 kpc), 20 CLASH clusters;
  data/external/clusters/xcop_profiles.csv : X-COP, 12 clusters: total NFW (M200, c200; Eckert et al. 2022, A&A 662, A123),
      hydrostatic masses and gas fractions at R500 and R200 (Eckert et al. 2019, A&A 621, A40), gas and stellar fractions at
      R2500/R500/R200 where published (Eckert et al. 2022, five clusters).
Values are parsed from the ar5iv HTML of the two papers (data/external/clusters/eckert2019.html, eckert2022.html; strings like
'1214-4+41214_{-4}^{+4}' are read as the first number, the error from the LaTeX part)."""
import os
import re
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
D = "data/external/clusters"
S = D


def first_num(s):
    m = re.search(r"[-+]?\d+\.?\d*", str(s).replace("−", "-"))
    return float(m.group(0)) if m else np.nan


def err_num(s):
    """symmetric error from 'a±b' or 'a-b+c' strings"""
    s = str(s).replace("−", "-")
    m = re.search(r"\\pm\s*(\d+(?:\.\d+)?)", s)      # LaTeX part of the cell: 'a±ba\pm b'
    if m:
        return float(m.group(1))
    m = re.search(r"_\{-([\d.]+)\}\^\{\+([\d.]+)\}", s)
    if m:
        return 0.5 * (float(m.group(1)) + float(m.group(2)))
    return np.nan


# ---- Tian 2020
def tsv(p):
    lines = [l for l in open(p, encoding="utf-8", errors="replace") if not l.startswith("#") and l.strip()]
    hdr = lines[0].rstrip("\n").split("\t"); rows = [l.rstrip("\n").split("\t") for l in lines[3:]]
    df = pd.DataFrame(rows, columns=hdr)
    for c in df.columns:
        try:
            df[c] = pd.to_numeric(df[c].str.strip())
        except Exception:
            df[c] = df[c].str.strip()
    return df


f2 = tsv(f"{D}/tian2020_clash_rar_fig2.tsv")
f2 = f2.rename(columns={"AName": "cluster", "Rad": "R_kpc", "log(gbar)": "log_gbar", "log(gtot)": "log_gtot",
                        "e_log(gbar)": "e_log_gbar", "e_log(gtot)": "e_log_gtot"})
f2[["cluster", "R_kpc", "log_gbar", "e_log_gbar", "log_gtot", "e_log_gtot"]].to_csv(f"{D}/tian2020_clash_rar.csv", index=False)
print("CLASH:", len(f2), "points,", f2.cluster.nunique(), "clusters")

# ---- X-COP
t22 = pd.read_html(f"{S}/eckert2022.html")
nfw = t22[16].copy(); nfw.columns = ["cluster", "z", "NH", "R500", "M500", "M200", "c200", "M200_DM", "c200_DM"]
nfw = nfw[nfw.cluster.notna() & (nfw.cluster.astype(str) != "Cluster") & (nfw.cluster.astype(str) != "nan")]
nfw = nfw[nfw.z.apply(lambda v: str(v) not in ("nan", "NaN") and re.search(r"\d", str(v)) is not None)]
X = pd.DataFrame(dict(cluster=nfw.cluster.astype(str), z=nfw.z.apply(first_num), R500_kpc=nfw.R500.apply(first_num),
                      M500_1e14=nfw.M500.apply(first_num), M200_1e14=nfw.M200.apply(first_num), e_M200=nfw.M200.apply(err_num),
                      c200=nfw.c200.apply(first_num), e_c200=nfw.c200.apply(err_num)))
stars = t22[19].copy(); stars.columns = ["cluster", "fgas2500", "fstar2500", "fbar2500", "fgas500", "fstar500", "fbar500", "fgas200", "fstar200", "fbar200"]
stars["cluster"] = stars.cluster.astype(str)
for c in stars.columns[1:]:
    stars[c] = stars[c].apply(first_num) / 100.0
t19 = pd.read_html(f"{S}/eckert2019.html")[10].copy()
t19.columns = ["cluster", "MHSE500", "MHSE200", "fgas500", "fgas200", "alpha500", "alpha200", "Mtot500", "Mtot200"]
t19 = t19.iloc[2:]
G = pd.DataFrame(dict(cluster=t19.cluster.astype(str), fgas500=t19.fgas500.apply(first_num), e_fgas500=t19.fgas500.apply(err_num),
                      fgas200=t19.fgas200.apply(first_num), e_fgas200=t19.fgas200.apply(err_num),
                      MHSE500=t19.MHSE500.apply(first_num), MHSE200=t19.MHSE200.apply(first_num)))
X = X.merge(G, on="cluster", how="left").merge(stars[["cluster", "fgas2500", "fstar2500", "fstar500", "fstar200"]], on="cluster", how="left")
X.to_csv(f"{D}/xcop_profiles.csv", index=False)
print(X.to_string(index=False))
