# -*- coding: utf-8 -*-
"""Parse the circular-velocity table of Mistele et al. 2024 (ApJL 969, L3; arXiv:2406.09685, plots/Vc-table.tex) into
data/external/lensing/mistele2024_vc.csv: sample (All/ETG/LTG), log10 R (kpc), V_c per baryonic-mass bin with statistical and
systematic errors. Bins: log M_b = 9.0-10.5, 10.5-10.8, 10.8-11.1, 11.1-11.5 (mean 10.10, 10.66, 10.96, 11.29)."""
import os
import re
import sys

import pandas as pd

src = sys.argv[1]
tex = open(src, encoding="utf-8").read()
tex = re.sub(r"%[^\n]*", "", tex)
tex = tex.replace("\\startdata", "")
rows = []
num = r"[-+]?\d+\.?\d*"
for line in tex.split("\\\\"):          # LaTeX row terminator
    line = " ".join(line.split())        # cells span several physical lines
    if "&" not in line:
        continue
    cells = [c.strip() for c in line.split("&")]
    if len(cells) < 6:
        continue
    sample = cells[0].strip()
    if sample not in ("All", "ETG", "LTG"):
        continue
    lr = float(re.search(num, cells[1]).group(0))
    row = dict(sample=sample, logR=lr, R_kpc=10 ** lr)
    for k in range(4):
        nums = re.findall(num, cells[2 + k].replace("\\pm", " "))
        if len(nums) >= 3:
            row[f"Vc_bin{k+1}"], row[f"e_stat_bin{k+1}"], row[f"e_sys_bin{k+1}"] = float(nums[0]), float(nums[1]), float(nums[2])
        else:
            row[f"Vc_bin{k+1}"] = row[f"e_stat_bin{k+1}"] = row[f"e_sys_bin{k+1}"] = float("nan")
    rows.append(row)
T = pd.DataFrame(rows)
os.makedirs("data/external/lensing", exist_ok=True)
T.to_csv("data/external/lensing/mistele2024_vc.csv", index=False)
print(T.to_string(index=False))
