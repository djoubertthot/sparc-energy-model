"""Build data/sparc_sample/sparc_enriched.csv from the downloaded catalogs.

SPARC Table 1 has no coordinates, no dust mass, no metallicity and no
stellar mass. Vflat = 0 means the flat velocity was not measured. Hubble-flow
distances (f_D = 1) use a coarse percentage error.

Adopted replacements:
- coordinates from SIMBAD, or from DustPedia when the name matches
- stellar mass from Marasco+2025 SED fits, else the Duey+2025 W1 color M/L
  applied to the SPARC 3.6 um luminosity, else z0MGS
- dust mass and SFR from DustPedia THEMIS (Nersesian+2019), else z0MGS SFR
- gas-phase 12+log(O/H) from De Vis+2019 (Pilyugin & Grebel 2016 S scale)
- Vflat from Lelli+2019 when SPARC left it at 0
- a Hubble-flow distance only when z0MGS or Marasco+2025 is tighter and
  excludes the SPARC value; a tighter error bar is kept even when the
  value stays

Rotation curves are rewritten from Rotmod_LTG.zip. Radii scale with distance
and the disk, bulge and gas velocities scale with the square root, which is
the SPARC mass-model scaling. Unchanged galaxies are copied byte for byte.
"""
import math
import re
import subprocess
import tempfile
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "data" / "catalogs"
ROTMOD = ROOT / "data" / "rotmod"
OUT = ROOT / "data" / "sparc_sample"
ZIP = CATALOG / "Rotmod_LTG.zip"
SIMBAD_CACHE = CATALOG / "simbad_coords.csv"

SPARC_ROW = re.compile(
    r"^\s*(?P<Galaxy>[A-Za-z0-9.+-]+)\s+"
    r"(?P<T>\d+)\s+"
    r"(?P<D>\d+(?:\.\d+)?)\s+"
    r"(?P<e_D>\d+(?:\.\d+)?)\s+"
    r"(?P<f_D>\d+)\s+"
    r"(?P<Inc>\d+(?:\.\d+)?)\s+"
    r"(?P<e_Inc>\d+(?:\.\d+)?)\s+"
    r"(?P<L36>\d+(?:\.\d+)?)\s+"
    r"(?P<e_L36>\d+(?:\.\d+)?)\s+"
    r"(?P<Reff>\d+(?:\.\d+)?)\s+"
    r"(?P<SBeff>\d+(?:\.\d+)?)\s+"
    r"(?P<Rdisk>\d+(?:\.\d+)?)\s+"
    r"(?P<SBdisk>\d+(?:\.\d+)?)\s+"
    r"(?P<MHI>\d+(?:\.\d+)?)\s+"
    r"(?P<RHI>\d+(?:\.\d+)?)\s+"
    r"(?P<Vflat>\d+(?:\.\d+)?)\s+"
    r"(?P<e_Vflat>\d+(?:\.\d+)?)\s+"
    r"(?P<Q>\d+)\s*"
    r"(?P<Ref>\S*)\s*$"
)
LN10 = math.log(10)


def field(line, start, end):
    if len(line) < start:
        return ""
    return line[start - 1 : end].strip()


def fnum(text):
    if text is None:
        return None
    text = str(text).strip()
    if text in {"", "-", "-99", "-99.0", "-99.00", "-9", "-9.0", "-9.00"}:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def name_keys(name):
    raw = re.sub(r"[\s_]+", "", str(name).upper()).replace("–", "-")
    keys = {raw}
    keys.add(re.sub(r"(?<=[A-Z])0+(?=\d)", "", raw))
    match = re.fullmatch(r"(NGC|IC|UGC|UGCA|DDO|PGC)0*(\d+)([A-Z]?)", raw)
    if match:
        keys.add(f"{match.group(1)}{int(match.group(2))}{match.group(3)}")
    match = re.fullmatch(r"ESO0*(\d+)-G?0*(\d+)", raw)
    if match:
        left, right = int(match.group(1)), int(match.group(2))
        keys.update({f"ESO{left}-{right}", f"ESO{left}-G{right}"})
    return {key for key in keys if key}


def read_text(path):
    return path.read_text(encoding="latin1")


def load_sparc():
    rows = []
    for line in read_text(CATALOG / "SPARC_Lelli2016c.mrt").splitlines():
        match = SPARC_ROW.match(line)
        if not match:
            continue
        row = match.groupdict()
        parsed = {"Galaxy": row["Galaxy"], "Ref": row["Ref"]}
        for key in row:
            if key in {"Galaxy", "Ref"}:
                continue
            parsed[key] = float(row[key]) if key not in {"T", "f_D", "Q"} else int(row[key])
        rows.append(parsed)
    if len(rows) != 175:
        raise SystemExit(f"Expected 175 SPARC galaxies, parsed {len(rows)}")
    return rows


def load_bulges():
    masses = {}
    for line in read_text(CATALOG / "Bulges.mrt").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        value = fnum(parts[-1])
        if value is None or not re.match(r"[A-Za-z]", parts[0]):
            continue
        masses[parts[0]] = value
    return masses


def load_wise():
    rows = {}
    for line in read_text(CATALOG / "wise_ii.table1").splitlines():
        galaxy = field(line, 1, 11)
        log_l = fnum(field(line, 12, 18))
        log_m = fnum(field(line, 19, 25))
        if not galaxy or log_l is None or log_m is None or log_m < 5:
            continue
        rows[galaxy] = {
            "logL": log_l,
            "logM": log_m,
            "e_logM": fnum(field(line, 26, 31)),
            "color": fnum(field(line, 32, 36)),
            "ML": fnum(field(line, 37, 41)),
        }
    return rows


def load_btfr():
    rows = {}
    for line in read_text(CATALOG / "BTFR_Lelli2019.mrt").splitlines():
        parts = line.split()
        if len(parts) < 12 or not re.match(r"[A-Za-z]", parts[0]):
            continue
        if fnum(parts[1]) is None or not 5 < fnum(parts[1]) < 13:
            continue
        vf = fnum(parts[5])
        rows[parts[0]] = {
            "Inc": fnum(parts[3]),
            "e_Inc": fnum(parts[4]),
            "Vf": None if vf in (None, 0.0) else vf,
            "e_Vf": fnum(parts[6]),
        }
    return rows


def load_marasco():
    rows = {}
    for line in read_text(CATALOG / "marasco2025_props.dat").splitlines():
        parts = line.split()
        if len(parts) < 14 or not re.match(r"[A-Za-z]", parts[0]):
            continue
        rows[parts[0]] = {
            "D": fnum(parts[1]),
            "e_D": fnum(parts[2]),
            "logMphot": fnum(parts[5]),
            "e_logMphot": fnum(parts[6]),
            "pval": fnum(parts[7]),
            "Vfl": fnum(parts[12]),
            "e_Vfl": fnum(parts[13]),
        }
    return rows


def index_rows(rows):
    index = {}
    ambiguous = set()
    for row in rows:
        for key in row["keys"]:
            if key in ambiguous:
                continue
            if key in index and index[key] is not row:
                ambiguous.add(key)
                index.pop(key, None)
                continue
            index[key] = row
    return index, ambiguous


def match_one(galaxy, index, ambiguous):
    found = []
    seen = set()
    for key in name_keys(galaxy):
        if key in ambiguous:
            return None, "ambiguous"
        row = index.get(key)
        if row is None or id(row) in seen:
            continue
        seen.add(id(row))
        found.append(row)
    if len(found) == 1:
        return found[0], "name"
    if len(found) > 1:
        return None, "ambiguous"
    return None, ""


def load_dustpedia():
    rows = []
    for line in read_text(CATALOG / "dustpedia_themis.dat").splitlines():
        name = field(line, 1, 23)
        if not name or name.lower() == "galaxy":
            continue
        rows.append({
            "name": name,
            "keys": name_keys(name),
            "RA": fnum(field(line, 25, 33)),
            "Dec": fnum(field(line, 35, 43)),
            "D": fnum(field(line, 50, 61)),
            "SFR": fnum(field(line, 63, 70)),
            "e_SFR": fnum(field(line, 72, 79)),
            "Mstar": fnum(field(line, 81, 89)),
            "Mdust": fnum(field(line, 205, 213)),
            "e_Mdust": fnum(field(line, 215, 223)),
            "chi2": fnum(field(line, 297, 302)),
        })
    return rows


def load_z0mgs():
    rows = []
    for line in read_text(CATALOG / "z0mgs_table4.dat").splitlines():
        keys = set()
        for start, end in ((9, 17), (19, 28), (30, 36)):
            keys.update(name_keys(field(line, start, end)))
        if not keys:
            continue
        rows.append({
            "keys": keys,
            "D": fnum(field(line, 38, 42)),
            # VizieR tags this column [Mpc], but the range 0.01-1.1 and the
            # ReadMe text are a log10 uncertainty. 0.04 dex on NGC 2403 is the
            # published floor; 0.04 Mpc would be a 1% error, which it is not.
            "e_dex": fnum(field(line, 44, 47)),
            "logM": fnum(field(line, 49, 53)),
            "e_logM": fnum(field(line, 55, 58)),
            "logSFR": fnum(field(line, 74, 78)),
            "e_logSFR": fnum(field(line, 80, 83)),
            "sfr_method": field(line, 85, 92),
        })
    return rows


def load_hi():
    rows = []
    for line in read_text(CATALOG / "devis2019_hi.dat").splitlines():
        name = field(line, 1, 23)
        if not name:
            continue
        rows.append({
            "keys": name_keys(name),
            "D": fnum(field(line, 45, 60)),
            "MHI": fnum(field(line, 126, 149)),
            "flag": field(line, 238, 238),
        })
    return rows


def load_oh():
    rows = []
    for line in read_text(CATALOG / "devis2019_metallicity.dat").splitlines():
        name = field(line, 1, 23)
        oh = fnum(field(line, 78, 94))
        if not name or oh is None or not 7.0 <= oh <= 9.6:
            continue
        rows.append({
            "keys": name_keys(name),
            "RA": fnum(field(line, 25, 33)),
            "Dec": fnum(field(line, 35, 43)),
            "OH": oh,
            "e_down": fnum(field(line, 96, 116)),
            "e_up": fnum(field(line, 118, 138)),
        })
    return rows


def simbad_aliases(name):
    aliases = []

    def add(alias):
        if alias and alias not in aliases:
            aliases.append(alias)

    match = re.fullmatch(r"(NGC|IC|UGC|UGCA|DDO|PGC)(\d+)([A-Za-z]?)", name)
    if match:
        prefix, number, suffix = match.group(1), int(match.group(2)), match.group(3)
        add(f"{prefix} {number}{suffix}")
        if prefix == "UGC":
            add(f"UGC {number:05d}")
        if prefix == "DDO":
            add(f"DDO {number:03d}")
        add(name)
    else:
        add(name)
    match = re.fullmatch(r"ESO(\d+)-G(\d+)", name)
    if match:
        left, right = int(match.group(1)), int(match.group(2))
        add(f"ESO {left:03d}-G{right:03d}")
        add(f"ESO {left:03d}-{right:03d}")
    if name.lower() == "camb":
        add("Camelopardalis B")
        add("Cam B")
    match = re.fullmatch(r"([DF])(\d+)-V(\d+)", name)
    if match:
        add(f"LSBC {name}")
        add(f"LSBC {match.group(1)}{match.group(2)}-V{int(match.group(3)):02d}")
    match = re.fullmatch(r"([DF])(\d+)-(\d+)", name)
    if match:
        add(f"LSBC {name}")
        add(f"LSBC {match.group(1)}{match.group(2)}-{int(match.group(3)):02d}")
    return aliases


def is_stellar(otype):
    # "Starburst" contains the letters star but is a galaxy class.
    if "*" in otype:
        return True
    return re.search(r"star(?!burst)", otype, re.IGNORECASE) is not None


def query_simbad(ident):
    script = (
        "output console=off script=off\n"
        'format object "%COO(d;A)|%COO(d;D)|%OTYPE(S)\\n"\n'
        f"query id {ident}\n"
    )
    # curl @- in parallel reads a shared stdin on Windows and swaps replies.
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="ascii") as handle:
        handle.write(script)
        path = handle.name
    try:
        proc = subprocess.run(
            ["curl.exe", "-fsS", "--max-time", "40", "-F", f"script=@{path}",
             "https://simbad.cds.unistra.fr/simbad/sim-script"],
            capture_output=True,
        )
    finally:
        Path(path).unlink(missing_ok=True)
    text = proc.stdout.decode("latin1", errors="replace")
    for line in text.splitlines():
        match = re.fullmatch(
            r"\s*([+-]?\d+(?:\.\d+)?)\|([+-]?\d+(?:\.\d+)?)\|(.*)\s*",
            line,
        )
        if not match:
            continue
        otype = match.group(3).strip()
        if is_stellar(otype):
            return None
        return float(match.group(1)), float(match.group(2)), ident, otype
    return None


def resolve_simbad(galaxy):
    for ident in simbad_aliases(galaxy):
        found = query_simbad(ident)
        if found:
            ra, dec, used, otype = found
            return {"Galaxy": galaxy, "RA": ra, "Dec": dec, "query": used, "otype": otype}
    return {"Galaxy": galaxy, "RA": None, "Dec": None, "query": "", "otype": ""}


def load_simbad(galaxies):
    cached = {}
    if SIMBAD_CACHE.exists():
        old = pd.read_csv(SIMBAD_CACHE)
        cached = {row["Galaxy"]: row for row in old.to_dict("records")}
    missing = [name for name in galaxies if name not in cached]
    if missing:
        print(f"SIMBAD lookup for {len(missing)} galaxies")
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(resolve_simbad, name): name for name in missing}
            done = 0
            for future in as_completed(futures):
                result = future.result()
                cached[result["Galaxy"]] = result
                done += 1
                if done % 25 == 0 or done == len(missing):
                    print(f"  {done}/{len(missing)}")
        pd.DataFrame(cached.values()).to_csv(SIMBAD_CACHE, index=False)
    return cached


def separation_arcsec(ra1, dec1, ra2, dec2):
    dra = (ra1 - ra2) * math.cos(math.radians((dec1 + dec2) / 2))
    return math.hypot(dra, dec1 - dec2) * 3600


def match_by_coordinate(ra, dec, rows, radius=20):
    if ra is None or dec is None:
        return None
    found = []
    for row in rows:
        if row.get("RA") is None or row.get("Dec") is None:
            continue
        if separation_arcsec(ra, dec, row["RA"], row["Dec"]) <= radius:
            found.append(row)
    if len(found) == 1:
        return found[0]
    return None


def frac_error(distance, error, dex=False):
    if distance is None or distance <= 0 or error is None or error <= 0:
        return None
    if dex:
        return error * LN10
    return error / distance


def consider_distance(old_d, old_e, new_d, new_frac, source):
    """Return a replacement for a Hubble-flow distance, or None.

    The central value moves only when the new measurement is tighter and
    excludes the SPARC value. A tighter error that still covers the SPARC
    value updates the uncertainty and leaves the curve distance unchanged.
    """
    if new_d is None or new_d <= 0 or new_frac is None or new_frac <= 0 or new_frac > 0.35:
        return None
    old_frac = old_e / old_d
    if new_frac > 0.8 * old_frac:
        return None
    relative = abs(new_d - old_d) / old_d
    if relative > 0.5:
        return "tension"
    if relative > max(0.05, new_frac):
        return {"kind": "value", "D": new_d, "e_D": new_d * new_frac, "source": source}
    return {"kind": "error", "D": old_d, "e_D": old_d * new_frac, "source": source}


def choose_distance(row, zrow, mrow):
    if row["f_D"] != 1:
        return None, []
    flags = []
    options = []
    if zrow and zrow["D"]:
        result = consider_distance(
            row["D"], row["e_D"], zrow["D"], frac_error(zrow["D"], zrow["e_dex"], dex=True),
            "z0MGS Leroy+2019",
        )
        if result == "tension":
            flags.append("distance_tension_z0mgs")
        elif result:
            options.append(result)
    if mrow and mrow["D"] and mrow["e_D"]:
        result = consider_distance(
            row["D"], row["e_D"], mrow["D"], mrow["e_D"] / mrow["D"],
            "Marasco+2025",
        )
        if result == "tension":
            flags.append("distance_tension_marasco")
        elif result:
            options.append(result)
    if len(options) >= 2 and abs(options[0]["D"] - options[1]["D"]) / row["D"] > 0.15:
        flags.append("distance_catalog_tension")
        options = [item for item in options if item["kind"] == "error"]
    if not options:
        return None, flags
    best = min(options, key=lambda item: item["e_D"] / item["D"])
    return best, flags


def log_mass_error(base_dex, d_cat, e_cat_dex, d_use, e_use):
    error = base_dex or 0.0
    if d_cat and d_use and abs(d_cat - d_use) / d_use > 0.01:
        sig_use = (e_use / d_use) / LN10 if e_use else 0.0
        error = math.hypot(error, 2 * math.hypot(sig_use, e_cat_dex or 0.0))
    return error or None


def choose_stellar_mass(l36, d_use, e_use, mrow, wrow, zrow):
    """Stellar mass on the adopted distance. L36 is in 1e9 Lsun at that distance."""
    luminosity = l36 * 1e9
    if luminosity <= 0:
        return None
    pval = 0.0 if not mrow or mrow["pval"] is None else mrow["pval"]
    # p = 0.000 means the SED was rejected, but it is still the only mass
    # Marasco published for a few dwarfs that have no WISE or z0MGS mass.
    rejected = pval < 0.001
    if mrow and mrow["logMphot"] and mrow["logMphot"] > 0 and (not rejected or wrow is None):
        scale = 0.0
        if mrow["D"] and abs(mrow["D"] - d_use) / d_use > 0.01:
            scale = 2 * math.log10(d_use / mrow["D"])
        mass = 10 ** (mrow["logMphot"] + scale)
        ratio = mass / luminosity
        if 0.15 <= ratio <= 2.0:
            dex = mrow["e_logMphot"] if mrow["e_logMphot"] not in (None, 0.0) else None
            e_cat = (mrow["e_D"] / mrow["D"]) / LN10 if mrow["D"] and mrow["e_D"] else 0.0
            return {
                "Mstar": mass,
                "e_logM": log_mass_error(dex or 0.0, mrow["D"], e_cat, d_use, e_use) if dex else None,
                "Y36": ratio,
                "source": "Marasco+2025 SED",
                "weak": (not rejected) and pval < 0.01,
                "rejected": rejected,
            }
    if wrow and wrow["ML"] and 0.15 <= wrow["ML"] <= 2.0:
        return {
            "Mstar": wrow["ML"] * luminosity,
            "e_logM": wrow["e_logM"],
            "Y36": wrow["ML"],
            "source": "Duey+2025 W1 M/L on L[3.6]",
            "weak": False,
            "rejected": False,
        }
    if zrow and zrow["logM"] and zrow["D"] and (zrow["e_logM"] or 1) <= 0.4:
        mass = 10 ** (zrow["logM"] + 2 * math.log10(d_use / zrow["D"]))
        ratio = mass / luminosity
        if 0.15 <= ratio <= 2.0:
            return {
                "Mstar": mass,
                "e_logM": log_mass_error(zrow["e_logM"], zrow["D"], zrow["e_dex"], d_use, e_use),
                "Y36": ratio,
                "source": "z0MGS Leroy+2019",
                "weak": False,
                "rejected": False,
            }
    return None


def scaled_dust(row, d_use):
    if row is None or row["Mdust"] is None or row["Mdust"] <= 0:
        return None
    if row["chi2"] is not None and row["chi2"] > 10:
        return None
    if not row["D"] or row["D"] <= 0:
        return None
    factor = (d_use / row["D"]) ** 2
    return row["Mdust"] * factor, (row["e_Mdust"] or 0) * factor


def choose_sfr(d_use, dust, zrow):
    if dust and dust["SFR"] and dust["SFR"] > 0 and dust["D"]:
        if dust["chi2"] is None or dust["chi2"] <= 10:
            factor = (d_use / dust["D"]) ** 2
            return dust["SFR"] * factor, (dust["e_SFR"] or None) and dust["e_SFR"] * factor, None, "DustPedia THEMIS Nersesian+2019"
    if zrow and zrow["logSFR"] is not None and zrow["D"]:
        log_sfr = zrow["logSFR"] + 2 * math.log10(d_use / zrow["D"])
        return 10 ** log_sfr, None, zrow["e_logSFR"], "z0MGS Leroy+2019"
    return None, None, None, ""


def curve_bytes():
    files = {}
    with zipfile.ZipFile(ZIP) as archive:
        for info in archive.infolist():
            if info.filename.lower().endswith("_rotmod.dat"):
                files[Path(info.filename).stem.replace("_rotmod", "")] = archive.read(info)
    return files


def scale_curve(raw, distance):
    text = raw.decode("ascii")
    match = re.search(r"Distance\s*=\s*([0-9.]+)", text)
    if not match:
        return raw, None
    old = float(match.group(1))
    factor = distance / old
    if abs(factor - 1) < 0.005:
        return raw, 1.0
    root = math.sqrt(factor)
    lines = []
    for line in text.splitlines(keepends=True):
        body = line.rstrip("\r\n")
        endl = line[len(body):]
        if body.startswith("#") or not body.strip():
            if "Distance" in body:
                body = re.sub(r"(Distance\s*=\s*)([0-9.]+)", rf"\g<1>{distance:.2f}", body, count=1)
            lines.append(body + endl)
            continue
        parts = body.split("\t")
        for index in (0, 3, 4, 5):
            if index >= len(parts) or not parts[index].strip():
                continue
            value = float(parts[index]) * (factor if index == 0 else root)
            parts[index] = f"{value:.2f}"
        lines.append("\t".join(parts) + endl)
    return "".join(lines).encode("ascii"), factor


def blank(value):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return ""
    return value


def build():
    sparc = load_sparc()
    bulges = load_bulges()
    wise = load_wise()
    btfr = load_btfr()
    marasco = load_marasco()
    dust_rows = load_dustpedia()
    z_rows = load_z0mgs()
    hi_rows = load_hi()
    oh_rows = load_oh()
    dust_index, dust_amb = index_rows(dust_rows)
    z_index, z_amb = index_rows(z_rows)
    hi_index, hi_amb = index_rows(hi_rows)
    oh_index, oh_amb = index_rows(oh_rows)
    curves = curve_bytes()
    simbad = load_simbad([row["Galaxy"] for row in sparc])

    enriched = []
    changes = []
    ROTMOD.mkdir(parents=True, exist_ok=True)

    for row in sparc:
        galaxy = row["Galaxy"]
        flags = []
        dust, dust_how = match_one(galaxy, dust_index, dust_amb)
        zrow, _ = match_one(galaxy, z_index, z_amb)
        hi, _ = match_one(galaxy, hi_index, hi_amb)
        oh, _ = match_one(galaxy, oh_index, oh_amb)
        mrow = marasco.get(galaxy)
        wrow = wise.get(galaxy)
        brow = btfr.get(galaxy)
        coord = simbad.get(galaxy, {})
        ra = coord.get("RA")
        dec = coord.get("Dec")
        if pd.isna(ra) or pd.isna(dec):
            ra, dec = None, None
        if ra is None and dust and dust.get("RA") is not None:
            ra, dec = dust["RA"], dust["Dec"]
            coord_source = "DustPedia"
        elif ra is not None:
            coord_source = "SIMBAD"
        else:
            coord_source = ""
            flags.append("coordinates_missing")
        if dust is None and ra is not None:
            dust = match_by_coordinate(ra, dec, dust_rows)
            dust_how = "coordinate" if dust else dust_how
        if oh is None and ra is not None:
            oh = match_by_coordinate(ra, dec, oh_rows)
        if dust_how == "ambiguous":
            flags.append("dustpedia_ambiguous")
        if dust and dust_how == "name" and ra is not None and dust.get("RA") is not None:
            offset = separation_arcsec(ra, dec, dust["RA"], dust["Dec"])
            if offset > 120:
                # A name match in DustPedia is the same galaxy. A SIMBAD
                # position this far away is a different object.
                flags.append("simbad_rejected")
                ra, dec, coord_source = dust["RA"], dust["Dec"], "DustPedia"
            elif offset > 60:
                flags.append("coordinate_tension_dustpedia")

        adopted, distance_flags = choose_distance(row, zrow, mrow)
        flags.extend(distance_flags)
        if adopted and adopted["kind"] == "value":
            d_use, e_use, d_source = adopted["D"], adopted["e_D"], adopted["source"]
            changes.append({
                "Galaxy": galaxy, "quantity": "distance",
                "previous": row["D"], "adopted": d_use, "unit": "Mpc",
                "source": d_source,
                "reason": "Hubble-flow distance excluded by a tighter catalog value",
            })
        elif adopted and adopted["kind"] == "error":
            d_use, e_use, d_source = row["D"], adopted["e_D"], adopted["source"] + " error"
            changes.append({
                "Galaxy": galaxy, "quantity": "e_distance",
                "previous": row["e_D"], "adopted": e_use, "unit": "Mpc",
                "source": adopted["source"],
                "reason": "tighter distance uncertainty; SPARC value kept",
            })
        else:
            d_use, e_use, d_source = row["D"], row["e_D"], "SPARC Lelli+2016"
        if row["f_D"] == 1 and e_use / d_use >= 0.2:
            flags.append("distance_uncertain")
        if dust and dust["D"] and abs(math.log10(dust["D"] / d_use)) > 0.15:
            flags.append("distance_tension_dustpedia")

        factor_table = d_use / row["D"]
        l36 = row["L36"] * factor_table ** 2
        e_l36 = row["e_L36"] * factor_table ** 2
        mhi = row["MHI"] * factor_table ** 2
        reff = row["Reff"] * factor_table
        rdisk = row["Rdisk"] * factor_table
        rhi = row["RHI"] * factor_table
        mhi_source = "SPARC Lelli+2016"

        inc, e_inc, inc_source = row["Inc"], row["e_Inc"], "SPARC Lelli+2016"
        if brow and brow["Inc"] and brow["e_Inc"] is not None:
            delta = abs(brow["Inc"] - row["Inc"])
            if delta > max(2.0, row["e_Inc"]) and brow["e_Inc"] <= row["e_Inc"]:
                inc, e_inc, inc_source = brow["Inc"], brow["e_Inc"], "Lelli+2019 BTFR"
                changes.append({
                    "Galaxy": galaxy, "quantity": "inclination",
                    "previous": row["Inc"], "adopted": inc, "unit": "deg",
                    "source": inc_source,
                    "reason": "later inclination excludes the SPARC value and has a smaller error",
                })
                flags.append("curve_not_recomputed_for_inclination")
        if inc < 40 or inc > 80 or e_inc >= 10:
            flags.append("deprojection_uncertain")

        vflat, e_vflat, v_source = row["Vflat"], row["e_Vflat"], "SPARC Lelli+2016"
        if row["Vflat"] == 0:
            filled = None
            if brow and brow["Vf"]:
                filled = brow["Vf"], brow["e_Vf"], "Lelli+2019 BTFR"
            elif mrow and mrow["Vfl"] and mrow["Vfl"] > 0:
                filled = mrow["Vfl"], mrow["e_Vfl"], "Marasco+2025"
            if filled:
                vflat, e_vflat, v_source = filled
                changes.append({
                    "Galaxy": galaxy, "quantity": "Vflat",
                    "previous": "", "adopted": vflat, "unit": "km/s",
                    "source": v_source, "reason": "SPARC sentinel 0 replaced with a measured flat velocity",
                })
            else:
                vflat, e_vflat, v_source = None, None, ""
                flags.append("Vflat_missing")
        if row["RHI"] == 0:
            rhi = None
            flags.append("RHI_missing")
        if row["Q"] == 3:
            flags.append("low_quality")

        if row["MHI"] == 0 and hi and hi["flag"] == "0" and hi["MHI"] and hi["D"]:
            mhi = hi["MHI"] / 1e9 * (d_use / hi["D"]) ** 2
            mhi_source = "De Vis+2019"
            changes.append({
                "Galaxy": galaxy, "quantity": "MHI",
                "previous": "", "adopted": mhi, "unit": "1e9 Msun",
                "source": mhi_source, "reason": "SPARC HI mass was 0",
            })
        elif hi and hi["flag"] == "0" and hi["MHI"] and hi["D"] and mhi:
            other = hi["MHI"] / 1e9 * (d_use / hi["D"]) ** 2
            if other > 0 and abs(math.log10(other / mhi)) > 0.3:
                flags.append("hi_mass_tension")

        stellar = choose_stellar_mass(l36, d_use, e_use, mrow, wrow, zrow)
        if stellar:
            changes.append({
                "Galaxy": galaxy, "quantity": "Mstar",
                "previous": 0.5 * l36 * 1e9, "adopted": stellar["Mstar"], "unit": "Msun",
                "source": stellar["source"],
                "reason": "constant M/L = 0.5 is imprecise; using a per-galaxy SED or color M/L",
            })
            if stellar["weak"]:
                flags.append("sed_fit_weak")
            if stellar["rejected"]:
                flags.append("sed_fit_rejected")
            if mrow and wrow and mrow["logMphot"] and wrow["logM"]:
                if abs((mrow["logMphot"]) - wrow["logM"]) > 0.3:
                    flags.append("mstar_tension")
        else:
            flags.append("Mstar_missing")

        dust_mass = scaled_dust(dust, d_use)
        if dust_mass:
            changes.append({
                "Galaxy": galaxy, "quantity": "Mdust",
                "previous": "", "adopted": dust_mass[0], "unit": "Msun",
                "source": "DustPedia THEMIS Nersesian+2019",
                "reason": "dust mass is absent from SPARC",
            })
        sfr, e_sfr, e_log_sfr, sfr_source = choose_sfr(d_use, dust, zrow)
        if sfr:
            changes.append({
                "Galaxy": galaxy, "quantity": "SFR",
                "previous": "", "adopted": sfr, "unit": "Msun/yr",
                "source": sfr_source, "reason": "SFR is absent from SPARC",
            })
        if oh:
            changes.append({
                "Galaxy": galaxy, "quantity": "12+logOH",
                "previous": "", "adopted": oh["OH"], "unit": "dex",
                "source": "De Vis+2019 PG16S",
                "reason": "gas-phase metallicity is absent from SPARC",
            })
        if ra is not None:
            changes.append({
                "Galaxy": galaxy, "quantity": "coordinates",
                "previous": "", "adopted": f"{ra:.5f},{dec:.5f}", "unit": "deg",
                "source": coord_source, "reason": "coordinates are absent from SPARC Table 1",
            })

        lbul = bulges.get(galaxy)
        raw = curves.get(galaxy)
        curve_factor = None
        if raw is None:
            flags.append("rotmod_missing")
        else:
            adjusted, curve_factor = scale_curve(raw, d_use)
            if curve_factor is None:
                flags.append("rotmod_distance_unreadable")
                adjusted = raw
            (ROTMOD / f"{galaxy}_rotmod.dat").write_bytes(adjusted)
            if curve_factor and abs(curve_factor - 1) >= 0.005:
                changes.append({
                    "Galaxy": galaxy, "quantity": "rotmod_distance_factor",
                    "previous": 1, "adopted": curve_factor, "unit": "",
                    "source": d_source,
                    "reason": "radii scale with D and Vgas, Vdisk, Vbul scale with sqrt(D)",
                })

        if lbul and lbul > 0.05 and raw is not None:
            text = raw.decode("ascii", errors="replace")
            vbul_max = 0.0
            for line in text.splitlines():
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    vbul_max = max(vbul_max, abs(float(parts[5])))
            if vbul_max == 0:
                flags.append("bulge_velocity_missing")

        enriched.append({
            "Galaxy": galaxy,
            "RA_deg": blank(ra),
            "Dec_deg": blank(dec),
            "coord_source": coord_source,
            "T": row["T"],
            "Q": row["Q"],
            "f_D": row["f_D"],
            "D_Mpc": d_use,
            "e_D_Mpc": e_use,
            "D_source": d_source,
            "Inc_deg": inc,
            "e_Inc_deg": e_inc,
            "Inc_source": inc_source,
            "Vflat_kms": blank(vflat),
            "e_Vflat_kms": blank(e_vflat),
            "Vflat_source": v_source,
            "L36_1e9_Lsun": l36,
            "e_L36_1e9_Lsun": e_l36,
            "MHI_1e9_Msun": mhi,
            "MHI_source": mhi_source,
            "Reff_kpc": blank(reff),
            "Rdisk_kpc": rdisk,
            "RHI_kpc": blank(rhi),
            "SBeff_Lsun_pc2": row["SBeff"],
            "SBdisk_Lsun_pc2": row["SBdisk"],
            "Lbul_1e9_Lsun": blank(lbul),
            "Mstar_Msun": blank(stellar["Mstar"] if stellar else None),
            "e_logMstar_dex": blank(stellar["e_logM"] if stellar else None),
            "Y36": blank(stellar["Y36"] if stellar else None),
            "Mstar_source": stellar["source"] if stellar else "",
            "Mdust_Msun": blank(dust_mass[0] if dust_mass else None),
            "e_Mdust_Msun": blank(dust_mass[1] if dust_mass else None),
            "Mdust_source": "DustPedia THEMIS Nersesian+2019" if dust_mass else "",
            "SFR_Msun_yr": blank(sfr),
            "e_SFR_Msun_yr": blank(e_sfr),
            "e_logSFR_dex": blank(e_log_sfr),
            "SFR_source": sfr_source,
            "OH": blank(oh["OH"] if oh else None),
            "e_OH_down": blank(oh["e_down"] if oh else None),
            "e_OH_up": blank(oh["e_up"] if oh else None),
            "OH_source": "De Vis+2019 PG16S" if oh else "",
            "sparc_D_Mpc": row["D"],
            "sparc_e_D_Mpc": row["e_D"],
            "sparc_Inc_deg": row["Inc"],
            "sparc_e_Inc_deg": row["e_Inc"],
            "sparc_Vflat_kms": row["Vflat"],
            "sparc_L36_1e9_Lsun": row["L36"],
            "sparc_MHI_1e9_Msun": row["MHI"],
            "sparc_Reff_kpc": row["Reff"],
            "sparc_Rdisk_kpc": row["Rdisk"],
            "sparc_RHI_kpc": row["RHI"],
            "curve_distance_factor": blank(curve_factor),
            "flags": ";".join(sorted(set(flags))),
            "Ref": row["Ref"],
        })

    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(enriched)
    log = pd.DataFrame(changes)
    frame.to_csv(OUT / "sparc_enriched.csv", index=False, float_format="%.6g")
    log.to_csv(OUT / "adjustment_log.csv", index=False, float_format="%.6g")
    return frame, log


def summarize(frame, log):
    print(f"galaxies {len(frame)}")
    print(f"coordinates {frame['RA_deg'].astype(str).ne('').sum()}")
    print("Mstar sources:")
    print(frame["Mstar_source"].replace("", "missing").value_counts().to_string())
    print(f"dust {frame['Mdust_source'].ne('').sum()}  sfr {frame['SFR_source'].ne('').sum()}  OH {frame['OH_source'].ne('').sum()}")
    moved = log[log["quantity"] == "distance"]
    tightened = log[log["quantity"] == "e_distance"]
    print(f"distance values replaced {len(moved)}")
    print(f"distance errors tightened {len(tightened)}")
    print(f"Vflat filled {(log['quantity'] == 'Vflat').sum()}")
    print(f"inclination replaced {(log['quantity'] == 'inclination').sum()}")
    print(f"curves rescaled {(log['quantity'] == 'rotmod_distance_factor').sum()}")
    flag_counts = {}
    for text in frame["flags"]:
        for flag in str(text).split(";"):
            if flag:
                flag_counts[flag] = flag_counts.get(flag, 0) + 1
    print("flags:")
    for flag, count in sorted(flag_counts.items(), key=lambda item: -item[1]):
        print(f"  {count:3} {flag}")
    ngc = frame.loc[frame["Galaxy"] == "NGC2403"].iloc[0]
    if abs(float(ngc["D_Mpc"]) - 3.16) > 1e-6:
        raise SystemExit("NGC2403 distance was not preserved")
    if abs(float(ngc["Vflat_kms"]) - 131.2) > 1e-6:
        raise SystemExit("NGC2403 Vflat was not preserved")
    if ngc["RA_deg"] != "" and abs(float(ngc["RA_deg"]) - 114.214) > 0.02:
        raise SystemExit(f"NGC2403 coordinates look wrong: {ngc['RA_deg']}")
    y36 = float(ngc["Y36"])
    if not 0.2 <= y36 <= 0.8:
        raise SystemExit(f"NGC2403 Y36 out of the expected SED range: {y36}")
    print(f"NGC2403 Y36 {y36:.3f} from {ngc['Mstar_source']}")
    placed = frame[frame["RA_deg"].astype(str).ne("")].copy()
    placed["RA_deg"] = placed["RA_deg"].astype(float)
    placed["Dec_deg"] = placed["Dec_deg"].astype(float)
    for i, left in placed.iterrows():
        for _, right in placed.loc[placed.index > i].iterrows():
            if separation_arcsec(left.RA_deg, left.Dec_deg, right.RA_deg, right.Dec_deg) < 5:
                raise SystemExit(f"coordinates collide: {left.Galaxy} and {right.Galaxy}")
    missing_coord = frame.loc[frame["RA_deg"].astype(str).eq(""), "Galaxy"].tolist()
    if missing_coord:
        print("no coordinates: " + ", ".join(missing_coord))


def main():
    required = [
        CATALOG / "SPARC_Lelli2016c.mrt",
        ZIP,
        CATALOG / "marasco2025_props.dat",
        CATALOG / "wise_ii.table1",
        CATALOG / "BTFR_Lelli2019.mrt",
        CATALOG / "dustpedia_themis.dat",
        CATALOG / "z0mgs_table4.dat",
        CATALOG / "devis2019_hi.dat",
        CATALOG / "devis2019_metallicity.dat",
        CATALOG / "Bulges.mrt",
    ]
    missing = [str(path.name) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing catalogs: " + ", ".join(missing) + ". Run scripts/download_catalogs.py")
    frame, log = build()
    summarize(frame, log)
    print(f"wrote {OUT / 'sparc_enriched.csv'}")
    print(f"wrote {OUT / 'adjustment_log.csv'}")


if __name__ == "__main__":
    main()
