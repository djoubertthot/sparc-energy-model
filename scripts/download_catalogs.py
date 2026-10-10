"""Download the SPARC tables and the catalogs used to fill missing values.

Raw files go to data/catalogs/ (gitignored). Rotation curves are extracted
to data/rotmod/. Re-run is safe: an existing file of the expected size is kept.
"""
import gzip
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "data" / "catalogs"
ROTMOD = ROOT / "data" / "rotmod"

SPARC = "https://astroweb.case.edu/SPARC"
CDS = "https://cdsarc.cds.unistra.fr/ftp"

FILES = {
    "SPARC_Lelli2016c.mrt": f"{SPARC}/SPARC_Lelli2016c.mrt",
    "MassModels_Lelli2016c.mrt": f"{SPARC}/MassModels_Lelli2016c.mrt",
    "Bulges.mrt": f"{SPARC}/Bulges.mrt",
    "wise_ii.table1": f"{SPARC}/wise_ii.table1",
    "BTFR_Lelli2019.mrt": f"{SPARC}/BTFR_Lelli2019.mrt",
    "Rotmod_LTG.zip": f"{SPARC}/Rotmod_LTG.zip",
    "sparc_table1.dat": f"{CDS}/J/AJ/152/157/table1.dat",
    "sparc_table2.dat": f"{CDS}/J/AJ/152/157/table2.dat",
    "sparc_refs.dat": f"{CDS}/J/AJ/152/157/refs.dat",
    "marasco2025_phot.dat": f"{CDS}/J/A+A/695/L23/tablec1.dat",
    "marasco2025_props.dat": f"{CDS}/J/A+A/695/L23/tablec2.dat",
    "dustpedia_themis.dat": f"{CDS}/J/A+A/624/A80/themis.dat",
    "z0mgs_table4.dat.gz": f"{CDS}/J/ApJS/244/24/table4.dat.gz",
    "devis2019_hi.dat": f"{CDS}/J/A+A/623/A5/dp-hi.dat",
    "devis2019_metallicity.dat": f"{CDS}/J/A+A/623/A5/dp-mg.dat",
}


def download(name, url):
    dest = CATALOG / name
    if dest.exists() and dest.stat().st_size > 0:
        print(f"keep    {name} ({dest.stat().st_size} bytes)")
        return dest
    print(f"get     {name}")
    data = None
    req = urllib.request.Request(url, headers={"User-Agent": "sparc-energy-model/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=180) as response:
            data = response.read()
    except Exception as exc:
        curl = shutil.which("curl") or shutil.which("curl.exe")
        if curl is None:
            print(f"FAILED  {name}: {exc}")
            return None
        print(f"retry   {name} via curl ({exc.__class__.__name__})")
        proc = subprocess.run(
            [curl, "-fsSL", "--retry", "2", "-o", str(dest), url],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0 or not dest.exists() or dest.stat().st_size < 100:
            print(f"FAILED  {name}: {proc.stderr.strip() or exc}")
            dest.unlink(missing_ok=True)
            return None
        print(f"saved   {name} ({dest.stat().st_size} bytes)")
        return dest
    if data is None or len(data) < 100:
        print(f"FAILED  {name}: response too small")
        return None
    dest.write_bytes(data)
    print(f"saved   {name} ({len(data)} bytes)")
    return dest


def inflate_gzip():
    for gz in CATALOG.glob("*.gz"):
        plain = gz.with_suffix("")
        if plain.suffix == ".dat" and plain.exists() and plain.stat().st_size > 0:
            continue
        if gz.suffix != ".gz":
            continue
        with gzip.open(gz, "rb") as src, plain.open("wb") as out:
            shutil.copyfileobj(src, out)
        print(f"unzip   {plain.name} ({plain.stat().st_size} bytes)")


def extract_rotmod(zip_path):
    if zip_path is None or not zip_path.exists():
        print("rotmod zip missing; curves will be built from the mass-model table")
        return 0
    ROTMOD.mkdir(parents=True, exist_ok=True)
    count = 0
    with zipfile.ZipFile(zip_path) as archive:
        for info in archive.infolist():
            if not info.filename.lower().endswith("_rotmod.dat"):
                continue
            target = ROTMOD / Path(info.filename).name
            target.write_bytes(archive.read(info))
            count += 1
    print(f"extracted {count} rotmod files to {ROTMOD}")
    return count


def main():
    CATALOG.mkdir(parents=True, exist_ok=True)
    ROTMOD.mkdir(parents=True, exist_ok=True)
    fetched = {}
    for name, url in FILES.items():
        fetched[name] = download(name, url)
    inflate_gzip()
    extract_rotmod(fetched.get("Rotmod_LTG.zip"))
    failed = [name for name, path in fetched.items() if path is None]
    print(f"\n{len(fetched) - len(failed)}/{len(fetched)} catalogs on disk")
    if failed:
        print("missing: " + ", ".join(failed))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
