import urllib.request
from pathlib import Path

def download_sparc():
    """Download SPARC data."""
    Path("data").mkdir(exist_ok=True)
    Path("data/rotmod").mkdir(exist_ok=True)
    
    urls = {
        "table": "http://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt",
        "rotmod_zip": "http://astroweb.case.edu/SPARC/Rotmod_LTG.zip"  # Update if needed
    }
    
    for name, url in urls.items():
        try:
            dest = f"data/{name}.dat" if name == "table" else f"data/{name}.zip"
            urllib.request.urlretrieve(url, dest)
            print(f"Downloaded {name} to {dest}")
        except Exception as e:
            print(f"Failed to download {name}: {e}")
            print("Manual download recommended from http://astroweb.case.edu/SPARC/")
    
    print("\n=== Other catalogs ===")
    print("CALIFA/PyCASSO: http://pycasso.ufsc.br/ (manual, large)")
    print("MaNGA FIREFLY: https://www.sdss.org/dr17/manga/manga-data/manga-firefly-value-added-catalog/ (manual, large)")
    print("Place any extracted age tables in data/ for later cross-matching.")

if __name__ == "__main__":
    download_sparc()
