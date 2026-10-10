# SPARC Energy-Mass Model

Radial integration model exploring energy accumulation as a contributor to apparent dynamical mass in galaxies.

## Installation
```bash
git clone https://github.com/djoubertthot/sparc-energy-model.git
cd sparc-energy-model
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Data
```bash
python scripts/download_sparc.py
```
Place rotmod files in `data/rotmod/` if not auto-downloaded.

## Running
```bash
python scripts/generate_rms_csv.py
python scripts/analyze_results.py
```

## Documentation
- THEORY.md : theory, parameters, methodology
- MEMORY.md : development history

## Contributing
This is an exploratory research repository. Contributions on age estimates or alternative amplification laws are welcome.
