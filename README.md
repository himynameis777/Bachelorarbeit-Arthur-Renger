README — BAC/code

Purpose
- Short guide to set up and run the baseline calibration (`simulation.py`) and analysis (`analysis.r`) in this folder.

Files
- `abm.nlogox` — NetLogo model (expected model file).
- `simulation.py` — Python script to run baseline calibration experiments.
- `analysis.r` — R script that reads `run_summary.csv` and `population_timeseries.csv` and writes plots/outputs to `plots/` and `outputs/` inside this folder.

Prerequisites
- Python 3.8+ with `pip3`.
- R (3.6+ or newer) with ability to install CRAN packages.
- NetLogo installed (a matching NetLogo version; update `--netlogo-home` or `NETLOGO_HOME` if necessary).
- On macOS: install command-line developer tools if not present:

```bash
xcode-select --install
# or if Xcode is installed
# sudo xcode-select --switch /Applications/Xcode.app
```

Python dependencies
```bash
pip3 install --user pandas pynetlogo
```

R dependencies
```bash
Rscript -e 'install.packages(c("tidyverse","lubridate","patchwork","viridis","broom","scales"))'
```

Environment variables / options
- `NETLOGO_HOME`: Path to NetLogo installation (optional). Example:

```bash
export NETLOGO_HOME=/Applications/NetLogo\ 7.0.3
```

Running
- Quick smoke test (Python):

```bash
python3 simulation.py --quick
```

- Full run (Python):

```bash
python3 simulation.py --repetitions 3 --max-ticks 2000
# or pass --netlogo-home /path/to/NetLogo
```

- Analysis (R):

```bash
# From the `code/` folder
Rscript analysis.r
# or from repo root
Rscript code/analysis.r
```

Outputs
- CSV outputs are written to `code/outputs`.
- Plots are written to `code/plots`.

Troubleshooting
- If the Python run fails early with an `xcode-select` message, the macOS developer tools installation is required (see commands above).
- If `pandas` or `pynetlogo` import fails, install the Python dependencies shown above.
- If R complains about missing packages, run the `Rscript -e 'install.packages(...)'` command above.
- If the NetLogo model is not found, ensure `abm.nlogox` is present in `code/` or pass `--model /path/to/model.nlogox`.

Contact
- If you want me to change defaults (model filename, output paths, or add a `--check-env` flag), tell me which change you prefer.