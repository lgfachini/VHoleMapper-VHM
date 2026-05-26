# VHoleMapper (VHM)

VHoleMapper is a Python toolkit for detecting, validating, and analyzing pi-holes on molecular electrostatic potential (ESP) surfaces.

It supports batch workflows where each molecular system has its own folder with an ESP surface file (`vtx.txt`) and an XYZ geometry.

## Features

- Automatic pi-hole candidate detection from local ESP maxima
- Local reference frame construction from four user-selected atoms
- Above/below plane searches with configurable axial and radial windows
- Optional Moving Least Squares (MLS) derivative validation
- Sector-based validation of pi-hole anisotropy
- CSV, TXT, XYZ, and SVG outputs

## Project Structure

```text
VHoleMapper-VHM/
  main.py                  # Backward-compatible runner
  pyproject.toml           # Packaging and test configuration
  requirements.txt
  config/
    settings.py            # User-editable run settings
    settings.example.py
  src/
    vhm/
      cli.py               # Command-line driver
      config.py            # Dataclasses and path resolution
      pipeline.py          # Workflow orchestration
      geometry.py          # Local frame and coordinate transforms
      io_utils.py          # XYZ and VTX parsing
      search.py            # Local maximum detection
      mls.py               # MLS fitting and derivative checks
      sector.py            # Sector validation
      outputs.py           # CSV/TXT/XYZ writers
      plotting.py          # SVG plots
      models.py            # Result/data models
  tests/
  data/
```

## Install

```bash
pip install -e .
```

For development:

```bash
pip install -e ".[dev]"
```

## Input Layout

Each system should live in its own folder:

```text
data/
  benzene/
    benzene.xyz
    vtx.txt
  f6benzene/
    f6benzene.xyz
    vtx.txt
```

If `AUTO_XYZ_FILENAME` is `None`, VHM expects the XYZ filename to match the folder name.

## Configure

Edit [config/settings.py](config/settings.py).

Automatic mode scans `DATA_DIR` and applies the same reference atoms to every folder:

```python
RUN_MODE = "auto"
DATA_DIR = Path("data")
AUTO_REFERENCE_ATOMS = (1, 3, 6, 10)
```

Manual mode lets each target use its own atoms and filenames:

```python
RUN_MODE = "manual"

MANUAL_TARGETS = (
    TargetSpec(
        folder=Path("data/benzene"),
        reference_atom_indices_1based=(1, 3, 6, 10),
        orientation_atom_index_1based=None,
    ),
)
```

Atom indices are 1-based.

## Run

From the project root:

```bash
python main.py
```

After installing the package:

```bash
vhm
```

## Outputs

For each system folder, VHM writes:

- `pi_hole_candidates.csv`: all analyzed candidates
- `pi_holes_validated.txt`: human-readable validated pi-hole report
- `molecule_with_piholes.xyz`: molecule plus dummy `X` atoms at validated pi-holes
- `pihole_plots/`: SVG sector plots

## Notes

- XYZ coordinates default to angstrom.
- VTX coordinates default to bohr and are converted to angstrom.
- ESP values are assumed to be in kcal/mol.
- `delta_thr` is the minimum accepted `V_hole,mean`.

## Citation

If you use this software, please cite:

*Quantitative Analysis of Metal-Centered pi-Holes in {TM(cyclen)}2+ Complexes*

Journal: *ACS Inorganic Chemistry*

Year: 2026

DOI: 10.1021/acs.inorgchem.6c01328

## License

This project is licensed under the GPL-3.0 License.
