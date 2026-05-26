# 🧲 VHoleMapper (VHM): A Pi-Hole Detection and Analysis Toolkit

VHoleMapper (VHM) is a Python-based toolkit for detecting, validating, and visualizing pi-holes on molecular electrostatic potential (ESP) surfaces.

It is designed for batch workflows where each molecular system is stored in its own folder with:

* an ESP surface file, usually `vtx.txt`
* an XYZ molecular geometry, usually named after the folder

VHM builds a local molecular reference frame, searches for ESP maxima above and below the selected plane, optionally applies a Moving Least Squares (MLS) derivative check, and validates candidates by comparing the candidate potential with surrounding angular sector averages.

---

## 📁 Structure

```text
VHoleMapper-VHM/
|-- config/
|   |-- settings.py          # edit your run configuration here
|   `-- settings.example.py
|-- data/                    # input systems and generated outputs
|   |-- benzene/
|   |   |-- benzene.xyz
|   |   `-- vtx.txt
|   `-- f6benzene/
|       |-- f6benzene.xyz
|       `-- vtx.txt
|-- src/vhm/
|   |-- cli.py               # command-line driver
|   |-- config.py            # dataclasses, validation, path resolution
|   |-- models.py            # LocalFrame, MlsFit, CandidateResult
|   |-- pipeline.py          # full analysis workflow
|   |-- geometry.py          # local frame construction and transforms
|   |-- io_utils.py          # XYZ and VTX parsers
|   |-- search.py            # axial filtering and local maxima search
|   |-- mls.py               # Moving Least Squares derivative analysis
|   |-- sector.py            # angular sector validation
|   |-- outputs.py           # CSV, TXT, and XYZ exporters
|   `-- plotting.py          # SVG pi-hole plots
|-- tests/
|-- main.py                  # shortcut entry: python main.py
|-- pyproject.toml
`-- requirements.txt
```

---

## ⚙️ Features

* Reads molecular geometries from `.xyz` files.
* Reads ESP surface points and potentials from `vtx.txt`.
* Defines a local coordinate frame from four user-selected reference atoms.
* Searches both sides of the reference plane for discrete local ESP maxima.
* Supports automatic batch mode and manual per-system target definitions.
* Optionally validates candidates with MLS gradient/Hessian criteria.
* Uses angular sector sampling to estimate:
  * `V_s,max`
  * `V_hole,mean`
  * `V_hole,max`
  * `V_hole,min`
  * `Delta V_hole`
  * anisotropy angle
* Exports CSV, TXT, XYZ, and SVG outputs.
* Runs from `config/settings.py` via `python main.py`, `python -m vhm`, or `vhm`.

---

## ▶️ How to Use

### 📦 1. Install Requirements

```bash
pip install -r requirements.txt
```

Or install the package in editable mode:

```bash
pip install -e .
```

For development and tests:

```bash
pip install -e ".[dev]"
```

### 📂 2. Add Your Input Data

Place each molecular system in its own folder under `data/`:

```text
data/
`-- molecule_name/
    |-- molecule_name.xyz
    `-- vtx.txt
```

If `AUTO_XYZ_FILENAME = None`, VHM expects the XYZ file to be named after the folder:

```text
data/benzene/benzene.xyz
data/f6benzene/f6benzene.xyz
```

If all folders use the same XYZ filename, set `AUTO_XYZ_FILENAME` in `config/settings.py`.

### 🚀 3. Configure and Run

Edit [config/settings.py](config/settings.py) at the project root.

Then run:

```bash
python main.py
```

Alternatives:

```bash
python -m vhm
vhm   # after pip install -e .
```

---

## 🧭 Run Modes

### 🔹 Automatic Mode

Automatic mode scans `DATA_DIR` and applies the same reference atoms to every selected folder.

```python
RUN_MODE = "auto"

DATA_DIR = Path("data")
AUTO_FOLDER_GLOB = "*"
AUTO_REFERENCE_ATOMS = (1, 3, 6, 10)
AUTO_ORIENTATION_ATOM = None
AUTO_XYZ_FILENAME = None
AUTO_VTX_FILENAME = "vtx.txt"
```

Useful when all systems share the same atom ordering and topology.

### 🔹 Manual Mode

Manual mode analyzes only the listed targets. Each target can have its own reference atoms, orientation atom, XYZ filename, and VTX filename.

```python
RUN_MODE = "manual"

MANUAL_TARGETS = (
    TargetSpec(
        folder=Path("data/benzene"),
        reference_atom_indices_1based=(1, 3, 6, 10),
        orientation_atom_index_1based=None,
        xyz_filename=None,
        vtx_filename="vtx.txt",
    ),
)
```

Atom indices are 1-based, meaning the first atom in the XYZ file is atom `1`.

---

## 🛠️ Main Configuration Options

All numerical analysis parameters are stored in the `ANALYSIS_OPTIONS` dictionary in [config/settings.py](config/settings.py).

### 📏 Coordinate Units

```python
"xyz_coord_unit": "angstrom",
"vtx_coord_unit": "bohr",
```

Accepted units are `"angstrom"` and `"bohr"`. Internally, VHM converts coordinates to angstrom.

### 🎯 Axial Search Region

```python
"max_distance_from_axis": 0.25,
"z_search_min": 0.8,
"z_search_max": 3.0,
```

These parameters define the cylindrical region around the local z axis where candidate pi-holes are searched.

### 🔍 Local Maximum Search

```python
"local_max_radius": 0.5,
"strict_local_max": True,
"min_local_neighbors": 8,
"require_positive_potential": False,
"top_n_candidates_per_side": None,
```

The code builds a local neighborhood around each surface point and keeps points whose ESP values are local maxima.

### 🧮 MLS Validation

```python
"do_mls_test": False,
"mls_radius": 1.50,
"mls_min_neighbors": 10,
"mls_sigma": 1.20,
"mls_include_candidate": False,
```

When enabled, the MLS test fits a local quadratic model around a candidate and evaluates tangential gradient and Hessian behavior.

Additional MLS criteria:

```python
"tangential_grad_norm_thr": 200.0,
"tangential_hess_eig_max_thr": 0.0,
"allow_semidefinite_tangential_maximum": True,
"use_tangential_stationary_shift_test": False,
"max_tangential_stationary_shift": 0.50,
```

### 🧩 Sector Validation

```python
"angular_sectors": 8,
"sector_inner_radius": 0.25,
"sector_outer_radius": 3.50,
"sector_half_height": 0.50,
"min_points_per_sector": 50,
"delta_thr": -100.0,
```

For each candidate, VHM samples surrounding ESP points in angular sectors around the candidate. A candidate is validated when:

```text
V_hole,mean >= delta_thr
```

### 📈 Plotting

```python
"plot_figsize": (8, 8),
"plot_atom_labels": True,
"plot_sector_ring_radius": None,
"plot_show_all_shell_points": False,
"plot_only_validated": True,
```

SVG plots show the reference atoms projected into the local plane, the candidate, angular sectors, and anisotropy markers.

---

## 🔬 Workflow

For each target folder, VHM performs:

1. Resolve input and output paths.
2. Read the XYZ geometry and VTX surface data.
3. Convert coordinates to angstrom.
4. Build the local frame from four reference atoms.
5. Transform surface points into local coordinates.
6. Search above and below the reference plane for local ESP maxima.
7. Optionally apply MLS derivative validation.
8. Run angular sector validation.
9. Write CSV, TXT, XYZ, and SVG outputs.
10. Print a terminal summary.

---

## 📊 Outputs

For each analyzed system, VHM writes outputs next to the input files:

```text
data/system_name/
|-- pi_hole_candidates.csv
|-- pi_holes_validated.txt
|-- molecule_with_piholes.xyz
`-- pihole_plots/
    `-- pihole_01_above_idx1234.svg
```

### `pi_hole_candidates.csv`

Contains all analyzed candidates, including rejected and validated candidates, with local/global coordinates, ESP values, MLS descriptors, sector descriptors, validation status, and rejection/validation reason.

### `pi_holes_validated.txt`

Human-readable report containing only validated pi-holes.

### `molecule_with_piholes.xyz`

Original molecule plus validated pi-holes written as dummy `X` atoms.

### `pihole_plots/`

SVG visualizations of selected candidates. By default, only validated candidates are plotted.

---

## 🖼️ Example Output

The included `data/benzene/` and `data/f6benzene/` folders demonstrate the expected input/output layout.

Example generated files include:

* `data/benzene/pi_hole_candidates.csv`
* `data/benzene/pi_holes_validated.txt`
* `data/benzene/molecule_with_piholes.xyz`
* `data/benzene/pihole_plots/pihole_01_above_idx4584.svg`

---

## 📚 Documentation

Each module contains docstrings for public functions. After installing the package, you can inspect modules with:

```bash
pydoc vhm.pipeline
pydoc vhm.sector
pydoc vhm.mls
```

Or explore the package through an IDE.

---

## 📌 Notes

* Reference atom indices use 1-based indexing.
* XYZ and VTX coordinates are converted to angstrom internally.
* ESP values are assumed to be in kcal/mol.
* The local z axis is defined from the fitted reference plane normal.
* `AUTO_ORIENTATION_ATOM` can flip the positive local z direction toward a chosen atom.
* `continue_on_error = False` stops on the first failed target. Set it to `True` to process remaining targets and report failures at the end.
* Output files are overwritten when the same system is rerun.

---

## 🧪 Applications

VHM is intended for workflows involving:

* pi-hole and sigma-hole analysis
* directional non-covalent interaction studies
* electrostatic potential surface analysis
* high-throughput comparison of related molecular systems
* transition-metal and main-group molecular ESP characterization

---

## 📄 Citation

If you use this software, please cite:

*Quantitative Analysis of Metal-Centered Pi-Holes in {TM(cyclen)}2+ Complexes*

Journal: *ACS Inorganic Chemistry*

Year: 2026

DOI: 10.1021/acs.inorgchem.6c01328

---

## 👤 Author

Lucas Gian Fachini - *PhD Candidate in Inorganic and Theoretical Chemistry*  
[GitHub: lgfachini](https://github.com/lgfachini)

---

## ⚖️ License

This project is licensed under the GPL-3 License. See [LICENSE](LICENSE).

---

## 💡 Acknowledgments

This project uses concepts from:

* electrostatic potential surface analysis
* pi-hole and sigma-hole theory
* directional non-covalent interactions
* local surface fitting and curvature analysis
* sector-based anisotropy descriptors

---

## ✅ Tests

```bash
pip install -e ".[dev]"
pytest
```
