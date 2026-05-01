# 🧲 VHoleMapper (VHM): A π-Hole Detection and Analysis Toolkit

VHoleMapper (VHM) is a Python-based toolkit for detecting, validating, and analyzing **π-holes** on molecular electrostatic potential (ESP) surfaces.

It is designed for high-throughput workflows and supports automated analysis of multiple molecular systems using ESP data (e.g., `vtx.txt`) and molecular geometries (`.xyz`).

---

## ⚙️ Features

* 🔍 **Automatic π-hole detection** via local maxima search
* 📐 **Local reference frame construction** from user-defined atoms
* 📊 **Sector-based validation** of π-hole anisotropy
* 🧮 Optional **MLS (Moving Least Squares)** derivative validation
* 🧭 **Tangential gradient and Hessian analysis**
* 📁 Batch processing of multiple molecular systems
* 📈 **2D planar visualization** of π-holes and anisotropy
* 📄 Output formats:

  * CSV (all candidates)
  * TXT (validated π-holes)
  * XYZ (molecule + dummy π-hole atoms)
  * SVG (sector plots)

---

## 📁 Project Structure

```
VHoleMapper/
├── main.py             # Entry point (user controls everything here)
├── config.py           # Global configuration and parameters
├── pipeline.py         # Main workflow orchestration
├── geometry.py         # Local frame construction and transformations
├── search.py           # Local maxima detection
├── mls.py              # MLS fitting and derivative analysis
├── sector.py           # Sector-based validation
├── plotting.py         # Visualization tools
├── outputs.py          # File export (CSV, TXT, XYZ)
├── io_utils.py         # File parsing (XYZ, VTX)
├── models.py           # Data structures
└── data/               # Input folders (one system per folder)
```

---

## ▶️ How to Use

### 1. Install Dependencies

```bash
pip install numpy scipy matplotlib
```

---

### 2. Prepare Input Data

Each system must be inside its own folder:

```
data/
└── system_name/
    ├── vtx.txt
    └── system_name.xyz
```

---

### 3. Configure `main.py`

Choose one mode:

#### 🔹 Automatic mode

```python
RUN_MODE = "auto"
DATA_DIR = Path("data")
```

Uses the same reference atoms for all folders.

#### 🔹 Manual mode

```python
RUN_MODE = "manual"

MANUAL_TARGETS = (
    TargetSpec(
        folder=Path("data/system"),
        reference_atom_indices_1based=(1,2,3,4),
        orientation_atom_index_1based=None,
    ),
)
```

Allows different definitions per system.

---

### 4. Run

```bash
python main.py
```

---

## 📊 Outputs

For each system:

* `pi_hole_candidates.csv` → all detected candidates
* `pi_holes_validated.txt` → validated π-holes
* `molecule_with_piholes.xyz` → geometry + dummy atoms
* `pihole_plots/` → SVG visualizations

---

## 📚 Methodology

The workflow consists of:

1. **Local frame definition** from 4 reference atoms
2. **Axial search** for ESP maxima
3. Optional **MLS derivative validation**
4. **Sector-based validation** using angular partitioning
5. Extraction of:

   * Vₛ,max
   * Vhole,mean
   * ΔVhole
   * anisotropy angle

---

## 📌 Notes

* Atom indices are **1-based**
* Coordinates are internally handled in **angstrom**
* ESP values are assumed in **kcal/mol**
* Sector validation is robust against charged systems

---
## 📄 Citation

If you happen to use this software, please cite:

Quantitative Analysis of Metal-Centered π-Holes in {TM(cyclen)}2+ Complexes
journal: Inorganic Chemistry
year: 2026
doi: 10.1021/acs.inorgchem.6c01328

---

## 👨‍🔬 Author

Lucas Gian Fachini
*PhD Candidate in Inorganic and Theoretical Chemistry*

GitHub: https://github.com/lgfachini

---

## 📄 License

This project is licensed under the GPL-3.0 License.

---

## 💡 Acknowledgments

Concepts used in this project include:

* Electrostatic potential analysis
* σ-hole and π-hole theory
* Directional non-covalent interactions
* Local surface fitting and curvature analysis

---
