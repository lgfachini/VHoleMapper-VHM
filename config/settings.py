from __future__ import annotations

from pathlib import Path

from vhm import TargetSpec


# Edit this file for your calculations.
#
# Run from the project root:
#   python main.py
#   python -m vhm


# ---------------------------------------------------------------------------
# Dataset selection
# ---------------------------------------------------------------------------

RUN_MODE = "auto"  # "auto" or "manual"

DATA_DIR = Path("data")
AUTO_FOLDER_GLOB = "*"
AUTO_REFERENCE_ATOMS = (1, 3, 6, 10)
AUTO_ORIENTATION_ATOM = None
AUTO_XYZ_FILENAME = None
AUTO_VTX_FILENAME = "vtx.txt"

MANUAL_TARGETS = (
    TargetSpec(
        folder=Path("data/example_folder"),
        reference_atom_indices_1based=(1, 2, 3, 4),
        orientation_atom_index_1based=None,
        xyz_filename=None,
        vtx_filename="vtx.txt",
    ),
)


# ---------------------------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------------------------

ANALYSIS_OPTIONS = {
    "xyz_coord_unit": "angstrom",
    "vtx_coord_unit": "bohr",
    "max_distance_from_axis": 0.25,
    "z_search_min": 0.8,
    "z_search_max": 3.0,
    "local_max_radius": 0.5,
    "strict_local_max": True,
    "min_local_neighbors": 8,
    "require_positive_potential": False,
    "top_n_candidates_per_side": None,
    "do_mls_test": False,
    "mls_radius": 1.50,
    "mls_min_neighbors": 10,
    "mls_sigma": 1.20,
    "mls_include_candidate": False,
    "tangential_grad_norm_thr": 200.0,
    "tangential_hess_eig_max_thr": 0.0,
    "allow_semidefinite_tangential_maximum": True,
    "use_tangential_stationary_shift_test": False,
    "max_tangential_stationary_shift": 0.50,
    "angular_sectors": 8,
    "sector_inner_radius": 0.25,
    "sector_outer_radius": 3.50,
    "sector_half_height": 0.50,
    "min_points_per_sector": 50,
    "delta_thr": -100.0,
    "plot_figsize": (8, 8),
    "plot_atom_labels": True,
    "plot_sector_ring_radius": None,
    "plot_show_all_shell_points": False,
    "plot_only_validated": True,
    "print_debug": True,
    "continue_on_error": False,
}
