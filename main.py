from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Tuple

from config import AnalysisConfig, TargetSpec
from pipeline import run_batch


# =============================================================================
# USER INPUT: DATASET SELECTION
# =============================================================================
#
# This section controls which folders are analyzed and how each folder is
# interpreted by the workflow.
#
# Two modes are available:
#
#   "auto"
#       Automatically scans DATA_DIR and analyzes all matching direct
#       subfolders. The same reference atoms and orientation atom are applied
#       to every folder. Use this mode when all folders have the same atom
#       ordering/topology.
#
#   "manual"
#       Analyzes only the folders listed in MANUAL_TARGETS. Each folder can
#       have its own reference atoms, orientation atom, and file names. Use
#       this mode when the dataset contains different molecular systems or
#       different atom orderings.
#
# All atom indices below use 1-based indexing, meaning that the first atom in
# the XYZ file is atom 1, not atom 0.
# =============================================================================

RUN_MODE = "auto"  # Allowed values: "auto" or "manual"


# =============================================================================
# AUTO MODE SETTINGS
# =============================================================================
#
# These options are used only when RUN_MODE = "auto".
#
# Expected folder layout:
#
#   DATA_DIR/
#       system_1/
#           vtx.txt
#           system_1.xyz
#       system_2/
#           vtx.txt
#           system_2.xyz
#
# If AUTO_XYZ_FILENAME is None, the workflow expects each XYZ file to be named
# after its folder, for example:
#
#   data/example/example.xyz
#
# If all folders use the same XYZ filename, set AUTO_XYZ_FILENAME to that name.
# =============================================================================

DATA_DIR = Path("data")

# Glob pattern used to select folders inside DATA_DIR.
# Examples:
#   "*"       -> all direct subfolders
#   "conf*"   -> only folders whose names start with "conf"
#   "*Cu*"    -> only folders whose names contain "Cu"
AUTO_FOLDER_GLOB = "*"

# Four atoms used to define the local reference plane.
# These atoms should be chosen so that they represent the plane or surface
# relative to which the axial search is performed.
AUTO_REFERENCE_ATOMS = (1, 3, 6, 10)

# Optional atom used to orient the positive local z direction.
# If None, the plane normal direction is kept as returned by the plane fit.
# If an atom index is provided, +z is oriented toward that atom.
AUTO_ORIENTATION_ATOM = None

# Optional XYZ filename override.
# Use None to assume <folder_name>.xyz.
AUTO_XYZ_FILENAME = None

# VTX filename expected inside each folder.
AUTO_VTX_FILENAME = "vtx.txt"


# =============================================================================
# MANUAL MODE SETTINGS
# =============================================================================
#
# These options are used only when RUN_MODE = "manual".
#
# Add one TargetSpec for each folder you want to analyze. This allows each
# target to have its own reference atoms, orientation atom, and file names.
#
# The example below can be edited, duplicated, or removed.
# =============================================================================

MANUAL_TARGETS = (
    TargetSpec(
        folder=Path("data/example_folder"),

        # Four atoms defining the local reference plane.
        reference_atom_indices_1based=(1, 2, 3, 4),

        # Optional atom used to orient +z. Use None if no orientation atom is needed.
        orientation_atom_index_1based=None,

        # Use None to assume <folder_name>.xyz, or provide a filename string.
        xyz_filename=None,

        # Name of the VTX file inside this folder.
        vtx_filename="vtx.txt",
    ),
)


# =============================================================================
# FOLDER DISCOVERY HELPERS
# =============================================================================
#
# The functions below convert the user selections above into TargetSpec objects.
# Usually, you do not need to edit this section.
# =============================================================================

def discover_folders(data_dir: Path, pattern: str = "*") -> Tuple[Path, ...]:
    """
    Return direct subfolders inside data_dir matching a glob pattern.

    Parameters
    ----------
    data_dir
        Root directory containing one subfolder per system.
    pattern
        Glob pattern used to filter folders.

    Returns
    -------
    tuple[Path, ...]
        Sorted tuple of matching folders.
    """
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    if not data_dir.is_dir():
        raise NotADirectoryError(f"DATA_DIR is not a directory: {data_dir}")

    folders = tuple(sorted(path for path in data_dir.glob(pattern) if path.is_dir()))

    if not folders:
        raise RuntimeError(f"No folders found in {data_dir} with pattern {pattern!r}")

    return folders


def build_auto_targets(
    folders: Iterable[Path],
    reference_atoms: Tuple[int, int, int, int],
    orientation_atom: Optional[int],
    xyz_filename: Optional[str],
    vtx_filename: str,
) -> Tuple[TargetSpec, ...]:
    """
    Create TargetSpec entries for folders discovered in auto mode.
    """
    return tuple(
        TargetSpec(
            folder=folder,
            reference_atom_indices_1based=reference_atoms,
            orientation_atom_index_1based=orientation_atom,
            xyz_filename=xyz_filename,
            vtx_filename=vtx_filename,
        )
        for folder in folders
    )


def selected_targets() -> Tuple[TargetSpec, ...]:
    """
    Resolve the analysis targets according to RUN_MODE.
    """
    if RUN_MODE == "auto":
        folders = discover_folders(DATA_DIR, AUTO_FOLDER_GLOB)
        return build_auto_targets(
            folders=folders,
            reference_atoms=AUTO_REFERENCE_ATOMS,
            orientation_atom=AUTO_ORIENTATION_ATOM,
            xyz_filename=AUTO_XYZ_FILENAME,
            vtx_filename=AUTO_VTX_FILENAME,
        )

    if RUN_MODE == "manual":
        if not MANUAL_TARGETS:
            raise RuntimeError("RUN_MODE is 'manual', but MANUAL_TARGETS is empty.")
        return tuple(MANUAL_TARGETS)

    raise ValueError("RUN_MODE must be either 'auto' or 'manual'.")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main() -> None:
    """
    Run the full pi-hole analysis workflow.
    
    ```
    Workflow:
    1. Select folders and atom definitions.
    2. Build AnalysisConfig with numerical parameters.
    3. Print the selected targets.
    4. Run the batch analysis.
    """
    targets = selected_targets()
    
    config = AnalysisConfig(
        targets=targets,
    
        # ---------------------------------------------------------------------
        # COORDINATE UNITS
        # ---------------------------------------------------------------------
        # Units of the input files:
        # - xyz_coord_unit: unit used in the XYZ geometry file
        # - vtx_coord_unit: unit used in the VTX surface file
        #
        # These must match how the files were generated.
        xyz_coord_unit="angstrom",
        vtx_coord_unit="bohr",
    
        # ---------------------------------------------------------------------
        # AXIAL SEARCH REGION (π-hole candidate region)
        # ---------------------------------------------------------------------
        # Defines the cylindrical region along the local z-axis where
        # candidates are searched.
        #
        # max_distance_from_axis:
        #   Maximum radial distance from the z-axis (Å)
        #
        # z_search_min / z_search_max:
        #   Defines the allowed axial range (distance from the reference plane)
        #
        # Smaller values → stricter, more localized search
        # Larger values → broader search, more candidates (but more noise)
        max_distance_from_axis=0.25,
        z_search_min=0.8,
        z_search_max=3.0,
    
        # ---------------------------------------------------------------------
        # DISCRETE LOCAL MAXIMUM DETECTION
        # ---------------------------------------------------------------------
        # Controls how maxima are identified on the ESP surface.
        #
        # local_max_radius:
        #   Radius used to define the neighborhood for comparison
        #
        # strict_local_max:
        #   True  → candidate must be strictly greater than all neighbors
        #   False → allows equal values
        #
        # min_local_neighbors:
        #   Minimum number of neighbors required to validate a maximum
        #
        # require_positive_potential:
        #   If True, only considers points with positive ESP
        #
        # top_n_candidates_per_side:
        #   Limits number of candidates per side (None = no limit)
        local_max_radius=0.5,
        strict_local_max=True,
        min_local_neighbors=8,
        require_positive_potential=False,
        top_n_candidates_per_side=None,
    
        # ---------------------------------------------------------------------
        # MLS (MOVING LEAST SQUARES) VALIDATION
        # ---------------------------------------------------------------------
        # Optional derivative-based validation of candidates.
        #
        # do_mls_test:
        #   Enables/disables MLS validation
        #
        # mls_radius:
        #   Neighborhood radius used for local fitting (Å)
        #
        # mls_min_neighbors:
        #   Minimum points required for a stable fit
        #
        # mls_sigma:
        #   Gaussian weighting parameter (controls locality)
        #
        # mls_include_candidate:
        #   Whether to include the candidate point in the fit
        do_mls_test=False,
        mls_radius=1.50,
        mls_min_neighbors=10,
        mls_sigma=1.20,
        mls_include_candidate=False,
    
        # ---------------------------------------------------------------------
        # TANGENTIAL DERIVATIVE CRITERIA (MLS)
        # ---------------------------------------------------------------------
        # Criteria applied to the gradient and Hessian in the tangent plane.
        #
        # tangential_grad_norm_thr:
        #   Maximum allowed magnitude of the tangential gradient
        #
        # tangential_hess_eig_max_thr:
        #   Threshold for maximum Hessian eigenvalue (if strict test used)
        #
        # allow_semidefinite_tangential_maximum:
        #   True  → allows flat maxima (semidefinite Hessian)
        #   False → requires strictly negative curvature
        #
        # use_tangential_stationary_shift_test:
        #   Enables additional check based on predicted stationary point
        #
        # max_tangential_stationary_shift:
        #   Maximum allowed displacement (Å)
        tangential_grad_norm_thr=200.0,
        tangential_hess_eig_max_thr=0.0,
        allow_semidefinite_tangential_maximum=True,
        use_tangential_stationary_shift_test=False,
        max_tangential_stationary_shift=0.50,
    
        # ---------------------------------------------------------------------
        # SECTOR VALIDATION (π-hole characterization)
        # ---------------------------------------------------------------------
        # Evaluates anisotropy of the ESP around the candidate.
        #
        # angular_sectors:
        #   Number of angular divisions around the candidate
        #
        # sector_inner_radius / sector_outer_radius:
        #   Radial range used to compute sector averages (Å)
        #
        # sector_half_height:
        #   Thickness of the sampling region along z (Å)
        #
        # min_points_per_sector:
        #   Minimum number of points required per sector
        #
        # delta_thr:
        #   Minimum V_hole,mean required for validation
        angular_sectors=8,
        sector_inner_radius=0.25,
        sector_outer_radius=3.50,
        sector_half_height=0.50,
        min_points_per_sector=50,
        delta_thr=-100.0,
    
        # ---------------------------------------------------------------------
        # PLOTTING SETTINGS
        # ---------------------------------------------------------------------
        # Controls visualization of results.
        #
        # plot_figsize:
        #   Size of generated figures
        #
        # plot_atom_labels:
        #   Show labels for reference atoms
        #
        # plot_sector_ring_radius:
        #   Radius used to draw sector annotations (None = auto)
        #
        # plot_show_all_shell_points:
        #   If True, plots all sampled surface points (can be heavy)
        #
        # plot_only_validated:
        #   If True, plots only validated π-holes
        plot_figsize=(8, 8),
        plot_atom_labels=True,
        plot_sector_ring_radius=None,
        plot_show_all_shell_points=False,
        plot_only_validated=True,
    
        # ---------------------------------------------------------------------
        # DEBUG / LOGGING
        # ---------------------------------------------------------------------
        # Enables verbose output during execution
        print_debug=True,
    )

    print("Targets selected for analysis:")
    for target in config.targets:
        print(
            f"  - {Path(target.folder)} | "
            f"reference atoms = {target.reference_atom_indices_1based} | "
            f"orientation atom = {target.orientation_atom_index_1based}"
        )
    
    run_batch(config)



if __name__ == "__main__":
    main()
