from __future__ import annotations

from pathlib import Path

from vhm import TargetSpec


# Edit this file for your calculations.
#
# Run from the project root:
#   python main.py
#   python -m vhm
#
# After installing with `pip install -e .`, you can also run:
#   vhm


# ---------------------------------------------------------------------------
# Dataset selection
# ---------------------------------------------------------------------------

# Select how targets are chosen.
#   "auto"   -> scan DATA_DIR and analyze matching direct subfolders.
#   "manual" -> analyze only the folders listed in MANUAL_TARGETS.
RUN_MODE = "auto"

# Root folder containing one subfolder per molecular system.
# Relative paths are resolved from the project root.
DATA_DIR = Path("data")

# Glob pattern used in auto mode to select folders inside DATA_DIR.
# Examples:
#   "*"       -> all direct subfolders
#   "conf*"   -> folders whose names start with "conf"
#   "*benzene" -> folders whose names end with "benzene"
AUTO_FOLDER_GLOB = "*"

# Four 1-based atom indices used to fit the reference plane and define the
# local frame. Choose atoms that represent the ring, face, or molecular plane
# relative to which the pi-hole search should be performed.
AUTO_REFERENCE_ATOMS = (1, 3, 6, 10)

# Optional 1-based atom index used to orient the positive local z direction.
# Set to None to keep the normal direction from the plane fit.
# Set to an atom index to point +z toward that atom.
AUTO_ORIENTATION_ATOM = None

# Optional XYZ filename override for auto mode.
#   None       -> expect <folder_name>.xyz inside each folder.
#   "geom.xyz" -> use the same filename in every selected folder.
AUTO_XYZ_FILENAME = None

# VTX surface filename expected inside each selected folder.
AUTO_VTX_FILENAME = "vtx.txt"

# Manual targets are used only when RUN_MODE = "manual".
# Use this mode when systems have different atom orderings, filenames, or
# reference atoms. Atom indices are 1-based, exactly as in the XYZ file.
MANUAL_TARGETS = (
    TargetSpec(
        # Folder containing the XYZ and VTX files.
        folder=Path("data/example_folder"),

        # Four atoms defining the local reference plane.
        reference_atom_indices_1based=(1, 2, 3, 4),

        # Optional atom used to orient +z. Use None if no orientation is needed.
        orientation_atom_index_1based=None,

        # Use None to expect <folder_name>.xyz, or provide a filename string.
        xyz_filename=None,

        # VTX surface filename inside this folder.
        vtx_filename="vtx.txt",
    ),
)


# ---------------------------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------------------------

ANALYSIS_OPTIONS = {
    # -----------------------------------------------------------------------
    # Coordinate units
    # -----------------------------------------------------------------------
    # Unit used in the XYZ geometry file. Accepted values: "angstrom", "bohr".
    "xyz_coord_unit": "angstrom",

    # Unit used in the VTX surface file. Accepted values: "angstrom", "bohr".
    # Coordinates are converted to angstrom internally.
    "vtx_coord_unit": "bohr",

    # -----------------------------------------------------------------------
    # Axial search region
    # -----------------------------------------------------------------------
    # Maximum radial distance from the local z axis, in angstrom.
    # Smaller values make the search more centered on the selected axis.
    "max_distance_from_axis": 0.25,

    # Minimum absolute z distance from the reference plane, in angstrom.
    # Points closer than this to the plane are ignored.
    "z_search_min": 0.8,

    # Maximum absolute z distance from the reference plane, in angstrom.
    # Points farther than this from the plane are ignored.
    "z_search_max": 3.0,

    # -----------------------------------------------------------------------
    # Discrete local maximum search
    # -----------------------------------------------------------------------
    # Neighbor radius, in angstrom, used to decide whether a surface point is a
    # local ESP maximum.
    "local_max_radius": 0.5,

    # True requires V(candidate) > all neighbor potentials.
    # False allows V(candidate) >= neighbor potentials.
    "strict_local_max": True,

    # Minimum number of neighboring surface points required inside
    # local_max_radius. Candidates with too few neighbors are ignored.
    "min_local_neighbors": 8,

    # If True, ignore candidates with non-positive ESP values.
    # Leave False when analyzing systems where meaningful pi-hole descriptors
    # may occur on a shifted or charged ESP scale.
    "require_positive_potential": False,

    # Optional cap on how many local maxima are kept per side.
    # None keeps every detected local maximum.
    "top_n_candidates_per_side": None,

    # -----------------------------------------------------------------------
    # Moving Least Squares (MLS) validation
    # -----------------------------------------------------------------------
    # Enable derivative-based validation around each candidate.
    # This is stricter and slower than sector-only validation.
    "do_mls_test": False,

    # Radius, in angstrom, of the local neighborhood used for MLS fitting.
    "mls_radius": 1.50,

    # Minimum number of neighboring points needed for a stable MLS fit.
    "mls_min_neighbors": 10,

    # Gaussian weighting width for MLS fitting, in angstrom.
    # Smaller values weight nearby points more strongly.
    "mls_sigma": 1.20,

    # If True, include the candidate point itself in the MLS fit.
    # False avoids forcing the fitted surface through the candidate value.
    "mls_include_candidate": False,

    # -----------------------------------------------------------------------
    # Tangential derivative criteria used by MLS
    # -----------------------------------------------------------------------
    # Maximum allowed tangential gradient norm at the candidate.
    "tangential_grad_norm_thr": 200.0,

    # Maximum allowed tangential Hessian eigenvalue when strict negative
    # curvature is required.
    "tangential_hess_eig_max_thr": 0.0,

    # True accepts negative semidefinite tangential curvature, allowing flat
    # maxima within numerical tolerance. False requires stricter curvature.
    "allow_semidefinite_tangential_maximum": True,

    # If True, require the stationary point predicted by the tangential
    # quadratic fit to be close to the candidate.
    "use_tangential_stationary_shift_test": False,

    # Maximum accepted distance, in angstrom, between the candidate and the
    # fitted tangential stationary point when the shift test is enabled.
    "max_tangential_stationary_shift": 0.50,

    # -----------------------------------------------------------------------
    # Sector validation
    # -----------------------------------------------------------------------
    # Number of angular sectors around each candidate.
    "angular_sectors": 8,

    # Inner radius, in angstrom, of the annular sector sampling region.
    "sector_inner_radius": 0.25,

    # Outer radius, in angstrom, of the annular sector sampling region.
    "sector_outer_radius": 3.50,

    # Half-height, in angstrom, of the slab sampled around the candidate along
    # local z. The full sampled thickness is 2 * sector_half_height.
    "sector_half_height": 0.50,

    # Minimum number of surface points required in a sector for its average to
    # be used. Sectors with fewer points are ignored.
    "min_points_per_sector": 50,

    # Minimum accepted V_hole,mean value, in kcal/mol.
    # A candidate validates when V_hole,mean >= delta_thr.
    "delta_thr": -100.0,

    # -----------------------------------------------------------------------
    # Plotting
    # -----------------------------------------------------------------------
    # Matplotlib figure size for each SVG plot, in inches.
    "plot_figsize": (8, 8),

    # Show labels for the four reference atoms in each planar plot.
    "plot_atom_labels": True,

    # Radius used to place sector labels and anisotropy vectors.
    # None uses the midpoint between sector_inner_radius and sector_outer_radius.
    "plot_sector_ring_radius": None,

    # If True, plot all surface points used in the sector shell. Useful for
    # debugging sector population, but can make dense plots visually heavy.
    "plot_show_all_shell_points": False,

    # If True, create SVG plots only for validated candidates.
    # False plots every analyzed candidate.
    "plot_only_validated": True,

    # -----------------------------------------------------------------------
    # Runtime behavior
    # -----------------------------------------------------------------------
    # Print detailed progress messages while processing each target.
    "print_debug": True,

    # False stops immediately on the first failed target.
    # True attempts all targets, then reports failures at the end.
    "continue_on_error": False,
}
