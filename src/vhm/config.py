from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Union

PathLike = Union[str, Path]


@dataclass(frozen=True)
class TargetSpec:
    """
    Description of one folder to analyze.

    User-editable fields:
    - folder: directory containing the input files.
    - reference_atom_indices_1based: four atoms used to define the reference plane.
      Use 1-based atom indices exactly as they appear in the XYZ file.
    - orientation_atom_index_1based: optional atom used only to orient the positive
      local z axis. Use None when no orientation atom is needed.
    - xyz_filename: optional XYZ filename override. If None, the code expects
      <folder_name>.xyz inside the folder.
    - vtx_filename: VTX filename inside the folder.
    """

    folder: PathLike
    reference_atom_indices_1based: Tuple[int, int, int, int]
    orientation_atom_index_1based: Optional[int] = None
    xyz_filename: Optional[str] = None
    vtx_filename: str = "vtx.txt"

    def __post_init__(self) -> None:
        object.__setattr__(self, "folder", Path(self.folder))
        object.__setattr__(
            self,
            "reference_atom_indices_1based",
            tuple(self.reference_atom_indices_1based),
        )
        if len(self.reference_atom_indices_1based) != 4:
            raise ValueError("reference_atom_indices_1based must contain exactly four atoms.")
        if any(index < 1 for index in self.reference_atom_indices_1based):
            raise ValueError("reference atom indices must use positive 1-based indexing.")
        if self.orientation_atom_index_1based is not None and self.orientation_atom_index_1based < 1:
            raise ValueError("orientation_atom_index_1based must be positive when set.")
        if not str(self.vtx_filename).strip():
            raise ValueError("vtx_filename must be non-empty.")


@dataclass(frozen=True)
class SigmaTargetSpec:
    """
    Description of one sigma-hole target.

    The bond is ordered. For bond_atom_indices_1based=(A, B), the sigma-hole is
    searched beyond atom B along the A -> B bond extension.
    """

    folder: PathLike
    bond_atom_indices_1based: Tuple[int, int]
    xyz_filename: Optional[str] = None
    vtx_filename: str = "vtx.txt"

    def __post_init__(self) -> None:
        object.__setattr__(self, "folder", Path(self.folder))
        object.__setattr__(self, "bond_atom_indices_1based", tuple(self.bond_atom_indices_1based))
        if len(self.bond_atom_indices_1based) != 2:
            raise ValueError("bond_atom_indices_1based must contain exactly two atoms.")
        if any(index < 1 for index in self.bond_atom_indices_1based):
            raise ValueError("bond atom indices must use positive 1-based indexing.")
        if self.bond_atom_indices_1based[0] == self.bond_atom_indices_1based[1]:
            raise ValueError("bond_atom_indices_1based must contain two different atoms.")
        if not str(self.vtx_filename).strip():
            raise ValueError("vtx_filename must be non-empty.")


AnalysisTarget = Union[TargetSpec, SigmaTargetSpec]


@dataclass(frozen=True)
class AnalysisConfig:
    """
    Global analysis parameters.

    Parameters users usually tweak:
    - targets: tuple of TargetSpec objects created in main.py.
    - xyz_coord_unit and vtx_coord_unit: coordinate units in the input files.
    - max_distance_from_axis, z_search_min, z_search_max: axial search region.
    - local_max_radius, strict_local_max, min_local_neighbors: local maximum test.
    - sector_inner_radius, sector_outer_radius, sector_half_height: sector shell.
    - delta_thr: minimum accepted V_hole,mean value.
    - plot_only_validated and plot_show_all_shell_points: plot behavior.

    Parameters users usually leave unchanged:
    - bohr_to_ang: unit conversion constant.
    - MLS parameters, unless the tangential derivative test is enabled.
    """

    # -------------------------------------------------------------------------
    # Targets
    # -------------------------------------------------------------------------
    targets: Tuple[AnalysisTarget, ...]
    hole_type: str = "pi"

    # -------------------------------------------------------------------------
    # Coordinate units
    # -------------------------------------------------------------------------
    # Accepted values: "angstrom" or "bohr".
    xyz_coord_unit: str = "angstrom"
    vtx_coord_unit: str = "bohr"
    bohr_to_ang: float = 0.529177210903

    # -------------------------------------------------------------------------
    # Axial search geometry
    # -------------------------------------------------------------------------
    # max_distance_from_axis: radial cutoff from the local z axis, in angstrom.
    # z_search_min/z_search_max: axial window, in angstrom, searched on both
    # sides of the reference plane.
    max_distance_from_axis: float = 0.7
    z_search_min: float = 0.5
    z_search_max: float = 3.5

    # -------------------------------------------------------------------------
    # Discrete local maximum search
    # -------------------------------------------------------------------------
    # local_max_radius: neighbor radius used to test whether a point is a local
    # maximum on the sampled surface.
    # strict_local_max: True requires V(candidate) > V(neighbors). False allows >=.
    # min_local_neighbors: minimum number of neighbors inside local_max_radius.
    # require_positive_potential: optional prefilter requiring V(candidate) > 0.
    # top_n_candidates_per_side: if set, keeps only the highest candidates per side.
    local_max_radius: float = 0.25
    strict_local_max: bool = True
    min_local_neighbors: int = 8
    require_positive_potential: bool = False
    top_n_candidates_per_side: Optional[int] = None

    # -------------------------------------------------------------------------
    # Local derivative test using MLS
    # -------------------------------------------------------------------------
    # do_mls_test: enable or disable the tangential MLS derivative test.
    # mls_radius: radius of the neighborhood used for the local quadratic fit.
    # mls_min_neighbors: minimum number of points required for the fit.
    # mls_sigma: Gaussian weighting width, in angstrom.
    # mls_include_candidate: include the candidate itself in the MLS fit.
    do_mls_test: bool = False
    mls_radius: float = 1.50
    mls_min_neighbors: int = 10
    mls_sigma: float = 1.20
    mls_include_candidate: bool = False

    # -------------------------------------------------------------------------
    # Tangential MLS criteria
    # -------------------------------------------------------------------------
    # tangential_grad_norm_thr: maximum accepted tangential gradient norm.
    # tangential_hess_eig_max_thr: maximum allowed Hessian eigenvalue when strict
    # negative curvature is required.
    # allow_semidefinite_tangential_maximum: True accepts negative semidefinite
    # tangential curvature within numerical tolerance.
    # use_tangential_stationary_shift_test: optionally require the fitted
    # stationary point to be close to the candidate.
    tangential_grad_norm_thr: float = 200.0
    tangential_hess_eig_max_thr: float = 0.0
    allow_semidefinite_tangential_maximum: bool = True
    use_tangential_stationary_shift_test: bool = False
    max_tangential_stationary_shift: float = 0.50

    # -------------------------------------------------------------------------
    # Sector validation
    # -------------------------------------------------------------------------
    # angular_sectors: number of angular sectors around each candidate.
    # sector_inner_radius/sector_outer_radius: radial shell around the candidate.
    # sector_half_height: half-height of the shell along local z.
    # min_points_per_sector: sectors with fewer points are ignored.
    # delta_thr: validation threshold for V_hole,mean, in kcal/mol.
    angular_sectors: int = 8
    sector_inner_radius: float = 0.50
    sector_outer_radius: float = 3.00
    sector_half_height: float = 0.50
    min_points_per_sector: int = 50
    delta_thr: float = -100.0

    # -------------------------------------------------------------------------
    # Plot controls
    # -------------------------------------------------------------------------
    # plot_only_validated: if True, only validated candidates are plotted.
    # plot_show_all_shell_points: useful for checking sector populations.
    plot_figsize: Tuple[float, float] = (8, 8)
    plot_atom_labels: bool = True
    plot_sector_ring_radius: Optional[float] = None
    plot_show_all_shell_points: bool = False
    plot_only_validated: bool = True

    # -------------------------------------------------------------------------
    # Runtime controls
    # -------------------------------------------------------------------------
    print_debug: bool = True
    continue_on_error: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "targets", tuple(self.targets))
        if not self.targets:
            raise ValueError("targets must contain at least one analysis target.")

        hole_type = self.hole_type.lower().strip()
        object.__setattr__(self, "hole_type", hole_type)
        if hole_type not in {"pi", "sigma"}:
            raise ValueError("hole_type must be either 'pi' or 'sigma'.")
        if hole_type == "pi" and not all(isinstance(target, TargetSpec) for target in self.targets):
            raise ValueError("pi-hole analysis requires TargetSpec targets.")
        if hole_type == "sigma" and not all(isinstance(target, SigmaTargetSpec) for target in self.targets):
            raise ValueError("sigma-hole analysis requires SigmaTargetSpec targets.")

        units = {"angstrom", "bohr"}
        if self.xyz_coord_unit.lower().strip() not in units:
            raise ValueError("xyz_coord_unit must be 'angstrom' or 'bohr'.")
        if self.vtx_coord_unit.lower().strip() not in units:
            raise ValueError("vtx_coord_unit must be 'angstrom' or 'bohr'.")
        if self.bohr_to_ang <= 0.0:
            raise ValueError("bohr_to_ang must be positive.")

        positive_fields = {
            "max_distance_from_axis": self.max_distance_from_axis,
            "z_search_min": self.z_search_min,
            "z_search_max": self.z_search_max,
            "local_max_radius": self.local_max_radius,
            "mls_radius": self.mls_radius,
            "mls_sigma": self.mls_sigma,
            "sector_inner_radius": self.sector_inner_radius,
            "sector_outer_radius": self.sector_outer_radius,
            "sector_half_height": self.sector_half_height,
            "max_tangential_stationary_shift": self.max_tangential_stationary_shift,
        }
        for name, value in positive_fields.items():
            if value <= 0.0:
                raise ValueError(f"{name} must be positive.")

        if self.z_search_min >= self.z_search_max:
            raise ValueError("z_search_min must be smaller than z_search_max.")
        if self.sector_inner_radius >= self.sector_outer_radius:
            raise ValueError("sector_inner_radius must be smaller than sector_outer_radius.")
        if self.min_local_neighbors < 1:
            raise ValueError("min_local_neighbors must be at least 1.")
        if self.mls_min_neighbors < 1:
            raise ValueError("mls_min_neighbors must be at least 1.")
        if self.angular_sectors < 2:
            raise ValueError("angular_sectors must be at least 2.")
        if self.min_points_per_sector < 1:
            raise ValueError("min_points_per_sector must be at least 1.")
        if self.top_n_candidates_per_side is not None and self.top_n_candidates_per_side < 1:
            raise ValueError("top_n_candidates_per_side must be positive when set.")
        if len(self.plot_figsize) != 2 or any(size <= 0.0 for size in self.plot_figsize):
            raise ValueError("plot_figsize must contain two positive values.")


@dataclass(frozen=True)
class MoleculePaths:
    """Resolved input and output paths for one target folder."""

    folder: Path
    vtx_file: Path
    xyz_file: Path
    output_all_candidates_csv: Path
    output_validated_txt: Path
    output_hole_xyz: Path
    output_plot_dir: Path


def make_paths(target: AnalysisTarget, hole_type: str = "pi") -> MoleculePaths:
    """
    Build input and output paths for one target.

    Default expected folder structure:
        target_folder/
            vtx.txt
            target_folder_name.xyz

    The XYZ and VTX filenames can be overridden in TargetSpec.
    """
    folder = Path(target.folder)
    if not folder.exists():
        raise FileNotFoundError(f"Target folder not found: {folder}")
    if not folder.is_dir():
        raise NotADirectoryError(f"Target path is not a directory: {folder}")

    vtx_file = folder / target.vtx_filename
    xyz_file = folder / target.xyz_filename if target.xyz_filename else folder / f"{folder.name}.xyz"

    if not vtx_file.exists():
        raise FileNotFoundError(f"VTX file not found: {vtx_file}")
    if not xyz_file.exists():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")

    if hole_type == "pi":
        candidates_csv = "pi_hole_candidates.csv"
        validated_txt = "pi_holes_validated.txt"
        hole_xyz = "molecule_with_piholes.xyz"
        plot_dir_name = "pihole_plots"
    else:
        candidates_csv = "sigma_hole_candidates.csv"
        validated_txt = "sigma_holes_validated.txt"
        hole_xyz = "molecule_with_sigmaholes.xyz"
        plot_dir_name = "sigmahole_plots"

    return MoleculePaths(
        folder=folder,
        vtx_file=vtx_file,
        xyz_file=xyz_file,
        output_all_candidates_csv=folder / candidates_csv,
        output_validated_txt=folder / validated_txt,
        output_hole_xyz=folder / hole_xyz,
        output_plot_dir=folder / plot_dir_name,
    )
