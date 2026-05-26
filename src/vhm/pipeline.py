from __future__ import annotations

from typing import List, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .config import AnalysisConfig, AnalysisTarget, MoleculePaths, SigmaTargetSpec, TargetSpec, make_paths
from .geometry import build_bond_frame, build_local_frame, to_local
from .io_utils import read_vtx, read_xyz
from .mls import mls_candidate_ok, mls_local_fit
from .models import CandidateResult, LocalFrame
from .outputs import format_optional, save_candidates_csv, save_molecule_with_piholes_xyz, save_validated_txt
from .plotting import save_hole_plots
from .search import find_local_maxima, get_side_mask, get_sigma_mask
from .sector import sector_validation


AnalysisReturn = Tuple[
    List[CandidateResult],
    List[str],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    LocalFrame,
    MoleculePaths,
]


def read_target_inputs(
    target: AnalysisTarget,
    config: AnalysisConfig,
) -> Tuple[MoleculePaths, List[str], np.ndarray, np.ndarray, np.ndarray]:
    """Resolve paths and read molecule/surface inputs for one target."""
    paths = make_paths(target, config.hole_type)
    molecule_name = paths.folder.name

    symbols, atom_coords = read_xyz(paths.xyz_file, config.xyz_coord_unit, config.bohr_to_ang)
    vtx_coords, potentials = read_vtx(paths.vtx_file, config.vtx_coord_unit, config.bohr_to_ang)

    if config.print_debug:
        print(f"\n[INFO] Molecule: {molecule_name}")
        print(f"[INFO] Folder: {paths.folder}")
        print(f"[INFO] XYZ atoms read: {len(symbols)}")
        print(f"[INFO] VTX points read: {len(vtx_coords)}")

    return paths, symbols, atom_coords, vtx_coords, potentials


def analyze_region_candidates(
    local_coords: np.ndarray,
    potentials: np.ndarray,
    tree: cKDTree,
    region_name: str,
    candidate_indices: np.ndarray,
    frame: LocalFrame,
    config: AnalysisConfig,
) -> List[CandidateResult]:
    """Analyze local maxima from one preselected axial/cylindrical region."""
    if config.print_debug:
        print(f"[INFO] Region {region_name}: {len(candidate_indices)} points in the search region.")

    maxima = find_local_maxima(local_coords, potentials, candidate_indices, config)

    if config.print_debug:
        print(f"[INFO] Region {region_name}: {len(maxima)} local maxima found.")

    results: List[CandidateResult] = []
    for point_index in maxima:
        mls_fit = None
        pre_reason = ""

        if config.do_mls_test:
            mls_fit = mls_local_fit(local_coords, potentials, tree, point_index, config)
            ok, mls_reason = mls_candidate_ok(mls_fit, config)
            pre_reason = mls_reason

            if not ok:
                if config.print_debug:
                    print(f"[INFO] Candidate {point_index} ({region_name}) rejected by tangential MLS: {mls_reason}")
                continue

        result = sector_validation(
            local_coords=local_coords,
            potentials=potentials,
            point_index=point_index,
            side=region_name,
            frame=frame,
            config=config,
            mls_fit=mls_fit,
            pre_reason=pre_reason,
        )
        results.append(result)

    return results


def analyze_pi_molecule(target: TargetSpec, config: AnalysisConfig) -> AnalysisReturn:
    """Run the complete pi-hole analysis for one target folder."""
    paths, symbols, atom_coords, vtx_coords, potentials = read_target_inputs(target, config)

    frame = build_local_frame(
        atom_coords=atom_coords,
        reference_indices_1based=target.reference_atom_indices_1based,
        orientation_atom_index_1based=target.orientation_atom_index_1based,
    )

    if config.print_debug:
        print(f"[INFO] Local frame origin (Angstrom): {frame.origin}")
        print("[INFO] Rotation matrix:")
        print(frame.rotation_matrix)

    local_coords = to_local(vtx_coords, frame)
    tree = cKDTree(local_coords)
    results: List[CandidateResult] = []

    for side in ("above", "below"):
        side_indices = np.where(get_side_mask(local_coords, side, config))[0]
        results.extend(
            analyze_region_candidates(
                local_coords=local_coords,
                potentials=potentials,
                tree=tree,
                region_name=side,
                candidate_indices=side_indices,
                frame=frame,
                config=config,
            )
        )

    return results, symbols, atom_coords, local_coords, potentials, frame, paths


def analyze_sigma_molecule(target: SigmaTargetSpec, config: AnalysisConfig) -> AnalysisReturn:
    """Run the complete sigma-hole analysis for one target folder."""
    paths, symbols, atom_coords, vtx_coords, potentials = read_target_inputs(target, config)

    frame = build_bond_frame(
        atom_coords=atom_coords,
        bond_indices_1based=target.bond_atom_indices_1based,
    )

    if config.print_debug:
        atom_a, atom_b = target.bond_atom_indices_1based
        print(f"[INFO] Sigma bond axis: atom {atom_a} -> atom {atom_b}")
        print(f"[INFO] Local frame origin at terminal atom {atom_b} (Angstrom): {frame.origin}")
        print("[INFO] Rotation matrix:")
        print(frame.rotation_matrix)

    local_coords = to_local(vtx_coords, frame)
    tree = cKDTree(local_coords)
    sigma_indices = np.where(get_sigma_mask(local_coords, config))[0]
    results = analyze_region_candidates(
        local_coords=local_coords,
        potentials=potentials,
        tree=tree,
        region_name="forward",
        candidate_indices=sigma_indices,
        frame=frame,
        config=config,
    )

    return results, symbols, atom_coords, local_coords, potentials, frame, paths


def analyze_molecule(target: AnalysisTarget, config: AnalysisConfig) -> AnalysisReturn:
    """Run the configured hole analysis for one target folder."""
    if config.hole_type == "sigma":
        if not isinstance(target, SigmaTargetSpec):
            raise TypeError("sigma-hole analysis requires SigmaTargetSpec targets.")
        return analyze_sigma_molecule(target, config)

    if not isinstance(target, TargetSpec):
        raise TypeError("pi-hole analysis requires TargetSpec targets.")
    return analyze_pi_molecule(target, config)


def write_outputs(
    target: AnalysisTarget,
    results: List[CandidateResult],
    symbols: List[str],
    atom_coords: np.ndarray,
    local_coords: np.ndarray,
    frame: LocalFrame,
    paths: MoleculePaths,
    config: AnalysisConfig,
) -> None:
    """Write all configured output files for one target."""
    hole_label = "sigma-hole" if config.hole_type == "sigma" else "pi-hole"
    xyz_label = "sigma-holes" if config.hole_type == "sigma" else "pi-holes"
    filename_prefix = "sigmahole" if config.hole_type == "sigma" else "pihole"

    save_candidates_csv(results, paths.output_all_candidates_csv)
    save_validated_txt(results, paths.output_validated_txt, hole_label=hole_label)
    save_molecule_with_piholes_xyz(
        symbols,
        atom_coords,
        results,
        paths.output_hole_xyz,
        only_validated=True,
        hole_label=xyz_label,
    )

    reference_atoms = target.reference_atom_indices_1based if isinstance(target, TargetSpec) else None
    bond_atoms = target.bond_atom_indices_1based if isinstance(target, SigmaTargetSpec) else None
    save_hole_plots(
        atom_symbols=symbols,
        atom_coords_global=atom_coords,
        frame=frame,
        reference_atom_indices_1based=reference_atoms,
        bond_atom_indices_1based=bond_atoms,
        local_coords_surface=local_coords,
        results=results,
        out_dir=paths.output_plot_dir,
        config=config,
        hole_label=hole_label,
        filename_prefix=filename_prefix,
    )


def print_summary(molecule_name: str, results: List[CandidateResult], hole_label: str = "pi-holes") -> None:
    """Print a concise terminal summary."""
    total = len(results)
    validated = [result for result in results if result.validated]

    print("\n" + "=" * 80)
    print(f"SUMMARY: {molecule_name}")
    print("=" * 80)
    print(f"Total candidates analyzed: {total}")
    print(f"Total validated {hole_label}: {len(validated)}")

    if validated:
        print(f"\nValidated {hole_label}:")
        for counter, result in enumerate(sorted(validated, key=lambda item: item.vs_max, reverse=True), start=1):
            lambda_max_tan = max(result.tangential_hess_eig1, result.tangential_hess_eig2)
            print(
                f"{counter:2d}) region={result.side:7s} | "
                f"z={result.z_local:8.4f} Angstrom | "
                f"Vs,max={result.vs_max:10.4f} | "
                f"|grad_tan|={result.tangential_grad_norm:10.4f} | "
                f"lambda_max_tan={lambda_max_tan:10.4f} | "
                f"shift_tan={format_optional(result.tangential_stationary_shift, '.4f'):>8s} | "
                f"Vhole,mean={result.vhole_mean:10.4f}"
            )

    print("=" * 80 + "\n")


def run_batch(config: AnalysisConfig) -> None:
    """Run the analysis for all targets listed in config.targets."""
    failures: List[Tuple[AnalysisTarget, Exception]] = []
    summary_label = "sigma-holes" if config.hole_type == "sigma" else "pi-holes"

    for target in config.targets:
        try:
            results, symbols, atom_coords, local_coords, _potentials, frame, paths = analyze_molecule(target, config)
            write_outputs(target, results, symbols, atom_coords, local_coords, frame, paths, config)
            print_summary(paths.folder.name, results, summary_label)
            print(f"[INFO] CSV saved to: {paths.output_all_candidates_csv}")
            print(f"[INFO] TXT saved to: {paths.output_validated_txt}")
            print(f"[INFO] XYZ saved to: {paths.output_hole_xyz}")
            print(f"[INFO] Plots saved to: {paths.output_plot_dir}")
        except Exception as exc:
            if not config.continue_on_error:
                raise
            failures.append((target, exc))
            print(f"[ERROR] Analysis failed for {target.folder}: {exc}")

    if failures:
        failed = ", ".join(str(target.folder) for target, _exc in failures)
        raise RuntimeError(f"Analysis failed for {len(failures)} target(s): {failed}")
