from __future__ import annotations

from typing import List, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .config import AnalysisConfig, MoleculePaths, TargetSpec, make_paths
from .geometry import build_local_frame, to_local
from .io_utils import read_vtx, read_xyz
from .mls import mls_candidate_ok, mls_local_fit
from .models import CandidateResult, LocalFrame
from .outputs import format_optional, save_candidates_csv, save_molecule_with_piholes_xyz, save_validated_txt
from .plotting import save_pihole_plots
from .search import find_local_maxima, get_side_mask
from .sector import sector_validation


def analyze_molecule(
    target: TargetSpec,
    config: AnalysisConfig,
) -> Tuple[List[CandidateResult], List[str], np.ndarray, np.ndarray, np.ndarray, LocalFrame, MoleculePaths]:
    """Run the complete pi-hole analysis for one target folder."""
    paths = make_paths(target)
    molecule_name = paths.folder.name

    symbols, atom_coords = read_xyz(paths.xyz_file, config.xyz_coord_unit, config.bohr_to_ang)
    vtx_coords, potentials = read_vtx(paths.vtx_file, config.vtx_coord_unit, config.bohr_to_ang)

    if config.print_debug:
        print(f"\n[INFO] Molecule: {molecule_name}")
        print(f"[INFO] Folder: {paths.folder}")
        print(f"[INFO] XYZ atoms read: {len(symbols)}")
        print(f"[INFO] VTX points read: {len(vtx_coords)}")

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
        side_mask = get_side_mask(local_coords, side, config)
        side_indices = np.where(side_mask)[0]

        if config.print_debug:
            print(f"[INFO] Side {side}: {len(side_indices)} points in the search region.")

        maxima = find_local_maxima(local_coords, potentials, side_indices, config)
        if maxima is None:
            maxima = []

        if config.print_debug:
            print(f"[INFO] Side {side}: {len(maxima)} local maxima found.")

        for point_index in maxima:
            mls_fit = None
            pre_reason = ""

            if config.do_mls_test:
                mls_fit = mls_local_fit(local_coords, potentials, tree, point_index, config)
                ok, mls_reason = mls_candidate_ok(mls_fit, config)
                pre_reason = mls_reason

                if not ok:
                    if config.print_debug:
                        print(f"[INFO] Candidate {point_index} ({side}) rejected by tangential MLS: {mls_reason}")
                    continue

            result = sector_validation(
                local_coords=local_coords,
                potentials=potentials,
                point_index=point_index,
                side=side,
                frame=frame,
                config=config,
                mls_fit=mls_fit,
                pre_reason=pre_reason,
            )
            results.append(result)

    return results, symbols, atom_coords, local_coords, potentials, frame, paths


def write_outputs(
    target: TargetSpec,
    results: List[CandidateResult],
    symbols: List[str],
    atom_coords: np.ndarray,
    local_coords: np.ndarray,
    frame: LocalFrame,
    paths: MoleculePaths,
    config: AnalysisConfig,
) -> None:
    """Write all configured output files for one target."""
    save_candidates_csv(results, paths.output_all_candidates_csv)
    save_validated_txt(results, paths.output_validated_txt)
    save_molecule_with_piholes_xyz(symbols, atom_coords, results, paths.output_pihole_xyz, only_validated=True)
    save_pihole_plots(
        atom_symbols=symbols,
        atom_coords_global=atom_coords,
        frame=frame,
        reference_atom_indices_1based=target.reference_atom_indices_1based,
        local_coords_surface=local_coords,
        results=results,
        out_dir=paths.output_plot_dir,
        config=config,
    )


def print_summary(molecule_name: str, results: List[CandidateResult]) -> None:
    """Print a concise terminal summary."""
    total = len(results)
    validated = [result for result in results if result.validated]

    print("\n" + "=" * 80)
    print(f"SUMMARY: {molecule_name}")
    print("=" * 80)
    print(f"Total candidates analyzed: {total}")
    print(f"Total validated pi-holes: {len(validated)}")

    if validated:
        print("\nValidated pi-holes:")
        for counter, result in enumerate(sorted(validated, key=lambda item: item.vs_max, reverse=True), start=1):
            lambda_max_tan = max(result.tangential_hess_eig1, result.tangential_hess_eig2)
            print(
                f"{counter:2d}) side={result.side:5s} | "
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
    failures: List[Tuple[TargetSpec, Exception]] = []

    for target in config.targets:
        try:
            results, symbols, atom_coords, local_coords, _potentials, frame, paths = analyze_molecule(target, config)
            write_outputs(target, results, symbols, atom_coords, local_coords, frame, paths, config)
            print_summary(paths.folder.name, results)
            print(f"[INFO] CSV saved to: {paths.output_all_candidates_csv}")
            print(f"[INFO] TXT saved to: {paths.output_validated_txt}")
            print(f"[INFO] XYZ saved to: {paths.output_pihole_xyz}")
            print(f"[INFO] Plots saved to: {paths.output_plot_dir}")
        except Exception as exc:
            if not config.continue_on_error:
                raise
            failures.append((target, exc))
            print(f"[ERROR] Analysis failed for {target.folder}: {exc}")

    if failures:
        failed = ", ".join(str(target.folder) for target, _exc in failures)
        raise RuntimeError(f"Analysis failed for {len(failures)} target(s): {failed}")
