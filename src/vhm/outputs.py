from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path
from typing import List, Optional

import numpy as np

from .models import CandidateResult


def save_candidates_csv(results: List[CandidateResult], out_csv: Path) -> None:
    """Save all candidate descriptors to a CSV file."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if not results:
        with out_csv.open("w", newline="", encoding="utf-8") as handle:
            handle.write("no_candidates\n")
        return

    fieldnames = list(asdict(results[0]).keys())
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(asdict(result))


def format_optional(value: Optional[float], fmt: str = ".6f") -> str:
    """Format optional numerical values for text reports."""
    if value is None:
        return "None"
    if isinstance(value, float) and np.isnan(value):
        return "nan"
    return format(value, fmt)


def save_validated_txt(results: List[CandidateResult], out_txt: Path) -> None:
    """Save a human-readable TXT report with validated pi-holes."""
    validated = [result for result in results if result.validated]
    out_txt.parent.mkdir(parents=True, exist_ok=True)

    with out_txt.open("w", encoding="utf-8") as handle:
        handle.write("VALIDATED PI-HOLES\n")
        handle.write("=" * 80 + "\n\n")

        if not validated:
            handle.write("No validated pi-hole found.\n")
            return

        for counter, result in enumerate(sorted(validated, key=lambda item: item.vs_max, reverse=True), start=1):
            handle.write(f"Candidate {counter}\n")
            handle.write("-" * 80 + "\n")
            handle.write(f"Side: {result.side}\n")
            handle.write(f"Point index: {result.point_index}\n")
            handle.write(f"Local coordinates (Angstrom): ({result.x_local:.6f}, {result.y_local:.6f}, {result.z_local:.6f})\n")
            handle.write(f"Global coordinates (Angstrom): ({result.x_global:.6f}, {result.y_global:.6f}, {result.z_global:.6f})\n")
            handle.write(f"V_s,max (kcal/mol): {result.vs_max:.6f}\n")
            handle.write(f"Number of MLS neighbors: {result.n_mls_neighbors}\n")
            handle.write(f"|grad V| 3D (kcal mol^-1 Angstrom^-1): {format_optional(result.grad_norm)}\n")
            handle.write(f"grad V 3D = ({format_optional(result.grad_x)}, {format_optional(result.grad_y)}, {format_optional(result.grad_z)})\n")
            handle.write(f"3D Hessian eigenvalues: ({format_optional(result.hess_eig1)}, {format_optional(result.hess_eig2)}, {format_optional(result.hess_eig3)})\n")
            handle.write(f"Shift to 3D stationary point (Angstrom): {format_optional(result.stationary_shift)}\n")
            handle.write(f"|grad V| tangential (kcal mol^-1 Angstrom^-1): {format_optional(result.tangential_grad_norm)}\n")
            handle.write(f"Tangential grad V = ({format_optional(result.tangential_grad_u)}, {format_optional(result.tangential_grad_v)})\n")
            handle.write(f"Tangential Hessian eigenvalues: ({format_optional(result.tangential_hess_eig1)}, {format_optional(result.tangential_hess_eig2)})\n")
            handle.write(f"Shift to tangential stationary point (Angstrom): {format_optional(result.tangential_stationary_shift)}\n")
            handle.write(f"V_hole,mean (kcal/mol): {format_optional(result.vhole_mean)}\n")
            handle.write(f"V_hole,max (kcal/mol): {format_optional(result.vhole_max)}\n")
            handle.write(f"V_hole,min (kcal/mol): {format_optional(result.vhole_min)}\n")
            handle.write(f"Delta V_hole (kcal/mol): {format_optional(result.delta_vhole)}\n")
            handle.write(f"Anisotropy angle (degrees): {format_optional(result.anisotropy_angle_deg)}\n")
            if result.min_sector_id is not None:
                handle.write(f"Minimum sector: S{result.min_sector_id + 1}\n")
            if result.max_sector_id is not None:
                handle.write(f"Maximum sector: S{result.max_sector_id + 1}\n")
            handle.write(f"Note: {result.reason}\n\n")


def save_molecule_with_piholes_xyz(
    atom_symbols: List[str],
    atom_coords_global: np.ndarray,
    results: List[CandidateResult],
    out_xyz: Path,
    only_validated: bool = True,
) -> None:
    """Write an XYZ file containing the molecule plus pi-holes as dummy X atoms."""
    selected = [result for result in results if result.validated] if only_validated else results
    natoms_total = len(atom_symbols) + len(selected)

    out_xyz.parent.mkdir(parents=True, exist_ok=True)
    with out_xyz.open("w", encoding="utf-8") as handle:
        handle.write(f"{natoms_total}\n")
        handle.write("Molecule + pi-holes as dummy atoms X\n")

        for symbol, coord in zip(atom_symbols, atom_coords_global):
            handle.write(f"{symbol:2s}  {coord[0]: .8f}  {coord[1]: .8f}  {coord[2]: .8f}\n")

        for result in selected:
            handle.write(f"X   {result.x_global: .8f}  {result.y_global: .8f}  {result.z_global: .8f}\n")
