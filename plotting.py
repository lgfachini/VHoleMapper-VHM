from __future__ import annotations

import math
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Wedge

from config import AnalysisConfig
from geometry import to_local
from models import CandidateResult, LocalFrame


def plot_reference_atoms(
    ax: plt.Axes,
    atom_symbols: List[str],
    atom_coords_global: np.ndarray,
    frame: LocalFrame,
    reference_atom_indices_1based: Tuple[int, int, int, int],
    config: AnalysisConfig,
) -> None:
    """Plot the four reference atoms projected onto the local reference plane."""
    reference_idx0 = [idx - 1 for idx in reference_atom_indices_1based]
    reference_symbols = [atom_symbols[idx] for idx in reference_idx0]
    reference_coords_global = atom_coords_global[reference_idx0]
    reference_coords_local = to_local(reference_coords_global, frame)

    ax.scatter(reference_coords_local[:, 0], reference_coords_local[:, 1], s=350, marker="o", color="green")

    if config.plot_atom_labels:
        for idx_1based, symbol, coord in zip(reference_atom_indices_1based, reference_symbols, reference_coords_local):
            ax.text(coord[0], coord[1], f"{symbol}{idx_1based}", fontsize=10, ha="center", va="center", color="black")


def plot_sector_geometry(ax: plt.Axes, result: CandidateResult, config: AnalysisConfig) -> None:
    """Plot sector rings, sector dividers, labels, and highlighted min/max sectors."""
    px = result.x_local
    py = result.y_local

    ring_radius = 0.5 * (config.sector_inner_radius + config.sector_outer_radius) if config.plot_sector_ring_radius is None else config.plot_sector_ring_radius
    sector_width = 360.0 / config.angular_sectors

    ax.add_patch(Circle((px, py), config.sector_inner_radius, fill=False, linestyle="--", linewidth=1.0))
    ax.add_patch(Circle((px, py), config.sector_outer_radius, fill=False, linestyle="--", linewidth=1.0))

    for sector_id in range(config.angular_sectors):
        theta1 = sector_id * sector_width
        theta2 = (sector_id + 1) * sector_width

        if sector_id in {result.min_sector_id, result.max_sector_id}:
            ax.add_patch(
                Wedge(
                    center=(px, py),
                    r=config.sector_outer_radius,
                    theta1=theta1,
                    theta2=theta2,
                    width=config.sector_outer_radius - config.sector_inner_radius,
                    alpha=0.25,
                )
            )

        angle = math.radians(theta1)
        x1 = px + config.sector_inner_radius * math.cos(angle)
        y1 = py + config.sector_inner_radius * math.sin(angle)
        x2 = px + config.sector_outer_radius * math.cos(angle)
        y2 = py + config.sector_outer_radius * math.sin(angle)
        ax.plot([x1, x2], [y1, y2], linewidth=0.8, color="black")

        label_angle = math.radians(theta1 + 0.5 * sector_width)
        lx = px + (ring_radius + 0.15) * math.cos(label_angle)
        ly = py + (ring_radius + 0.15) * math.sin(label_angle)
        ax.text(lx, ly, f"S{sector_id + 1}", fontsize=8, ha="center", va="center")

    angle = math.radians(360.0)
    x1 = px + config.sector_inner_radius * math.cos(angle)
    y1 = py + config.sector_inner_radius * math.sin(angle)
    x2 = px + config.sector_outer_radius * math.cos(angle)
    y2 = py + config.sector_outer_radius * math.sin(angle)
    ax.plot([x1, x2], [y1, y2], linewidth=0.8, color="black")


def plot_anisotropy_angle(ax: plt.Axes, result: CandidateResult, config: AnalysisConfig) -> None:
    """Plot anisotropy vectors and Delta V_angle annotation."""
    if result.min_sector_id is None or result.max_sector_id is None:
        return

    px = result.x_local
    py = result.y_local

    ring_radius = 0.5 * (config.sector_inner_radius + config.sector_outer_radius) if config.plot_sector_ring_radius is None else config.plot_sector_ring_radius
    sector_width = 360.0 / config.angular_sectors

    min_angle_deg = (result.min_sector_id + 0.5) * sector_width
    max_angle_deg = (result.max_sector_id + 0.5) * sector_width

    min_angle = math.radians(min_angle_deg)
    max_angle = math.radians(max_angle_deg)

    min_x = px + ring_radius * math.cos(min_angle)
    min_y = py + ring_radius * math.sin(min_angle)
    max_x = px + ring_radius * math.cos(max_angle)
    max_y = py + ring_radius * math.sin(max_angle)

    ax.plot([px, min_x], [py, min_y], linewidth=2.0)
    ax.plot([px, max_x], [py, max_y], linewidth=2.0)

    theta_start = min(min_angle_deg, max_angle_deg)
    theta_end = max(min_angle_deg, max_angle_deg)
    if theta_end - theta_start > 180.0:
        theta_start, theta_end = theta_end, theta_start + 360.0

    ax.add_patch(Wedge(center=(px, py), r=0.65 * config.sector_inner_radius, theta1=theta_start, theta2=theta_end, width=0.03, alpha=0.5))

    if result.anisotropy_angle_deg is not None:
        ax.text(px - 0.15, py - 0.45, f"$\\Delta V_{{angle}} = {result.anisotropy_angle_deg:.1f}^\\circ$", fontsize=8, ha="left", va="top")


def plot_shell_points(ax: plt.Axes, result: CandidateResult, local_coords_surface: np.ndarray, config: AnalysisConfig) -> None:
    """Optionally plot all shell points used in sector validation."""
    if not config.plot_show_all_shell_points:
        return

    dx = local_coords_surface[:, 0] - result.x_local
    dy = local_coords_surface[:, 1] - result.y_local
    dz = local_coords_surface[:, 2] - result.z_local
    radial_xy = np.sqrt(dx ** 2 + dy ** 2)

    shell_mask = (
        (radial_xy >= config.sector_inner_radius)
        & (radial_xy <= config.sector_outer_radius)
        & (np.abs(dz) <= config.sector_half_height)
    )

    if result.side == "above":
        shell_mask &= local_coords_surface[:, 2] > 0.0
    else:
        shell_mask &= local_coords_surface[:, 2] < 0.0

    shell_points = local_coords_surface[shell_mask]
    if len(shell_points) > 0:
        ax.scatter(shell_points[:, 0], shell_points[:, 1], s=10, alpha=0.5)


def plot_pihole_planar_view(
    atom_symbols: List[str],
    atom_coords_global: np.ndarray,
    frame: LocalFrame,
    reference_atom_indices_1based: Tuple[int, int, int, int],
    result: CandidateResult,
    out_svg: Path,
    local_coords_surface: np.ndarray,
    config: AnalysisConfig,
) -> None:
    """Create a 2D plot in the reference plane."""
    fig, ax = plt.subplots(figsize=config.plot_figsize)

    plot_reference_atoms(ax, atom_symbols, atom_coords_global, frame, reference_atom_indices_1based, config)

    ax.scatter([result.x_local], [result.y_local], s=180, marker="x")
    ax.text(result.x_local, result.y_local, "    pi-hole", fontsize=10, va="bottom", ha="left")

    plot_sector_geometry(ax, result, config)
    plot_anisotropy_angle(ax, result, config)
    plot_shell_points(ax, result, local_coords_surface, config)

    delta_for_title = result.delta_vhole if result.delta_vhole is not None else float("nan")
    ax.set_xlabel("x_local (Å)")
    ax.set_ylabel("y_local (Å)")
    ax.set_title(
        f"{result.side} | point {result.point_index} | "
        f"Vs,max = {result.vs_max:.2f} kcal/mol | "
        f"Delta Vhole = {delta_for_title:.2f}"
    )
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)

    out_svg.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_svg, format="svg")
    plt.close(fig)


def save_pihole_plots(
    atom_symbols: List[str],
    atom_coords_global: np.ndarray,
    frame: LocalFrame,
    reference_atom_indices_1based: Tuple[int, int, int, int],
    local_coords_surface: np.ndarray,
    results: List[CandidateResult],
    out_dir: Path,
    config: AnalysisConfig,
) -> None:
    """Save SVG plots for selected candidates."""
    out_dir.mkdir(parents=True, exist_ok=True)
    selected = [result for result in results if result.validated] if config.plot_only_validated else results

    for counter, result in enumerate(selected, start=1):
        filename = f"pihole_{counter:02d}_{result.side}_idx{result.point_index}.svg"
        plot_pihole_planar_view(
            atom_symbols=atom_symbols,
            atom_coords_global=atom_coords_global,
            frame=frame,
            reference_atom_indices_1based=reference_atom_indices_1based,
            result=result,
            out_svg=out_dir / filename,
            local_coords_surface=local_coords_surface,
            config=config,
        )
