from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from .config import AnalysisConfig
from .geometry import to_global
from .models import CandidateResult, LocalFrame, MlsFit


def wrap_angle_deg(angle: float) -> float:
    """Wrap an angle to the [0, 360) interval."""
    return angle % 360.0


def smallest_angular_difference_deg(a: float, b: float) -> float:
    """Return the smallest angular difference between two angles in degrees."""
    diff = abs((a - b) % 360.0)
    return min(diff, 360.0 - diff)


def nan_mls_values() -> Tuple[np.ndarray, float, np.ndarray, Optional[float], np.ndarray, float, np.ndarray, Optional[float], Optional[int]]:
    """Return placeholder values for candidates without an MLS fit."""
    return (
        np.array([np.nan, np.nan, np.nan], dtype=float),
        float("nan"),
        np.array([np.nan, np.nan, np.nan], dtype=float),
        None,
        np.array([np.nan, np.nan], dtype=float),
        float("nan"),
        np.array([np.nan, np.nan], dtype=float),
        None,
        None,
    )


def extract_mls_values(fit: Optional[MlsFit]):
    """Extract scalar/vector MLS descriptors or return NaN placeholders."""
    if fit is None:
        return nan_mls_values()

    return (
        fit.grad_xyz,
        fit.grad_norm,
        fit.eigvals_3d,
        fit.stationary_shift,
        fit.grad_tan,
        fit.tangential_grad_norm,
        fit.eigvals_tan,
        fit.tangential_stationary_shift,
        fit.n_neighbors,
    )


def make_candidate_result(
    side: str,
    point_index: int,
    cand_local: np.ndarray,
    cand_global: np.ndarray,
    v_cand: float,
    mls_fit: Optional[MlsFit],
    validated: bool,
    reason: str,
    sector_values: Optional[dict] = None,
) -> CandidateResult:
    """Create a CandidateResult object with optional sector descriptors."""
    (
        grad_xyz,
        grad_norm,
        eigvals_3d,
        stationary_shift,
        grad_tan,
        tangential_grad_norm,
        eigvals_tan,
        tangential_stationary_shift,
        n_neighbors,
    ) = extract_mls_values(mls_fit)

    values = sector_values or {}

    return CandidateResult(
        side=side,
        point_index=point_index,
        x_local=float(cand_local[0]),
        y_local=float(cand_local[1]),
        z_local=float(cand_local[2]),
        x_global=float(cand_global[0]),
        y_global=float(cand_global[1]),
        z_global=float(cand_global[2]),
        vs_max=float(v_cand),
        grad_norm=float(grad_norm),
        grad_x=float(grad_xyz[0]),
        grad_y=float(grad_xyz[1]),
        grad_z=float(grad_xyz[2]),
        hess_eig1=float(eigvals_3d[0]),
        hess_eig2=float(eigvals_3d[1]),
        hess_eig3=float(eigvals_3d[2]),
        stationary_shift=None if stationary_shift is None else float(stationary_shift),
        tangential_grad_norm=float(tangential_grad_norm),
        tangential_grad_u=float(grad_tan[0]),
        tangential_grad_v=float(grad_tan[1]),
        tangential_hess_eig1=float(eigvals_tan[0]),
        tangential_hess_eig2=float(eigvals_tan[1]),
        tangential_stationary_shift=(None if tangential_stationary_shift is None else float(tangential_stationary_shift)),
        n_mls_neighbors=n_neighbors,
        vhole_mean=values.get("vhole_mean"),
        vhole_max=values.get("vhole_max"),
        vhole_min=values.get("vhole_min"),
        delta_vhole=values.get("delta_vhole"),
        mean_sector_potential=values.get("mean_sector_potential"),
        min_sector_potential=values.get("min_sector_potential"),
        max_sector_potential=values.get("max_sector_potential"),
        min_sector_id=values.get("min_sector_id"),
        max_sector_id=values.get("max_sector_id"),
        anisotropy_angle_deg=values.get("anisotropy_angle_deg"),
        validated=validated,
        reason=reason,
    )


def compute_sector_means(
    local_coords: np.ndarray,
    potentials: np.ndarray,
    candidate_local: np.ndarray,
    side: str,
    config: AnalysisConfig,
) -> np.ndarray:
    """Compute mean potential in each angular sector around a candidate."""
    dx = local_coords[:, 0] - candidate_local[0]
    dy = local_coords[:, 1] - candidate_local[1]
    dz = local_coords[:, 2] - candidate_local[2]

    radial_xy = np.sqrt(dx ** 2 + dy ** 2)
    angles = np.degrees(np.arctan2(dy, dx))
    angles = np.array([wrap_angle_deg(angle) for angle in angles])

    shell_mask = (
        (radial_xy >= config.sector_inner_radius)
        & (radial_xy <= config.sector_outer_radius)
        & (np.abs(dz) <= config.sector_half_height)
    )

    if side == "above":
        shell_mask &= local_coords[:, 2] > 0.0
    elif side == "below":
        shell_mask &= local_coords[:, 2] < 0.0
    else:
        raise ValueError("side must be either 'above' or 'below'.")

    shell_indices = np.where(shell_mask)[0]
    if len(shell_indices) == 0:
        return np.array([], dtype=float)

    sector_width = 360.0 / config.angular_sectors
    sector_means: List[float] = []

    for sector_id in range(config.angular_sectors):
        angle_start = sector_id * sector_width
        angle_end = (sector_id + 1) * sector_width

        idx = shell_indices[(angles[shell_indices] >= angle_start) & (angles[shell_indices] < angle_end)]

        if len(idx) < config.min_points_per_sector:
            sector_means.append(float("nan"))
        else:
            sector_means.append(float(np.mean(potentials[idx])))

    return np.array(sector_means, dtype=float)


def sector_validation(
    local_coords: np.ndarray,
    potentials: np.ndarray,
    point_index: int,
    side: str,
    frame: LocalFrame,
    config: AnalysisConfig,
    mls_fit: Optional[MlsFit] = None,
    pre_reason: str = "",
) -> CandidateResult:
    """Validate a candidate by comparing it to sector-averaged surrounding potentials."""
    candidate_local = local_coords[point_index]
    candidate_global = to_global(candidate_local.reshape(1, 3), frame)[0]
    v_cand = float(potentials[point_index])

    sector_means = compute_sector_means(local_coords, potentials, candidate_local, side, config)

    if len(sector_means) == 0:
        return make_candidate_result(
            side, point_index, candidate_local, candidate_global, v_cand, mls_fit,
            False, f"{pre_reason} No points found in the sector shell.".strip()
        )

    valid_mask = ~np.isnan(sector_means)
    if np.sum(valid_mask) < 2:
        return make_candidate_result(
            side, point_index, candidate_local, candidate_global, v_cand, mls_fit,
            False, f"{pre_reason} Insufficient number of sectors with data.".strip()
        )

    valid_sector_means = sector_means[valid_mask]
    mean_sector_potential = float(np.mean(valid_sector_means))
    min_sector_potential = float(np.min(valid_sector_means))
    max_sector_potential = float(np.max(valid_sector_means))

    min_sector_id = int(np.where(sector_means == min_sector_potential)[0][0])
    max_sector_id = int(np.where(sector_means == max_sector_potential)[0][0])

    vhole_mean = float(v_cand - mean_sector_potential)
    vhole_max = float(v_cand - min_sector_potential)
    vhole_min = float(v_cand - max_sector_potential)
    delta_vhole = float(max_sector_potential - min_sector_potential)

    sector_width = 360.0 / config.angular_sectors
    min_center_angle = (min_sector_id + 0.5) * sector_width
    max_center_angle = (max_sector_id + 0.5) * sector_width
    anisotropy_angle = float(smallest_angular_difference_deg(min_center_angle, max_center_angle))

    validated = vhole_mean >= config.delta_thr
    reason = "Validated." if validated else f"V_hole,mean < DELTA_THR ({config.delta_thr:.3f})."
    if pre_reason:
        reason = f"{pre_reason} {reason}"

    sector_values = {
        "vhole_mean": vhole_mean,
        "vhole_max": vhole_max,
        "vhole_min": vhole_min,
        "delta_vhole": delta_vhole,
        "mean_sector_potential": mean_sector_potential,
        "min_sector_potential": min_sector_potential,
        "max_sector_potential": max_sector_potential,
        "min_sector_id": min_sector_id,
        "max_sector_id": max_sector_id,
        "anisotropy_angle_deg": anisotropy_angle,
    }

    return make_candidate_result(
        side, point_index, candidate_local, candidate_global, v_cand, mls_fit,
        validated, reason, sector_values
    )
