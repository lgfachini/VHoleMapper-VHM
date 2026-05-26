from __future__ import annotations

from typing import List

import numpy as np
from scipy.spatial import cKDTree

from .config import AnalysisConfig


def get_side_mask(local_coords: np.ndarray, side: str, config: AnalysisConfig) -> np.ndarray:
    """Return a boolean mask for points in the axial search region on one side."""
    z = local_coords[:, 2]
    radial = np.sqrt(local_coords[:, 0] ** 2 + local_coords[:, 1] ** 2)

    if side == "above":
        return (
            (z >= config.z_search_min)
            & (z <= config.z_search_max)
            & (radial <= config.max_distance_from_axis)
        )

    if side == "below":
        return (
            (z <= -config.z_search_min)
            & (z >= -config.z_search_max)
            & (radial <= config.max_distance_from_axis)
        )

    raise ValueError("side must be either 'above' or 'below'.")


def find_local_maxima(
    local_coords: np.ndarray,
    potentials: np.ndarray,
    candidate_indices: np.ndarray,
    config: AnalysisConfig,
) -> List[int]:
    """Find discrete local maxima inside a preselected candidate region."""
    if len(candidate_indices) == 0:
        return []

    sub_coords = local_coords[candidate_indices]
    sub_potentials = potentials[candidate_indices]

    tree = cKDTree(sub_coords)
    maxima_local_indices: List[int] = []

    for local_i, point in enumerate(sub_coords):
        neighbors = tree.query_ball_point(point, r=config.local_max_radius)
        neighbors = [idx for idx in neighbors if idx != local_i]

        if len(neighbors) < config.min_local_neighbors:
            continue

        my_potential = sub_potentials[local_i]
        neighbor_potentials = sub_potentials[neighbors]

        if config.require_positive_potential and my_potential <= 0.0:
            continue

        if config.strict_local_max:
            is_maximum = bool(np.all(my_potential > neighbor_potentials))
        else:
            is_maximum = bool(np.all(my_potential >= neighbor_potentials))

        if is_maximum:
            maxima_local_indices.append(local_i)

    maxima_global_indices = [int(candidate_indices[idx]) for idx in maxima_local_indices]

    if (
        config.top_n_candidates_per_side is not None
        and len(maxima_global_indices) > config.top_n_candidates_per_side
    ):
        maxima_global_indices = sorted(
            maxima_global_indices,
            key=lambda idx: potentials[idx],
            reverse=True,
        )[: config.top_n_candidates_per_side]

    return maxima_global_indices
