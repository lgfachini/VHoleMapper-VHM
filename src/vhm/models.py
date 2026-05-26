from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class LocalFrame:
    """Local coordinate frame built from the user-defined reference atoms."""

    origin: np.ndarray
    rotation_matrix: np.ndarray


@dataclass
class MlsFit:
    """Results from the local MLS fit."""

    t1: np.ndarray
    t2: np.ndarray
    normal: np.ndarray
    grad_uvw: np.ndarray
    grad_xyz: np.ndarray
    grad_norm: float
    hessian_uvw: np.ndarray
    eigvals_3d: np.ndarray
    stationary_shift: Optional[float]
    grad_tan: np.ndarray
    tangential_grad_norm: float
    hessian_tan: np.ndarray
    eigvals_tan: np.ndarray
    tangential_stationary_shift: Optional[float]
    n_neighbors: int


@dataclass
class CandidateResult:
    """Computed descriptors for one pi-hole candidate."""

    side: str
    point_index: int
    x_local: float
    y_local: float
    z_local: float
    x_global: float
    y_global: float
    z_global: float
    vs_max: float

    grad_norm: Optional[float]
    grad_x: Optional[float]
    grad_y: Optional[float]
    grad_z: Optional[float]
    hess_eig1: Optional[float]
    hess_eig2: Optional[float]
    hess_eig3: Optional[float]
    stationary_shift: Optional[float]

    tangential_grad_norm: Optional[float]
    tangential_grad_u: Optional[float]
    tangential_grad_v: Optional[float]
    tangential_hess_eig1: Optional[float]
    tangential_hess_eig2: Optional[float]
    tangential_stationary_shift: Optional[float]

    n_mls_neighbors: Optional[int]

    vhole_mean: Optional[float]
    vhole_max: Optional[float]
    vhole_min: Optional[float]
    delta_vhole: Optional[float]
    mean_sector_potential: Optional[float]
    min_sector_potential: Optional[float]
    max_sector_potential: Optional[float]
    min_sector_id: Optional[int]
    max_sector_id: Optional[int]
    anisotropy_angle_deg: Optional[float]

    validated: bool
    reason: str
