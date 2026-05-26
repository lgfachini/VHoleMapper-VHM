from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .config import AnalysisConfig
from .geometry import normalize
from .models import MlsFit


def gaussian_weights(distances: np.ndarray, sigma: float) -> np.ndarray:
    """Return Gaussian weights for distance-weighted fitting."""
    sigma = max(float(sigma), 1.0e-10)
    return np.exp(-0.5 * (distances / sigma) ** 2)


def build_quadratic_design_matrix(dxyz: np.ndarray) -> np.ndarray:
    """Build the design matrix for a quadratic 3D polynomial."""
    dx = dxyz[:, 0]
    dy = dxyz[:, 1]
    dz = dxyz[:, 2]

    return np.column_stack(
        [
            np.ones(len(dxyz)),
            dx,
            dy,
            dz,
            0.5 * dx * dx,
            dx * dy,
            dx * dz,
            0.5 * dy * dy,
            dy * dz,
            0.5 * dz * dz,
        ]
    )


def weighted_least_squares(A: np.ndarray, y: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Solve a weighted least-squares problem."""
    sqrt_weights = np.sqrt(weights)
    Aw = A * sqrt_weights[:, None]
    yw = y * sqrt_weights
    coeffs, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    return coeffs


def weighted_pca_basis(dxyz: np.ndarray, weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute a weighted PCA tangent-normal basis."""
    weight_sum = np.sum(weights)
    if weight_sum <= 1.0e-16:
        raise ValueError("The weighted PCA weight sum is too small.")

    covariance = (dxyz * weights[:, None]).T @ dxyz / weight_sum
    eigvals, eigvecs = np.linalg.eigh(covariance)

    order = np.argsort(eigvals)
    normal = normalize(eigvecs[:, order[0]])
    t1 = normalize(eigvecs[:, order[2]])
    t2 = normalize(eigvecs[:, order[1]])

    if np.dot(np.cross(t1, t2), normal) < 0.0:
        t2 = -t2

    return t1, t2, normal


def project_to_surface_basis(
    dxyz: np.ndarray,
    t1: np.ndarray,
    t2: np.ndarray,
    normal: np.ndarray,
) -> np.ndarray:
    """Project Cartesian displacements onto the local tangent-normal basis."""
    u = dxyz @ t1
    v = dxyz @ t2
    w = dxyz @ normal
    return np.column_stack([u, v, w])


def mls_local_fit(
    local_coords: np.ndarray,
    potentials: np.ndarray,
    tree: cKDTree,
    point_index: int,
    config: AnalysisConfig,
) -> Optional[MlsFit]:
    """Run a local MLS fit around one candidate."""
    center = local_coords[point_index]
    neighbors = tree.query_ball_point(center, r=config.mls_radius)

    if not config.mls_include_candidate:
        neighbors = [idx for idx in neighbors if idx != point_index]

    if len(neighbors) < config.mls_min_neighbors:
        return None

    neighbor_indices = np.array(neighbors, dtype=int)
    dxyz = local_coords[neighbor_indices] - center
    y = potentials[neighbor_indices] - potentials[point_index]

    distances = np.linalg.norm(dxyz, axis=1)
    weights = gaussian_weights(distances, config.mls_sigma)

    try:
        t1, t2, normal = weighted_pca_basis(dxyz, weights)
    except Exception:
        return None

    duvw = project_to_surface_basis(dxyz, t1, t2, normal)
    design_matrix = build_quadratic_design_matrix(duvw)

    try:
        coeffs = weighted_least_squares(design_matrix, y, weights)
    except np.linalg.LinAlgError:
        return None

    grad_uvw = np.array([coeffs[1], coeffs[2], coeffs[3]], dtype=float)
    hessian_uvw = np.array(
        [
            [coeffs[4], coeffs[5], coeffs[6]],
            [coeffs[5], coeffs[7], coeffs[8]],
            [coeffs[6], coeffs[8], coeffs[9]],
        ],
        dtype=float,
    )

    try:
        eigvals_3d = np.linalg.eigvalsh(hessian_uvw)
    except np.linalg.LinAlgError:
        return None

    grad_tan = grad_uvw[:2].copy()
    hessian_tan = hessian_uvw[:2, :2].copy()

    try:
        eigvals_tan = np.linalg.eigvalsh(hessian_tan)
    except np.linalg.LinAlgError:
        return None

    try:
        tangential_stationary_delta = -np.linalg.solve(hessian_tan, grad_tan)
        tangential_stationary_shift = float(np.linalg.norm(tangential_stationary_delta))
    except np.linalg.LinAlgError:
        tangential_stationary_shift = None

    grad_xyz = grad_uvw[0] * t1 + grad_uvw[1] * t2 + grad_uvw[2] * normal
    grad_norm = float(np.linalg.norm(grad_xyz))

    try:
        stationary_delta_uvw = -np.linalg.solve(hessian_uvw, grad_uvw)
        stationary_shift = float(np.linalg.norm(stationary_delta_uvw))
    except np.linalg.LinAlgError:
        stationary_shift = None

    return MlsFit(
        t1=t1,
        t2=t2,
        normal=normal,
        grad_uvw=grad_uvw,
        grad_xyz=grad_xyz,
        grad_norm=grad_norm,
        hessian_uvw=hessian_uvw,
        eigvals_3d=eigvals_3d,
        stationary_shift=stationary_shift,
        grad_tan=grad_tan,
        tangential_grad_norm=float(np.linalg.norm(grad_tan)),
        hessian_tan=hessian_tan,
        eigvals_tan=eigvals_tan,
        tangential_stationary_shift=tangential_stationary_shift,
        n_neighbors=len(neighbor_indices),
    )


def mls_candidate_ok(fit: Optional[MlsFit], config: AnalysisConfig) -> Tuple[bool, str]:
    """Check whether an MLS fit satisfies the tangential derivative criteria."""
    if fit is None:
        return False, "MLS fit unavailable or insufficient local neighborhood."

    if fit.tangential_grad_norm > config.tangential_grad_norm_thr:
        return False, (
            "Tangential gradient is too large "
            f"({fit.tangential_grad_norm:.4f} > {config.tangential_grad_norm_thr:.4f})."
        )

    max_tan_eig = float(np.max(fit.eigvals_tan))

    if config.allow_semidefinite_tangential_maximum:
        if max_tan_eig > 1.0e-8:
            return False, (
                "Tangential Hessian is not negative semidefinite "
                f"(lambda_max,tan = {max_tan_eig:.4f})."
            )
    else:
        if max_tan_eig >= config.tangential_hess_eig_max_thr:
            return False, (
                "Tangential Hessian is not sufficiently negative "
                f"(lambda_max,tan = {max_tan_eig:.4f}; "
                f"threshold = {config.tangential_hess_eig_max_thr:.4f})."
            )

    if config.use_tangential_stationary_shift_test:
        if fit.tangential_stationary_shift is None:
            return False, "Could not estimate the tangential stationary point."
        if fit.tangential_stationary_shift > config.max_tangential_stationary_shift:
            return False, (
                "Tangential stationary point is too far "
                f"({fit.tangential_stationary_shift:.4f} > "
                f"{config.max_tangential_stationary_shift:.4f} Angstrom)."
            )

    return True, "Tangential MLS test passed."
