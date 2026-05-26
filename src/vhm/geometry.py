from __future__ import annotations

from typing import Iterable, Optional, Tuple

import numpy as np

from .models import LocalFrame


def normalize(vec: np.ndarray) -> np.ndarray:
    """Return a normalized vector."""
    norm = np.linalg.norm(vec)
    if norm < 1.0e-14:
        raise ValueError("Cannot normalize a near-zero vector.")
    return vec / norm


def fit_plane_from_points(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fit a least-squares plane and return its centroid and normal vector."""
    centroid = points.mean(axis=0)
    centered = points - centroid
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = normalize(vh[-1])
    return centroid, normal


def validate_atom_indices(atom_coords: np.ndarray, indices_1based: Iterable[int], label: str) -> None:
    """Validate that 1-based atom indices are inside the XYZ atom range."""
    natoms = len(atom_coords)
    bad = [idx for idx in indices_1based if idx < 1 or idx > natoms]
    if bad:
        raise IndexError(f"Invalid {label} atom indices {bad}; XYZ contains {natoms} atoms.")


def build_local_frame(
    atom_coords: np.ndarray,
    reference_indices_1based: Iterable[int],
    orientation_atom_index_1based: Optional[int],
) -> LocalFrame:
    """Build a right-handed local frame from four reference atoms."""
    reference_indices_1based = tuple(reference_indices_1based)
    if len(reference_indices_1based) != 4:
        raise ValueError("Exactly four reference atoms are required.")

    validate_atom_indices(atom_coords, reference_indices_1based, "reference")
    if orientation_atom_index_1based is not None:
        validate_atom_indices(atom_coords, (orientation_atom_index_1based,), "orientation")

    reference_idx0 = [idx - 1 for idx in reference_indices_1based]
    reference_points = atom_coords[reference_idx0]
    origin, z_axis = fit_plane_from_points(reference_points)

    if orientation_atom_index_1based is not None:
        orientation_point = atom_coords[orientation_atom_index_1based - 1]
        orientation_vec = orientation_point - origin
        if np.dot(orientation_vec, z_axis) < 0.0:
            z_axis = -z_axis

    first_ref_vec = atom_coords[reference_idx0[0]] - origin
    first_ref_proj = first_ref_vec - np.dot(first_ref_vec, z_axis) * z_axis
    x_axis = normalize(first_ref_proj)

    y_axis = normalize(np.cross(z_axis, x_axis))
    x_axis = normalize(np.cross(y_axis, z_axis))

    rotation_matrix = np.vstack([x_axis, y_axis, z_axis])
    return LocalFrame(origin=origin, rotation_matrix=rotation_matrix)


def build_bond_frame(
    atom_coords: np.ndarray,
    bond_indices_1based: Tuple[int, int],
) -> LocalFrame:
    """
    Build a right-handed local frame for sigma-hole analysis.

    For bond_indices_1based=(A, B), the origin is atom B and local +z points
    from atom A to atom B. Positive local z therefore follows the bond
    extension beyond atom B.
    """
    if len(bond_indices_1based) != 2:
        raise ValueError("Exactly two bond atoms are required.")

    validate_atom_indices(atom_coords, bond_indices_1based, "bond")
    atom_a_idx0, atom_b_idx0 = (idx - 1 for idx in bond_indices_1based)
    origin = atom_coords[atom_b_idx0]
    z_axis = normalize(atom_coords[atom_b_idx0] - atom_coords[atom_a_idx0])

    helper = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(helper, z_axis))) > 0.90:
        helper = np.array([0.0, 1.0, 0.0], dtype=float)

    x_axis = normalize(helper - np.dot(helper, z_axis) * z_axis)
    y_axis = normalize(np.cross(z_axis, x_axis))
    x_axis = normalize(np.cross(y_axis, z_axis))

    rotation_matrix = np.vstack([x_axis, y_axis, z_axis])
    return LocalFrame(origin=origin, rotation_matrix=rotation_matrix)


def to_local(coords: np.ndarray, frame: LocalFrame) -> np.ndarray:
    """Transform global coordinates to the local frame."""
    return (coords - frame.origin) @ frame.rotation_matrix.T


def to_global(local_coords: np.ndarray, frame: LocalFrame) -> np.ndarray:
    """Transform local coordinates to the global frame."""
    return local_coords @ frame.rotation_matrix + frame.origin
