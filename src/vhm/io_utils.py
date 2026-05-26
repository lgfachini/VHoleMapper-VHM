from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import numpy as np


def convert_coords_unit(coords: np.ndarray, unit: str, bohr_to_ang: float) -> np.ndarray:
    """Convert coordinates to angstrom."""
    unit = unit.lower().strip()
    if unit == "angstrom":
        return coords.copy()
    if unit == "bohr":
        return coords * bohr_to_ang
    raise ValueError(f"Unknown coordinate unit: {unit}")


def read_xyz(xyz_file: Path, coord_unit: str, bohr_to_ang: float) -> Tuple[List[str], np.ndarray]:
    """Read an XYZ file and return atom symbols and coordinates in angstrom."""
    if not xyz_file.exists():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")

    with xyz_file.open("r", encoding="utf-8") as handle:
        lines = [line.rstrip() for line in handle if line.strip()]

    try:
        natoms = int(lines[0])
    except Exception as exc:
        raise ValueError("The first XYZ line must contain the number of atoms.") from exc

    atom_lines = lines[2:2 + natoms]
    if len(atom_lines) != natoms:
        raise ValueError("The number of XYZ atom lines is inconsistent with the header.")

    symbols: List[str] = []
    coords: List[List[float]] = []

    for line in atom_lines:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ line: {line}")
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    coords_array = np.array(coords, dtype=float)
    return symbols, convert_coords_unit(coords_array, coord_unit, bohr_to_ang)


def read_vtx(vtx_file: Path, coord_unit: str, bohr_to_ang: float) -> Tuple[np.ndarray, np.ndarray]:
    """Read a VTX surface file and return coordinates in angstrom plus potentials."""
    if not vtx_file.exists():
        raise FileNotFoundError(f"VTX file not found: {vtx_file}")

    with vtx_file.open("r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]

    try:
        npoints = int(lines[0])
    except Exception as exc:
        raise ValueError("The first VTX line must contain the number of points.") from exc

    data_lines = lines[1:]
    if len(data_lines) != npoints:
        raise ValueError(
            f"Inconsistent number of VTX points. Expected {npoints}, found {len(data_lines)}."
        )

    coords: List[List[float]] = []
    potentials: List[float] = []

    for line in data_lines:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid VTX line: {line}")
        x, y, z = map(float, parts[:3])
        pot_kcal = float(parts[-1])
        coords.append([x, y, z])
        potentials.append(pot_kcal)

    coords_array = np.array(coords, dtype=float)
    potentials_array = np.array(potentials, dtype=float)
    return convert_coords_unit(coords_array, coord_unit, bohr_to_ang), potentials_array
