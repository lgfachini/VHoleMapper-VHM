from pathlib import Path

import pytest

from vhm import AnalysisConfig, TargetSpec


def target() -> TargetSpec:
    return TargetSpec(
        folder=Path("data/benzene"),
        reference_atom_indices_1based=(1, 2, 3, 4),
    )


def test_target_requires_four_reference_atoms():
    with pytest.raises(ValueError):
        TargetSpec(
            folder=Path("data/benzene"),
            reference_atom_indices_1based=(1, 2, 3),
        )


def test_target_rejects_zero_based_indices():
    with pytest.raises(ValueError):
        TargetSpec(
            folder=Path("data/benzene"),
            reference_atom_indices_1based=(0, 1, 2, 3),
        )


def test_config_rejects_invalid_units():
    with pytest.raises(ValueError):
        AnalysisConfig(targets=(target(),), xyz_coord_unit="nanometers")


def test_config_rejects_inverted_search_window():
    with pytest.raises(ValueError):
        AnalysisConfig(targets=(target(),), z_search_min=3.0, z_search_max=1.0)


def test_config_rejects_invalid_sector_count():
    with pytest.raises(ValueError):
        AnalysisConfig(targets=(target(),), angular_sectors=1)
