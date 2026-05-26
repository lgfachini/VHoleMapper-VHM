from pathlib import Path
from types import SimpleNamespace

import pytest

from vhm import AnalysisConfig, SigmaTargetSpec, TargetSpec
from vhm.cli import selected_targets


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


def test_sigma_target_rejects_same_atom_twice():
    with pytest.raises(ValueError):
        SigmaTargetSpec(
            folder=Path("data/benzene"),
            bond_atom_indices_1based=(1, 1),
        )


def test_sigma_config_requires_sigma_targets():
    with pytest.raises(ValueError):
        AnalysisConfig(targets=(target(),), hole_type="sigma")


def test_sigma_config_accepts_sigma_target():
    sigma_target = SigmaTargetSpec(
        folder=Path("data/benzene"),
        bond_atom_indices_1based=(1, 2),
    )
    config = AnalysisConfig(targets=(sigma_target,), hole_type="sigma")

    assert config.hole_type == "sigma"


def test_manual_targets_resolve_relative_folders_from_project_root():
    settings = SimpleNamespace(
        RUN_MODE=" Manual ",
        HOLE_TYPE="pi",
        MANUAL_TARGETS=(target(),),
    )

    selected = selected_targets(settings)

    assert selected[0].folder.is_absolute()
    assert selected[0].folder.name == "benzene"


def test_manual_sigma_targets_resolve_relative_folders_from_project_root():
    settings = SimpleNamespace(
        RUN_MODE="manual",
        HOLE_TYPE="sigma",
        MANUAL_TARGETS=(
            SigmaTargetSpec(
                folder=Path("data/benzene"),
                bond_atom_indices_1based=(1, 2),
            ),
        ),
    )

    selected = selected_targets(settings)

    assert selected[0].folder.is_absolute()
    assert selected[0].folder.name == "benzene"
