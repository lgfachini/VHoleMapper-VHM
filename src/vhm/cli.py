from __future__ import annotations

import argparse
import sys
from importlib import import_module
from pathlib import Path
from typing import Iterable, Optional, Tuple

from . import __version__
from .config import AnalysisConfig, TargetSpec


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_settings():
    root = project_root()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return import_module("config.settings")


def resolve_project_path(path: Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return project_root() / path


def discover_folders(data_dir: Path, pattern: str = "*") -> Tuple[Path, ...]:
    """Return direct subfolders inside data_dir matching a glob pattern."""
    data_dir = resolve_project_path(data_dir)
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")
    if not data_dir.is_dir():
        raise NotADirectoryError(f"DATA_DIR is not a directory: {data_dir}")

    folders = tuple(sorted(path for path in data_dir.glob(pattern) if path.is_dir()))
    if not folders:
        raise RuntimeError(f"No folders found in {data_dir} with pattern {pattern!r}")
    return folders


def build_auto_targets(
    folders: Iterable[Path],
    reference_atoms: Tuple[int, int, int, int],
    orientation_atom: Optional[int],
    xyz_filename: Optional[str],
    vtx_filename: str,
) -> Tuple[TargetSpec, ...]:
    """Create TargetSpec entries for folders discovered in auto mode."""
    return tuple(
        TargetSpec(
            folder=folder,
            reference_atom_indices_1based=reference_atoms,
            orientation_atom_index_1based=orientation_atom,
            xyz_filename=xyz_filename,
            vtx_filename=vtx_filename,
        )
        for folder in folders
    )


def selected_targets(settings) -> Tuple[TargetSpec, ...]:
    """Resolve analysis targets from the user-editable settings module."""
    run_mode = getattr(settings, "RUN_MODE", "auto")
    if run_mode == "auto":
        folders = discover_folders(settings.DATA_DIR, getattr(settings, "AUTO_FOLDER_GLOB", "*"))
        return build_auto_targets(
            folders=folders,
            reference_atoms=settings.AUTO_REFERENCE_ATOMS,
            orientation_atom=getattr(settings, "AUTO_ORIENTATION_ATOM", None),
            xyz_filename=getattr(settings, "AUTO_XYZ_FILENAME", None),
            vtx_filename=getattr(settings, "AUTO_VTX_FILENAME", "vtx.txt"),
        )

    if run_mode == "manual":
        manual_targets = tuple(getattr(settings, "MANUAL_TARGETS", ()))
        if not manual_targets:
            raise RuntimeError("RUN_MODE is 'manual', but MANUAL_TARGETS is empty.")
        return manual_targets

    raise ValueError("RUN_MODE must be either 'auto' or 'manual'.")


def config_from_settings(settings) -> AnalysisConfig:
    kwargs = dict(getattr(settings, "ANALYSIS_OPTIONS", {}))
    return AnalysisConfig(targets=selected_targets(settings), **kwargs)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vhm",
        description="VHoleMapper: detect, validate, and analyze pi-holes from ESP surface data.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    try:
        parser.parse_args(argv)
        settings = load_settings()
        config = config_from_settings(settings)
        from .pipeline import run_batch

        print("Targets selected for analysis:")
        for target in config.targets:
            print(
                f"  - {Path(target.folder)} | "
                f"reference atoms = {target.reference_atom_indices_1based} | "
                f"orientation atom = {target.orientation_atom_index_1based}"
            )

        run_batch(config)
    except (OSError, RuntimeError, ValueError, IndexError) as exc:
        parser.exit(2, f"{parser.prog}: error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
