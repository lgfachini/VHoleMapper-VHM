"""VHoleMapper: pi-hole and sigma-hole detection and analysis tools."""

from .config import AnalysisConfig, AnalysisTarget, MoleculePaths, SigmaTargetSpec, TargetSpec, make_paths

__version__ = "0.1.0"

__all__ = [
    "AnalysisConfig",
    "AnalysisTarget",
    "CandidateResult",
    "LocalFrame",
    "MlsFit",
    "MoleculePaths",
    "SigmaTargetSpec",
    "TargetSpec",
    "__version__",
    "analyze_molecule",
    "make_paths",
    "run_batch",
    "write_outputs",
]


def __getattr__(name: str):
    if name in {"CandidateResult", "LocalFrame", "MlsFit"}:
        from . import models

        return getattr(models, name)
    if name in {"analyze_molecule", "run_batch", "write_outputs"}:
        from . import pipeline

        return getattr(pipeline, name)
    raise AttributeError(f"module 'vhm' has no attribute {name!r}")
