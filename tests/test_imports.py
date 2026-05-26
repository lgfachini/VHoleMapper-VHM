def test_package_imports():
    import vhm

    assert vhm.__version__


def test_public_api_exports():
    from vhm import AnalysisConfig, SigmaTargetSpec, TargetSpec, run_batch

    assert AnalysisConfig is not None
    assert TargetSpec is not None
    assert SigmaTargetSpec is not None
    assert callable(run_batch)


def test_settings_loads():
    from importlib import import_module

    settings = import_module("config.settings")
    assert settings.RUN_MODE in {"auto", "manual"}
    assert settings.HOLE_TYPE in {"pi", "sigma"}
    assert hasattr(settings, "ANALYSIS_OPTIONS")
