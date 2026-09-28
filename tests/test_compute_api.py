"""Tests that the calculations work without matplotlib.

DataBrowser draws its QC views with pyqtgraph, so it must be able to use these
calculations without matplotlib. The figures still need matplotlib.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

SRC = str(Path(__file__).resolve().parent.parent / "src")


def _run_without_matplotlib(code: str) -> subprocess.CompletedProcess:
    """Runs code in a new interpreter in which importing matplotlib fails."""
    guard = "import sys; sys.modules['matplotlib'] = None\n"
    env = {**os.environ, "PYTHONPATH": SRC}
    return subprocess.run([sys.executable, "-c", guard + code], capture_output=True,
                          text=True, env=env)


def test_calculations_import_without_matplotlib():
    result = _run_without_matplotlib(
        "import qualitymetrics\n"
        "from qualitymetrics.ksdata import KilosortResults, unit_summary_table\n"
        "from qualitymetrics.metrics import noise_cutoff_parts\n"
        "import qualitymetrics.raw\n"
        "print('ok')\n")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "ok"


def test_figures_still_need_matplotlib():
    """Checks that the guard works, so the test above can't pass by accident."""
    result = _run_without_matplotlib("import qualitymetrics.plots\n")
    assert result.returncode != 0


def test_use_lab_style_is_still_importable_from_the_package():
    from qualitymetrics import use_lab_style
    from qualitymetrics.style import use_lab_style as direct

    assert use_lab_style is direct


def test_unit_summary_table_keeps_its_old_import_path():
    from qualitymetrics import plots
    from qualitymetrics.ksdata import unit_summary_table

    assert plots.unit_summary_table is unit_summary_table
