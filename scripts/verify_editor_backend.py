"""Run editor backend acceptance checks and save a reviewable JSON report.

Usage: .\.venv\Scripts\python.exe -E -B scripts\verify_editor_backend.py
No HTTP server, game controller, or LAMMPS worker is started or stopped.
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
import io
import json
from pathlib import Path
import platform
import sys
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from rdkit import rdBase


class RecordedResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cases = []

    def addSuccess(self, test):
        super().addSuccess(test)
        self.cases.append({"name": test.id(), "passed": True})

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.cases.append({"name": test.id(), "passed": False, "error": self._exc_info_to_string(err, test)})

    def addError(self, test, err):
        super().addError(test, err)
        self.cases.append({"name": test.id(), "passed": False, "error": self._exc_info_to_string(err, test)})


def main():
    from tests.test_editor_service import EditorServiceTests
    start = time.perf_counter()
    output = io.StringIO()
    runner = unittest.TextTestRunner(stream=output, verbosity=2, resultclass=RecordedResult)
    result = runner.run(unittest.defaultTestLoader.loadTestsFromTestCase(EditorServiceTests))
    report = {"format": "molecular-car-editor-backend-check-v1",
              "checkedAt": datetime.now(timezone(timedelta(hours=9))).isoformat(),
              "passed": result.wasSuccessful(), "testCount": result.testsRun,
              "seconds": round(time.perf_counter() - start, 3),
              "versions": {"python": platform.python_version(), "rdkit": rdBase.rdkitVersion, "numpy": np.__version__},
              "scope": "RDKit file/geometry service, stable indexed atom identity, fixed-point UFF, stereo roundtrip, fragments and actual race preparation on a private molecule. Does not prove browser UI or LAMMPS loading.",
              "checks": result.cases, "log": output.getvalue()}
    destination = ROOT / "outputs" / "editor-check" / "backend-report.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(output.getvalue())
    print(f"Report: {destination}")
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
