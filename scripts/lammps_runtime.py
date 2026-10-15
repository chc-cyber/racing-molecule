"""Load one matching Windows LAMMPS installation in a worker process."""

from contextlib import contextmanager
import importlib
import json
import os
from pathlib import Path
import struct
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]


@contextmanager
def open_lammps(installation=None, expected_version=None):
    if sys.version_info < (3, 10) or struct.calcsize("P") != 8:
        raise RuntimeError("Use a 64-bit Python 3.10 or newer interpreter.")
    if os.name != "nt":
        raise RuntimeError("This initial integration supports Windows installations.")

    config_path = PROJECT_ROOT / ".local" / "lammps.json"
    config = json.loads(config_path.read_text(encoding="utf-8-sig")) if config_path.exists() else {}
    selected = installation or os.environ.get("MOLECULAR_CAR_LAMMPS_HOME") or config.get("installation")
    if not selected:
        raise RuntimeError("Run scripts/setup.ps1 or provide --lammps-home.")
    root = Path(selected).resolve()
    module_directory = root / "Python"
    binary_directory = root / "bin"
    library_path = binary_directory / "liblammps.dll"
    for required in (module_directory / "lammps" / "__init__.py", library_path):
        if not required.is_file():
            raise FileNotFoundError(f"Missing LAMMPS component: {required}")

    # These settings apply only to this Python process, never the user's registry.
    os.environ["LAMMPSDLLPATH"] = str(binary_directory)
    os.environ["LAMMPS_POTENTIALS"] = str(root / "Potentials")
    os.environ["LAMMPS_PLUGIN_PATH"] = ""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.path.insert(0, str(module_directory))
    directory_handle = os.add_dll_directory(str(binary_directory))
    simulation = None
    try:
        module = importlib.import_module("lammps")
        actual_module = Path(module.__file__).resolve()
        if not actual_module.is_relative_to(module_directory):
            raise RuntimeError(f"Another LAMMPS module was already loaded: {actual_module}")
        simulation = module.lammps(cmdargs=["-log", "none", "-screen", "none", "-nocite"])
        version = simulation.version()
        # Explicit selections can have their own pin. Reuse the local pin only
        # when the selected installation actually matches that configuration.
        expected = expected_version
        if expected is None and config.get("installation"):
            if root == Path(config["installation"]).resolve():
                expected = config.get("expected_version")
        if expected and version != expected:
            raise RuntimeError(f"Expected LAMMPS {expected}, loaded {version}.")
        yield simulation, {
            "python_version": sys.version.split()[0],
            "python_executable": sys.executable,
            "python_bits": struct.calcsize("P") * 8,
            "lammps_version": version,
            "lammps_module": str(actual_module),
            "lammps_library": str(library_path),
        }
    finally:
        if simulation is not None:
            simulation.close()
        directory_handle.close()
        sys.path.remove(str(module_directory))
