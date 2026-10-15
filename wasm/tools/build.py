"""Build a single-threaded, float64 LAMMPS/RDKit browser module."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / ".local" / "wasm"


def environment():
    env = os.environ.copy()
    sdk = LOCAL / "emsdk"
    env["EMSDK"] = str(sdk)
    env["EM_CONFIG"] = str(sdk / ".emscripten")
    env["EMSDK_PYTHON"] = sys.executable
    env["EM_CACHE"] = str(sdk / "upstream" / "emscripten" / "cache")
    env["EMCC_CORES"] = "4"
    tools = LOCAL / "buildtools"
    env["PATH"] = os.pathsep.join([str(tools / "cmake" / "data" / "bin"),
        str(tools / "ninja" / "data" / "bin"), str(tools / "bin"),
        str(sdk / "upstream" / "bin"), str(sdk / "upstream" / "emscripten"), env["PATH"]])
    return env


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--configure-only", action="store_true")
    parser.add_argument("--link-only", action="store_true", help="Reuse already-built dependency libraries")
    args = parser.parse_args()
    env = environment()
    cmake = LOCAL / "buildtools" / "cmake" / "data" / "bin" / "cmake.exe"
    ninja = LOCAL / "buildtools" / "ninja" / "data" / "bin" / "ninja.exe"
    if not ninja.is_file():
        ninja = LOCAL / "buildtools" / "bin" / "ninja.exe"
    sdk = LOCAL / "emsdk"
    subprocess.run([sys.executable, str(sdk / "upstream" / "emscripten" / "embuilder.py"), "build", "zlib"], env=env, check=True)
    toolchain = sdk / "upstream" / "emscripten" / "cmake" / "Modules" / "Platform" / "Emscripten.cmake"
    common = ["-G", "Ninja", f"-DCMAKE_MAKE_PROGRAM={ninja.as_posix()}",
        f"-DCMAKE_TOOLCHAIN_FILE={toolchain.as_posix()}", "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_CXX_STANDARD=17", "-DCMAKE_CXX_FLAGS=-O3 -fexceptions -ffp-contract=off",
        "-DCMAKE_C_FLAGS=-O3 -fexceptions -ffp-contract=off",
        f"-DPython3_EXECUTABLE={Path(sys.executable).as_posix()}",
        f"-DPython_EXECUTABLE={Path(sys.executable).as_posix()}",
        f"-DPYTHON_EXECUTABLE={Path(sys.executable).as_posix()}"]
    configurations = [
        ("lammps", LOCAL / "lammps" / "cmake", ["-DBUILD_MPI=OFF", "-DBUILD_OMP=OFF",
            "-DBUILD_SHARED_LIBS=OFF", "-DLAMMPS_EXCEPTIONS=ON", "-DLAMMPS_SIZES=smallbig",
            "-DBUILD_TOOLS=OFF", "-DBUILD_DOC=OFF"], ["lammps"]),
        ("rdkit", LOCAL / "rdkit", [
            "-DRDK_BUILD_PYTHON_WRAPPERS=OFF", "-DRDK_BUILD_CPP_TESTS=OFF",
            "-DRDK_BUILD_MINIMAL_LIB=ON", "-DRDK_BUILD_THREADSAFE_SSS=OFF",
            "-DRDK_USE_BOOST_SERIALIZATION=OFF", "-DRDK_USE_BOOST_IOSTREAMS=OFF",
            "-DRDK_USE_BOOST_STACKTRACE=OFF", "-DRDK_BUILD_DESCRIPTORS3D=OFF",
            "-DRDK_BUILD_FREETYPE_SUPPORT=OFF", "-DRDK_INSTALL_COMIC_FONTS=OFF",
            "-DRDK_BUILD_COORDGEN_SUPPORT=OFF", "-DRDK_BUILD_MAEPARSER_SUPPORT=OFF",
            "-DRDK_BUILD_CHEMDRAW_SUPPORT=OFF", "-DRDK_BUILD_MOLINTERCHANGE_SUPPORT=OFF",
            "-DRDK_BUILD_PUBCHEMSHAPE_SUPPORT=OFF", "-DRDK_USE_URF=OFF",
            "-DRDK_OPTIMIZE_POPCNT=OFF", "-DRDK_TEST_MMFF_COMPLIANCE=OFF",
            f"-DBoost_DIR={(ROOT / 'wasm' / 'cmake' / 'boost').as_posix()}",
            f"-DFETCHCONTENT_SOURCE_DIR_CATCH2={(LOCAL / 'catch2').as_posix()}",
            f"-DFETCHCONTENT_SOURCE_DIR_BETTER_ENUMS={(LOCAL / 'better-enums').as_posix()}",
            f"-DMOLECULAR_CAR_BOOST_DIR={(LOCAL / 'boost').as_posix()}"],
            ["FileParsers", "ForceFieldHelpers", "PartialCharges", "DistGeomHelpers"]),
    ]
    for name, source, flags, targets in configurations:
        if args.link_only:
            continue
        build = LOCAL / f"{name}-build"
        print(f"Configuring {name}", flush=True)
        subprocess.run([str(cmake), "-S", str(source), "-B", str(build), *common, *flags], env=env, check=True)
        if not args.configure_only:
            subprocess.run([str(cmake), "--build", str(build), "--target", *targets, "--parallel", str(args.jobs)], env=env, check=True)
    if args.configure_only:
        return
    output = ROOT / "wasm" / "dist" / "browser" / "runtime"
    output.mkdir(parents=True, exist_ok=True)
    empp = sdk / "upstream" / "emscripten" / "em++.py"
    libraries = list((LOCAL / "rdkit" / "lib").glob("*.a"))
    libraries += list((LOCAL / "rdkit-build" / "lib").glob("*.a"))
    if not libraries:
        raise RuntimeError("Missing RDKit static libraries")
    command = [sys.executable, str(empp), "-O3", "-std=c++20", "-fexceptions", "-ffp-contract=off", "-sUSE_ZLIB=1",
        "-DLAMMPS_SMALLBIG", "-DRDK_BUILD_STATIC_LIBS", "-DBOOST_DISABLE_THREADS",
        "-I", str(LOCAL / "lammps" / "src"), "-I", str(LOCAL / "rdkit" / "Code"),
        "-I", str(LOCAL / "rdkit-build" / "Code"), "-I", str(LOCAL / "boost"),
        str(ROOT / "wasm" / "src" / "potentials.cpp"), str(ROOT / "wasm" / "src" / "engine.cpp"),
        str(LOCAL / "lammps-build" / "liblammps.a"),
        "-Wl,--start-group", *map(str, libraries), "-Wl,--end-group", "--bind", "--no-entry",
        "-sMODULARIZE=1", "-sEXPORT_ES6=1", "-sENVIRONMENT=web,worker,node",
        "-sALLOW_MEMORY_GROWTH=1", "-sINITIAL_MEMORY=67108864", "-sMAXIMUM_MEMORY=536870912",
        "-sDISABLE_EXCEPTION_CATCHING=0", "-sEXPORTED_RUNTIME_METHODS=getExceptionMessage,decrementExceptionRefcount",
        "-sEXPORT_NAME=createMolecularCar", "-o", str(output / "molecular-car.mjs")]
    sources = sorted((ROOT / "wasm" / "src").glob("*"))
    source_hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources if p.is_file()}
    print("Linking combined browser MD module", flush=True)
    subprocess.run(command, env=env, check=True)
    for name, digest in source_hashes.items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"Source changed during compilation; rebuild before using this module: {name}")
    pins = json.loads((ROOT / "wasm" / "dependencies.json").read_text())
    sources = sorted((ROOT / "wasm" / "src").glob("*"))
    (output / "build.json").write_text(json.dumps({"dependencies": pins,
        "source_archives": json.loads((ROOT / "wasm" / "source-checksums.json").read_text()),
        "own_sources": source_hashes,
        "precision": "float64", "threads": 1, "fast_math": False,
        "python_runtime": False, "reference_model": "native UFF / frozen Au bonded-saip-game-v1 / fixed-gradient-race-v1",
        "validation": "Run tools/validate.mjs; build completion alone does not confirm numerical agreement."}, indent=2) + "\n")
    print(f"Built {output}", flush=True)


if __name__ == "__main__":
    main()
