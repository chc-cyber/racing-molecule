"""Real LAMMPS integration; RDKit supplies UFF forces through fix external.

All coordinates are Angstrom, time is femtoseconds, energy is kcal/mol.
The explicit surface/driver is a game model, not a validated Au/STM model.
"""

import time
import numpy as np

from scripts.lammps_runtime import open_lammps
from .molecule import prepare_molecule
from .external_fields import ExternalFieldEvaluator
from .relaxation import relax_joint_geometry, JointRelaxationError
from .track import TRACK
from .start_pose import relaxation_heading_report
from .element_models import SUPPORTED_RACE_ELEMENTS, element_support_metadata
from .fields import (energy_and_forces, GRADIENT_BASE_FIELD, GRADIENT_FIELD_SLOPE,
                     GRADIENT_POLARIZATION_GAIN, ADSORPTION_DEPTH,
                     ADSORPTION_HEIGHT, ADSORPTION_ALPHA, SURFACE_SPACING, gradient_coefficients)

DEFAULT_SETTINGS = {"temperature": 150.0, "damping_fs": 350.0, "corrugation": 0.08, "lattice_scale": 1.5,
                    "surface": "fourier-au"}
SURFACE_MODES = ("continuum", "periodic-au", "fourier-au")
FROZEN_AU_MODES = ("periodic-au", "fourier-au")
SURFACE_LABELS = {"periodic-au": "冻结 Au 周期势", "fourier-au": "冻结 Au 代理势 v2"}
DEFAULT_TIP = {"x": 6.0, "y": 0.0, "z": 7.0, "strength": 0.6}
DEFAULT_DRIVE = {"mode": "gradient", "strength": 0.6, "growth": 1.}
TIMESTEP_FS = 1.0
FRAME_STEPS = 50  # Preserve the 50 fs animation sampling interval.
MAX_STEPS = 5_000_000  # 5 ns: also allow useful runs at weaker drive strengths.
DEFAULT_RANDOM_SEED = 98245165
WORLD_BOUNDS = ((-110.0, 1110.0), (-350.0, 350.0), (-40.0, 120.0))


def validate_surface_mode(mode):
    if not isinstance(mode, str) or mode not in SURFACE_MODES:
        raise ValueError("表面模型必须是 continuum、periodic-au 或 fourier-au。")
    return mode


def periodic_surface_error_message(error):
    message = str(error)
    if "more than three normal neighbours" in message:
        return "局部近邻超过 3 个，超出当前 SAIP 法向模型的支持范围。" + message
    if "zero normal magnitude" in message:
        return "局部平面退化，近邻可能共线，无法确定 SAIP 法向。" + message
    return message


def drive_and_walls_energy_and_forces(positions, charges, alphas, tip, drive):
    """Keep the existing field and track walls, removing only the Morse surface.

    Corrugation is explicitly zero here. Subtracting the known Morse term from
    the original field function preserves both drive modes and the soft walls;
    subtracting a zero-drive call would accidentally remove those walls too.
    Coordinates remain absolute: neither the race nor its driver is periodic.
    """
    energy, forces = energy_and_forces(positions, charges, alphas, tip, 0.0, drive=drive)
    height = np.asarray(positions)[:, 2] - ADSORPTION_HEIGHT
    exponential = np.exp(-ADSORPTION_ALPHA * height)
    energy -= float(np.sum(ADSORPTION_DEPTH * (exponential**2 - 2.0 * exponential)))
    forces[:, 2] -= 2.0 * ADSORPTION_DEPTH * ADSORPTION_ALPHA * (exponential**2 - exponential)
    return energy, forces


class Simulation:
    def __init__(self, content, filename, settings=None, *, thermostat=True,
                 random_seed=DEFAULT_RANDOM_SEED, drive=None, external_accelerator="auto", start_pose=None):
        self.settings = dict(DEFAULT_SETTINGS, **(settings or {}))
        from .lattice_corrugation import validate_lattice_scale
        self.settings["lattice_scale"] = validate_lattice_scale(self.settings.get("lattice_scale", DEFAULT_SETTINGS["lattice_scale"]))
        validate_surface_mode(self.settings["surface"])
        self.surface = None
        if isinstance(random_seed, bool) or not isinstance(random_seed, int) or not 1 <= random_seed <= 900_000_000:
            raise ValueError("随机种子必须是 1–900000000 之间的整数。")
        self.random_seed = random_seed
        self.velocity_seed = 17329 if random_seed == DEFAULT_RANDOM_SEED else 1 + ((random_seed ^ 17329) % 900_000_000)
        self.drive = dict(DEFAULT_DRIVE)
        if drive:
            self.drive.update(drive)
        if self.drive["mode"] != "gradient" or not 0 <= self.drive["strength"] <= 1:
            raise ValueError("驱动模式或强度不正确。")
        self.callback_error = None
        self.sequence = 0
        self.prepared = prepare_molecule(content, filename, start_pose=start_pose)
        self.masses = np.array([a["mass"] for a in self.prepared.atoms])
        self.charges = np.array([a["charge"] for a in self.prepared.atoms])
        self.alphas = np.array([a["polarizability"] for a in self.prepared.atoms])
        self.n = len(self.masses)
        if np.ptp(self.prepared.positions, axis=0).max() > 45:
            raise ValueError("分子跨度超过当前赛车模型支持的 45 Å，请使用更小的模型。")
        if self.settings["surface"] in FROZEN_AU_MODES:
            mode = self.settings["surface"]
            label = SURFACE_LABELS[mode]
            unsupported = sorted({atom["element"] for atom in self.prepared.atoms} - set(SUPPORTED_RACE_ELEMENTS))
            if unsupported:
                raise ValueError(label + "缺少当前结构中这些元素的表面参数：" + "、".join(unsupported) + "。")
            if mode == "fourier-au":
                from .fourier_atlas import FourierAtlasPotential as create_surface
            else:
                from .periodic_au import create_periodic_surface as create_surface
            try:
                self.surface = create_surface(self.prepared.atoms, self.prepared.positions,
                    bonds=self.prepared.bonds, normal_policy="bonded-saip-game-v1")
                from .lattice_corrugation import LatticeCorrugation
                self.surface = LatticeCorrugation(self.surface, self.settings["lattice_scale"])
            except ValueError as error:
                raise ValueError(label + "无法准备当前结构：" + periodic_surface_error_message(error) +
                                 " 可选择连续表面继续使用该分子。") from error
        self._prepare_force_evaluation(external_accelerator)
        self.context = open_lammps()
        self.lmp, self.runtime = self.context.__enter__()
        self.runtime["forcefield"] = "RDKit UFF / LAMMPS integration"
        self.runtime["start_pose"] = self.prepared.start_pose
        self.runtime["elements"] = element_support_metadata(
            (a["element"] for a in self.prepared.atoms), self.settings["surface"])
        self.runtime["molecular_force_validation"] = getattr(self.prepared, "forcefield_validation", {})
        self.runtime["external_field_kernel"] = self.field_evaluator.metadata()
        self.runtime.update(model="uff-gradient-v3", timestep_fs=TIMESTEP_FS,
                            mass_policy="standard-element-masses-isotopes-ignored-v1",
                            langevin_seed=self.random_seed, velocity_seed=self.velocity_seed,
                            gradient={"base_field": GRADIENT_BASE_FIELD,
                                      "field_slope": GRADIENT_FIELD_SLOPE,
                                      "growth": self._growth,
                                      "actual_base_and_slope": list(gradient_coefficients(self._drive_strength, self._growth)),
                                      "polarization_game_gain": GRADIENT_POLARIZATION_GAIN})
        if self.surface is None:
            self.runtime["surface"] = {"mode": "continuum", "model": "effective-morse-triangular-v1",
                "energy_units": "kcal/mol", "coordinate_units": "angstrom",
                "adsorption_depth": ADSORPTION_DEPTH, "adsorption_height": ADSORPTION_HEIGHT,
                "adsorption_alpha": ADSORPTION_ALPHA, "lattice_spacing_angstrom": SURFACE_SPACING}
        else:
            mode = self.settings["surface"]
            self.runtime["model"] = f"uff-{mode}-bonded-normal-gradient-v2"
            self.runtime["surface"] = dict(self.surface.metadata(), mode=mode)
            # Per-evaluation proxy counters belong to frames, not immutable
            # preparation metadata copied to the HTTP process and replay.
            self.runtime["surface"].pop("last_evaluation_counts", None)
            self.runtime["forcefield"] = f"RDKit UFF / frozen Au {mode} surface / LAMMPS integration"
            if any(a["element"] not in ("C", "H") for a in self.prepared.atoms):
                self.runtime["model"] = f"uff-{mode}-hetero-lj-bonded-normal-gradient-v1"
                self.prepared.notes.append("C/H 使用原金表面 SAIP 游戏适配；其他元素使用非反应的 UFF 范德华混合金晶格势，尚未标定真实吸附或金属配位。")
        self.drive_enabled = False
        self.thermostat = thermostat
        try:
            self._initialize()
            if getattr(self.surface, "last_evaluation_counts", None) is not None:
                self.runtime["surface"]["preparation_evaluation_counts"] = dict(self.surface.last_evaluation_counts)
        except BaseException:
            self.close()
            raise

    @property
    def molecule(self):
        return {"name": self.prepared.name, "atoms": self.prepared.atoms,
                "bonds": self.prepared.bonds, "notes": self.prepared.notes}

    def _initialize(self):
        lmp = self.lmp
        # Never let the integration masses disagree with preparation/COM masses,
        # even if a future preparation path supplies distinct masses per element.
        atom_types = [(a["element"], float(a["mass"])) for a in self.prepared.atoms]
        types = list(dict.fromkeys(atom_types))
        lmp.commands_list([
            "units real", "atom_style atomic", "boundary f f f",
            "atom_modify map array sort 0 0.0",
            "region world block " + " ".join(str(v) for bounds in WORLD_BOUNDS for v in bounds) + " units box",
            f"create_box {len(types)} world",
        ])
        for i, (_, mass) in enumerate(types, 1):
            lmp.command(f"mass {i} {mass:.17g}")
        count = lmp.create_atoms(self.n, list(range(1, self.n + 1)),
                                 [types.index(atom_type) + 1 for atom_type in atom_types],
                                 self.prepared.positions.ravel().tolist())
        if count != self.n:
            raise RuntimeError("LAMMPS 未能创建全部赛车原子。")
        self.runtime.update(lammps_atoms=int(lmp.get_natoms()), integrated_atoms=self.n,
                            mobile_surface_atoms=0)
        lmp.commands_list([
            "pair_style zero 10.0", "pair_coeff * *",
            "fix molecular all external pf/callback 1 1",
            "compute internal_temperature all temp/com",
            "thermo_style custom step temp pe ke etotal",
            "thermo 1000000", "timestep 0.5",  # Keep preparation at its validated step.
        ])
        lmp.set_fix_external_callback("molecular", self._forces)
        try:
            self.runtime["preparation"] = relax_joint_geometry(
                lmp, self.n, self._forces, self._check_callback)
        except JointRelaxationError as error:
            self.runtime["preparation"] = error.diagnostics
            raise
        self.runtime["start_pose"]["after_relaxation"] = relaxation_heading_report(
            self.coordinates(), self.prepared.start_pose)
        lmp.command("reset_timestep 0")
        lmp.command(f"timestep {TIMESTEP_FS}")
        if self.settings["temperature"] > 0:
            lmp.command(f"velocity all create {self.settings['temperature']} {self.velocity_seed} mom yes rot yes dist gaussian")
        else:
            lmp.command("velocity all set 0 0 0")
        lmp.command("fix integrate all nve")
        if self.thermostat:
            self._set_bath()
        self.drive_enabled = True
        # A new run setup refreshes external forces after minimization/settings.
        lmp.command("run 0")
        self._check_callback()
        self.initial_center = self.center()
        race_length = TRACK["finish_x"] - TRACK["start_x"]
        self.runtime["race"] = {"rule": "longitudinal-displacement-v1",
                                "start_x": float(self.initial_center[0]),
                                "finish_x": float(self.initial_center[0]) + race_length,
                                "length_angstrom": race_length}
        self._cache_field_controls()
        self.started_at = time.perf_counter()

    def _forces(self, caller, step, nlocal, tags, positions, output):
        # ctypes callbacks cannot propagate exceptions through the native stack.
        # Save the exception and abort at the next batch boundary instead.
        try:
            if nlocal != self.n:
                raise RuntimeError("首版模拟需要单进程 LAMMPS。")
            # Also support isolated callback checks constructed without __init__.
            if not hasattr(self, "field_evaluator"):
                self._prepare_force_evaluation()
            # Existing local scripts also edit the public control dictionaries.
            # Revalidate only when their scalar values change, without rebuilding
            # effective dictionaries or checking fixed response arrays each step.
            if self._field_control_signature() != self._cached_field_signature:
                self._cache_field_controls()
            if self._last_tags is None or not np.array_equal(tags, self._last_tags):
                ids = np.asarray(tags, dtype=np.int64)-1
                if not np.array_equal(np.sort(ids), self._expected_ids):
                    raise RuntimeError("模拟原子标识不完整或重复，无法保持分子拓扑。")
                self._ids = ids
                self._last_tags = np.array(tags, copy=True)
                self._identity_order = np.array_equal(ids, self._expected_ids)
            ids = self._ids
            if self._identity_order:
                ordered = positions
            else:
                ordered = self._ordered
                ordered[ids] = positions
            flat = ordered.ravel().tolist()
            molecular_energy = self.prepared.forcefield.CalcEnergy(flat)
            forces = self._molecular_forces
            np.negative(np.asarray(self.prepared.forcefield.CalcGrad(flat)).reshape(self.n, 3), out=forces)
            strength = self._drive_strength if self.drive_enabled else 0.
            field_energy, field_forces = self.field_evaluator.energy_and_forces(
                ordered, 1, strength, 0., 0.,
                0., self._corrugation, self._growth)
            if self.surface is not None:
                surface_energy, surface_forces = self.surface.energy_and_forces(ordered)
                field_energy += surface_energy
                forces += surface_forces
            energy = float(molecular_energy + field_energy)
            forces += field_forces
            if not np.isfinite(energy) or not np.isfinite(forces).all():
                raise RuntimeError("模拟出现非有限受力，请降低驱动强度或整理输入结构。")
            if self._identity_order:
                np.copyto(output, forces)
            else:
                np.take(forces, ids, axis=0, out=output)
            self.lmp.fix_external_set_energy_global("molecular", energy)
        except Exception as error:
            self.callback_error = (periodic_surface_error_message(error)
                                   if self.surface is not None else str(error))
            output[:] = 0

    def _prepare_force_evaluation(self, accelerator="auto"):
        self.field_evaluator = ExternalFieldEvaluator(self.charges, self.alphas,
            surface_enabled=self.surface is None, accelerator=accelerator)
        self._ordered = np.empty((self.n, 3), dtype=np.float64)
        self._molecular_forces = np.empty((self.n, 3), dtype=np.float64)
        self._expected_ids = np.arange(self.n, dtype=np.int64)
        self._last_tags = None
        self._ids = self._expected_ids
        self._identity_order = True
        self._cache_field_controls()

    def _field_control_signature(self):
        mode = self.drive["mode"]
        return (mode, self.drive["strength"], self.drive.get("growth", 1.), self.settings["corrugation"])

    def _cache_field_controls(self):
        mode, strength = self.drive["mode"], self.drive["strength"]
        if (mode != "gradient" or isinstance(strength, bool)
                or not np.isfinite(strength) or not 0 <= strength <= 1):
            raise ValueError("驱动模式或强度不正确。")
        self._drive_mode = 1 if mode == "gradient" else 0
        self._drive_strength = float(strength)
        self._growth = self.drive.get("growth", 1.)
        gradient_coefficients(self._drive_strength, self._growth)
        self._growth = float(self._growth)
        self._corrugation = float(self.settings["corrugation"])
        if not np.isfinite(self._corrugation):
            raise ValueError("表面起伏必须是有限数值。")
        self._cached_field_signature = self._field_control_signature()

    def _check_callback(self):
        if self.callback_error:
            raise RuntimeError(self.callback_error)

    def coordinates(self):
        raw = self.lmp.gather_atoms("x", 1, 3)
        return np.ctypeslib.as_array(raw, shape=(self.n * 3,)).copy().reshape(self.n, 3)

    def center(self):
        return np.average(self.coordinates(), axis=0, weights=self.masses)

    def _set_bath(self):
        self.lmp.command(f"fix bath all langevin {self.settings['temperature']} {self.settings['temperature']} {self.settings['damping_fs']} {self.random_seed} tally yes")

    def update(self, changes):
        if "surface" in changes and changes["surface"] != self.settings["surface"]:
            raise ValueError("更换表面模型需要重置并重新准备分子。")
        if "corrugation" in changes and self.surface is not None:
            raise ValueError("冻结 Au 表面势的起伏由金晶格决定；请在连续表面模式调整起伏。")
        if "lattice_scale" in changes:
            from .lattice_corrugation import validate_lattice_scale
            scale = validate_lattice_scale(changes["lattice_scale"])
            if self.surface is None and scale != 1.:
                raise ValueError("晶格起伏倍率仅适用于冻结 Au 表面。")
            self.settings["lattice_scale"] = scale
            if self.surface is not None:
                self.surface.scale = scale
        if "tip" in changes or "autopilot" in changes:
            raise ValueError("针尖控制已移除。")
        if "drive" in changes:
            candidate = dict(self.drive, **changes["drive"])
            if candidate["mode"] != "gradient":
                raise ValueError("仅支持固定广域梯度场。")
            self.drive = candidate
        reset_bath = any(key in changes for key in ("temperature", "damping_fs"))
        for key in DEFAULT_SETTINGS:
            if key in changes:
                self.settings[key] = changes[key]
        if hasattr(self, "field_evaluator"):
            self._cache_field_controls()
        if reset_bath and self.thermostat:
            self.lmp.command("unfix bath")
            self._set_bath()

    def effective_drive(self, step):
        drive = dict(self.drive)
        if not self.drive_enabled:
            drive["strength"] = 0.0
        return drive

    def advance(self, steps=FRAME_STEPS):
        self.lmp.command(f"run {steps} post no")
        self._check_callback()
        self.sequence += 1
        return self.frame()

    def frame(self):
        positions = self.coordinates()
        raw_v = self.lmp.gather_atoms("v", 1, 3)
        velocities = np.ctypeslib.as_array(raw_v, shape=(self.n * 3,)).reshape(self.n, 3)
        center = np.average(positions, axis=0, weights=self.masses)
        velocity = np.average(velocities, axis=0, weights=self.masses)
        step = int(self.lmp.get_thermo("step"))
        if not np.isfinite(positions).all() or not np.isfinite(velocities).all():
            raise RuntimeError("模拟出现无效坐标或速度。")
        if (any(positions[:, axis].min() <= low or positions[:, axis].max() >= high
                for axis, (low, high) in enumerate(WORLD_BOUNDS))
                or positions[:, 2].min() < -5 or positions[:, 2].max() > 100):
            raise RuntimeError("赛车离开了原型的有效模拟区域，请重置并降低场强。")
        x = float(center[0])
        distance = x - float(self.initial_center[0])
        race_length = TRACK["finish_x"] - TRACK["start_x"]
        return {"seq": self.sequence, "step": step, "time_fs": step * TIMESTEP_FS,
                "positions": positions.ravel().tolist(), "center": center.tolist(),
                "speed_angstrom_ps": float(np.linalg.norm(velocity) * 1000),
                "energy": float(self.lmp.get_thermo("etotal")),
                "internal_temperature": float(self.lmp.extract_compute("internal_temperature", 0, 0)),
                "surface_mode": self.settings["surface"],
                "surface_evaluation": (dict(self.surface.last_evaluation_counts)
                                       if getattr(self.surface, "last_evaluation_counts", None) is not None else None),
                "drive": dict(self.drive),
                "effective_drive": self.effective_drive(step),
                "progress": min(1.0, max(0.0, distance / race_length)),
                "distance_angstrom": distance,
                "race_start_x": float(self.initial_center[0]),
                "race_length_angstrom": race_length}

    def close(self):
        if getattr(self, "context", None) is not None:
            self.context.__exit__(None, None, None)
            self.context = None


def worker_main(content, filename, settings, commands, events, initial_drive=None, start_pose=None):
    simulation = None
    try:
        events.put({"kind": "status", "status": "preparing"})
        simulation = Simulation(content, filename, settings, drive=initial_drive, start_pose=start_pose)
        events.put({"kind": "ready", "molecule": simulation.molecule,
                    "runtime": simulation.runtime, "frame": simulation.frame(),
                    "control": {"drive": dict(simulation.drive)},
                    "settings": simulation.settings})
        running = False
        playback_limit_fs = float("inf")
        computed_time_fs = 0.0
        while True:
            if not running or computed_time_fs >= playback_limit_fs:
                command = commands.get()
                pending = [command]
            else:
                pending = []
            while not commands.empty():
                try:
                    pending.append(commands.get_nowait())
                except Exception:
                    break
            for command in pending:
                action = command.get("action")
                if action == "stop":
                    return
                if action == "start":
                    running = True
                    events.put({"kind": "status", "status": "running"})
                elif action == "playback":
                    playback_limit_fs = command["until_fs"]
                elif action == "pause":
                    running = False
                    events.put({"kind": "status", "status": "paused"})
            if not running or computed_time_fs >= playback_limit_fs:
                continue
            began = time.perf_counter()
            frame = simulation.advance()
            computed_time_fs = frame["time_fs"]
            duration = time.perf_counter() - began
            events.put({"kind": "frame", "frame": frame,
                        "stats": {"steps_per_second": FRAME_STEPS / max(duration, 1e-9)}})
            if frame["progress"] >= 1 or frame["step"] >= MAX_STEPS:
                running = False
                events.put({"kind": "status", "status": "finished",
                            "reason": "finish_line" if frame["progress"] >= 1 else "time_limit"})
            # Run freely until the bounded playback window is full. Waiting
            # between batches never changes integration or random-number use.
    except Exception as error:
        events.put({"kind": "error", "error": str(error)})
    finally:
        if simulation is not None:
            simulation.close()
