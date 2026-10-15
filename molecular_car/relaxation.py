"""Deterministic joint geometry relaxation before race velocities and bath.

The objective and forces come from the existing LAMMPS external callback:
UFF + the selected surface + track walls, with the drive disabled. A native
minimizer's return is never sufficient: the full 3N Cartesian force two-norm
must independently meet the stated tolerance before a race may start.
"""
from __future__ import annotations

import time
import ctypes
import numpy as np
from .preparation_bfgs import BFGS_OPTIONS, PreparationBudgetStop, dense_bfgs

PROTOCOL_VERSION = "joint-force-relaxation-v2"
GLOBAL_FORCE_TOLERANCE = 1e-4  # kcal/mol/angstrom, two-norm of all 3N forces
WALL_TIME_LIMIT_SECONDS = 30.
FORCE_EVALUATION_LIMIT = 25_000
# Generous numerical admission bounds, not chemical validity criteria. They
# reject pathological parameter/geometry output rather than spend the entire
# minimization budget on billion-kcal energies and million-kcal/A forces.
INITIAL_MAX_ATOM_FORCE = 1e5  # kcal/mol/angstrom
INITIAL_MAX_ABSOLUTE_ENERGY = 1e8  # kcal/mol
# A short damped descent handles initial adsorption, then Newton's method
# resolves the disparate stiff UFF and soft adsorption/translation scales.
# CG is an alternative if Newton's line search stalls; a final Newton pass
# can then finish the remaining low-curvature directions. No random dynamics.
STAGES = (("fire", 350, 3500), ("hftn", 500, 10000),
          ("cg", 4000, 12000), ("hftn", 500, 5000))
STOP_REASON_LABELS = {
    "wall-time-budget": "准备等待时间已达到上限",
    "force-evaluation-budget": "准备阶段的求力计算额度已用完",
    "pathological-initial-energy-or-force": "当前构象或化学环境产生异常大的能量或受力",
    "all-minimizer-stages-exhausted-without-force-convergence": "已尝试全部几何优化阶段，残余力仍未达到要求",
}


class JointRelaxationError(RuntimeError):
    """No race is permitted when the joint conservative force fails its gate."""
    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        super().__init__(
            "赛车与表面的联合几何松弛未通过准备检查，不能开始比赛。"
            f"全局残余力 {diagnostics['global_force_norm_kcal_mol_A']:.6g} "
            f"kcal/mol/Å，要求 ≤ {diagnostics['global_force_tolerance_kcal_mol_A']:.6g}；"
            f"停止原因：{STOP_REASON_LABELS.get(diagnostics['termination_reason'], '残余力仍未达到要求')}。"
            "请整理初始构象或改用其他表面。")


def _force_diagnostics(lmp, atom_count):
    raw = lmp.gather_atoms("f", 1, 3)
    forces = np.ctypeslib.as_array(raw, shape=(atom_count * 3,)).reshape(atom_count, 3)
    energy = float(lmp.get_thermo("pe"))
    if not np.isfinite(forces).all() or not np.isfinite(energy):
        raise RuntimeError("联合几何松弛出现非有限能量或受力。")
    norm = float(np.linalg.norm(forces))
    return {"global_force_norm_kcal_mol_A": norm,
            "maximum_atom_force_kcal_mol_A": float(np.linalg.norm(forces, axis=1).max()),
            "rms_atom_force_kcal_mol_A": norm / np.sqrt(atom_count),
            "maximum_force_component_kcal_mol_A": float(np.abs(forces).max()),
            "potential_energy_kcal_mol": energy}


def relax_joint_geometry(lmp, atom_count, force_callback, check_callback, *,
                         force_tolerance=GLOBAL_FORCE_TOLERANCE,
                         wall_time_limit=WALL_TIME_LIMIT_SECONDS,
                         evaluation_limit=FORCE_EVALUATION_LIMIT, stages=STAGES,
                         bfgs_rescue=True):
    """Minimize the unchanged conservative game objective, or fail explicitly.

    Native ``timer timeout`` stops at the next outer minimizer iteration.
    The additional callback budget limits total true force evaluations across
    stages and uses LAMMPS's safe soft-stop API. Both controls are removed
    before thermal velocities, NVE integration or Langevin are installed.
    Neither the race timestep nor the force/model parameters are changed.
    """
    if atom_count <= 0 or not np.isfinite(force_tolerance) or force_tolerance <= 0:
        raise ValueError("A positive atom count and finite force tolerance are required.")
    if not np.isfinite(wall_time_limit) or wall_time_limit <= 0 or evaluation_limit < 1:
        raise ValueError("Relaxation budgets must be positive and finite.")
    started = time.perf_counter()
    evaluations = 0
    optimization_evaluations = 0
    verification_evaluations = 0
    verifying = False
    budget_stop = None
    records = []
    report = {"protocol_version": PROTOCOL_VERSION,
              "method": "FIRE / Hessian-free truncated Newton / optional CG / budgeted curvature-initialized BFGS rescue; joint UFF/surface/walls, zero drive, before thermal velocities and Langevin",
              "energy_tolerance": 0.,
              "global_force_tolerance_kcal_mol_A": float(force_tolerance),
              "force_criterion": "two-norm of the full 3N Cartesian conservative force vector; independently checked after every optimization stage",
              "iteration_limit": sum(stage[1] for stage in stages) + (BFGS_OPTIONS["maximum_iterations"] if bfgs_rescue else 0),
              "evaluation_limit": int(evaluation_limit),
              "evaluation_limit_qualification": "Limits initial, native minimizer, BFGS Hessian-probe and line-search force calls together. One additional final-coordinate verification per completed stage is counted separately; a reached budget never admits a race.",
              "wall_time_limit_seconds": float(wall_time_limit),
              "initial_admission_bounds": {"maximum_atom_force_kcal_mol_A": INITIAL_MAX_ATOM_FORCE,
                                           "maximum_absolute_energy_kcal_mol": INITIAL_MAX_ABSOLUTE_ENERGY,
                                           "qualification": "Generous numerical geometry/parameter guards, not a calibration or chemical validity claim."},
              "wall_time_limit_qualification": "Native soft stop at the next outer iteration; BFGS checks before each force call and after final verification. A running force evaluation or matrix operation may finish first.",
              "bfgs_rescue_options": dict(BFGS_OPTIONS) if bfgs_rescue else None,
              "stages": records}

    def counted_force(*args):
        nonlocal evaluations, optimization_evaluations, verification_evaluations, budget_stop
        evaluations += 1
        if verifying:
            verification_evaluations += 1
        else:
            optimization_evaluations += 1
        force_callback(*args)
        if not verifying and optimization_evaluations >= evaluation_limit and budget_stop is None:
            budget_stop = "force-evaluation-budget"
            lmp.force_timeout()

    try:
        lmp.set_fix_external_callback("molecular", counted_force)
        lmp.command("timer timeout off every 10")
        # A zero-step refresh is an energy/force evaluation, not MD warming.
        # Evaluate before enabling a native timer: an expired timer must never
        # skip this call and leave uninitialized/old forces at the gate.
        lmp.command("run 0")
        check_callback()
        state = _force_diagnostics(lmp, atom_count)
        report["initial"] = dict(state)
        if (state["maximum_atom_force_kcal_mol_A"] > INITIAL_MAX_ATOM_FORCE
                or abs(state["potential_energy_kcal_mol"]) > INITIAL_MAX_ABSOLUTE_ENERGY):
            report.update(state)
            report.update(iterations=0, force_callback_evaluations=evaluations,
                          optimization_force_callback_evaluations=optimization_evaluations,
                          final_verification_force_callback_evaluations=verification_evaluations,
                          wall_seconds=time.perf_counter() - started, force_tolerance_met=False,
                          accepted=False, final_coordinate_force_verified=True,
                          termination_reason="pathological-initial-energy-or-force")
            raise JointRelaxationError(report)
        for style, max_iterations, max_evaluations in stages:
            if budget_stop is not None or time.perf_counter() - started >= wall_time_limit:
                budget_stop = budget_stop or "wall-time-budget"
                break
            if state["global_force_norm_kcal_mol_A"] <= force_tolerance:
                break
            before_step = int(lmp.get_thermo("step"))
            before_evaluations = evaluations
            stage_started = time.perf_counter()
            before = dict(state)
            lmp.command(f"min_style {style}")
            # HFTN ignores min_modify; do not imply that its trust-region
            # Newton steps obey FIRE/CG's 0.1 Å displacement option.
            if style != "hftn":
                lmp.command("min_modify dmax 0.1")
            remaining = max(1, evaluation_limit - optimization_evaluations)
            lmp.command(f"timer timeout {max(1e-12, wall_time_limit - (time.perf_counter() - started)):.12g} every 1")
            lmp.command(f"minimize 0.0 {force_tolerance:.12g} {max_iterations} {min(max_evaluations, remaining)}")
            check_callback()
            # Newton Hessian probes perturb positions internally. Refresh
            # at final coordinates even after a soft stop, with the timer off.
            # This one verification is counted and cannot re-enable acceptance
            # after a budget stop. It does not advance coordinates or MD time.
            lmp.command("timer timeout off every 10")
            verifying = True
            try:
                lmp.command("run 0")
                check_callback()
            finally:
                verifying = False
            state = _force_diagnostics(lmp, atom_count)
            iterations = int(lmp.get_thermo("step")) - before_step
            elapsed = time.perf_counter() - started
            if budget_stop is None and elapsed >= wall_time_limit:
                budget_stop = "wall-time-budget"
            force_met = state["global_force_norm_kcal_mol_A"] <= force_tolerance
            passed = force_met and budget_stop is None
            reason = ("independently-verified-force-tolerance" if passed else budget_stop
                      or ("iteration-limit" if iterations >= max_iterations else "native-stopped-without-force-convergence"))
            records.append({"minimizer": style, "maximum_iterations": max_iterations,
                            "maximum_native_force_evaluations": min(max_evaluations, remaining),
                            "iterations": iterations, "force_callback_evaluations": evaluations - before_evaluations,
                            "wall_seconds": time.perf_counter() - stage_started,
                            "initial": before, "final": dict(state),
                            "force_tolerance_met": bool(force_met), "accepted": bool(passed),
                            "final_coordinate_force_verified": True, "termination_reason": reason,
                            "native_stop_qualification": "Actual force gate is independently checked. A native early stop without force convergence is not accepted; detailed native stop text is available when LAMMPS logging is enabled."})
        if (bfgs_rescue and budget_stop is None
                and state["global_force_norm_kcal_mol_A"] > force_tolerance):
            if time.perf_counter() - started >= wall_time_limit:
                budget_stop = "wall-time-budget"
            elif evaluation_limit - optimization_evaluations > 6 * atom_count + 1:
                before, before_evaluations = dict(state), evaluations
                stage_started = time.perf_counter()
                raw = lmp.gather_atoms("x", 1, 3)
                checkpoint = np.ctypeslib.as_array(raw, shape=(atom_count * 3,)).reshape(atom_count, 3).copy()

                def scatter(x):
                    lmp.scatter_atoms("x", 1, 3, (ctypes.c_double * (atom_count * 3))(*x.ravel()))

                def evaluate(x):
                    nonlocal budget_stop
                    if time.perf_counter() - started >= wall_time_limit:
                        budget_stop = budget_stop or "wall-time-budget"
                    if budget_stop is not None:
                        raise PreparationBudgetStop(budget_stop)
                    scatter(x)
                    lmp.command("run 0")
                    check_callback()
                    if budget_stop is not None:
                        raise PreparationBudgetStop(budget_stop)
                    raw_force = lmp.gather_atoms("f", 1, 3)
                    force = np.ctypeslib.as_array(raw_force, shape=(atom_count * 3,)).reshape(atom_count, 3).copy()
                    return float(lmp.get_thermo("pe")), force

                try:
                    checkpoint, rescue = dense_bfgs(checkpoint, evaluate, force_tolerance)
                finally:
                    # Every probe and rejected trial temporarily changes native
                    # coordinates. Restore the accepted checkpoint even on error.
                    lmp.command("timer timeout off every 10")
                    scatter(checkpoint)
                    verifying = True
                    try:
                        lmp.command("run 0")
                        check_callback()
                    finally:
                        verifying = False
                state = _force_diagnostics(lmp, atom_count)
                if time.perf_counter() - started >= wall_time_limit:
                    budget_stop = budget_stop or "wall-time-budget"
                met = state["global_force_norm_kcal_mol_A"] <= force_tolerance
                records.append({"minimizer": "curvature-initialized-bfgs",
                                "maximum_iterations": BFGS_OPTIONS["maximum_iterations"],
                                "iterations": rescue["iterations"],
                                "force_callback_evaluations": evaluations - before_evaluations,
                                "wall_seconds": time.perf_counter() - stage_started,
                                "initial": before, "final": dict(state),
                                "force_tolerance_met": bool(met), "accepted": bool(met and budget_stop is None),
                                "final_coordinate_force_verified": True,
                                "termination_reason": budget_stop or rescue["termination_reason"],
                                "optimizer": rescue})
            else:
                report["bfgs_rescue_skipped"] = "insufficient-shared-evaluation-budget-for-full-curvature-probe"
        if time.perf_counter() - started >= wall_time_limit:
            budget_stop = budget_stop or "wall-time-budget"
        report.update(state)
        report.update(iterations=sum(record["iterations"] for record in records),
                      force_callback_evaluations=evaluations,
                      optimization_force_callback_evaluations=optimization_evaluations,
                      final_verification_force_callback_evaluations=verification_evaluations,
                      wall_seconds=time.perf_counter() - started,
                      final_coordinate_force_verified=True,
                      force_tolerance_met=bool(state["global_force_norm_kcal_mol_A"] <= force_tolerance))
        report["accepted"] = report["force_tolerance_met"] and budget_stop is None
        report["termination_reason"] = ("independently-verified-force-tolerance" if report["accepted"]
                                         else budget_stop or "all-minimizer-stages-exhausted-without-force-convergence")
        if not report["accepted"]:
            raise JointRelaxationError(report)
        return report
    finally:
        # Expired native timers otherwise silently skip later run commands.
        try:
            lmp.command("timer timeout off every 10")
        finally:
            lmp.set_fix_external_callback("molecular", force_callback)
