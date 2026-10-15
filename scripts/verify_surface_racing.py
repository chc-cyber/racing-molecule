"""Reproducible real-MD races on continuous, frozen Au and mobile Au surfaces.

Run --stage pilot first. Each trial runs in its own interpreter. --stage races
resumes completed per-trial JSON/CSV automatically; --stage nve checks the
new mixed native SAIP/EAM plus Python UFF/field energy accounting. No game API
is called and only outputs/surface-racing-check is written.
"""

from __future__ import annotations

import os
for _variable in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_variable]="1"

import argparse
import csv
from datetime import datetime,timezone
import json
import hashlib
from pathlib import Path
import subprocess
import sys
import time
import traceback
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
OUTPUT=ROOT/"outputs"/"surface-racing-check"
SEEDS=(98245165,2718281,1234567)
MODELS=("continuous","frozen_au","mobile_au")


def source_fingerprint():
    paths=("scripts/verify_surface_racing.py","molecular_car/surface_racing.py",
        "molecular_car/fields.py","molecular_car/molecule.py","molecular_car/au_surface.py",
        "scripts/lammps_runtime.py","examples/benzene.mol","outputs/au-surface-check/au111-relaxed.xyz")
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    temporary.replace(path)


def trial_name(config):
    return f"{config['stage']}-{config['model']}-s{config['strength']:.2f}-seed{config['seed']}-steps{config['steps']}"


def columns_and_row(sample, initial):
    names=["step","time_ps","com_x_A","com_y_A","com_z_A","dx_A","dy_A","dz_A",
        "vx_A_ps","vy_A_ps","vz_A_ps","energy_eV","pe_eV","ke_eV","car_temperature_K","Au_temperature_K",
        "carbon_plane_rms_A","carbon_tilt_deg","bond_min_A","bond_max_A","min_z_A","max_z_A","max_car_force_eV_A",
        "top_Au_mean_z_A","COM_height_above_top_Au_A"]
    displacement=[a-b for a,b in zip(sample["center_angstrom"],initial["center_angstrom"])]
    values=[sample["step"],sample["time_ps"],*sample["center_angstrom"],*displacement,
        *sample["velocity_angstrom_per_ps"],sample["total_energy_eV"],sample["potential_energy_eV"],sample["kinetic_energy_eV"],
        sample["car_temperature_K"],sample["Au_temperature_K"],sample["carbon_plane_rms_angstrom"],
        sample["carbon_plane_tilt_deg"],sample["bond_min_angstrom"],sample["bond_max_angstrom"],
        sample["min_atom_z_angstrom"],sample["max_atom_z_angstrom"],sample["max_car_force_eV_per_angstrom"],
        sample["top_Au_mean_z_angstrom"],sample["COM_height_above_top_Au_angstrom"]]
    for atom,position in enumerate(sample["positions"],1):
        names.extend(f"atom{atom:02d}_{axis}_A" for axis in "xyz")
        values.extend(position)
    return names,values


def actual_blas_threads():
    """Read bundled OpenBLAS's runtime, not the NumPy build metadata."""
    import ctypes
    import numpy as np
    candidates=list((Path(np.__file__).parent.parent/"numpy.libs").glob("*openblas*.dll"))
    for path in candidates:
        try:
            library=ctypes.CDLL(str(path))
            function=getattr(library,"openblas_get_num_threads64_")
            function.restype=ctypes.c_int
            return {"threads":function(),"library":str(path)}
        except (OSError,AttributeError):
            pass
    return {"threads":None,"note":"Runtime API unavailable; requested environment limit is 1."}


def worker(config):
    import numpy as np
    from molecular_car.surface_racing import SurfaceRace,EV_TO_KCAL_MOL,DT_PS
    name=trial_name(config)
    directory=OUTPUT/"trials"/name
    directory.mkdir(parents=True,exist_ok=True)
    summary_path=directory/"summary.json"
    summary={"status":"running","configuration":config,"trial_id":name,"source_sha256":source_fingerprint()}
    save(summary_path,summary)
    simulation=None
    try:
        simulation=SurfaceRace(config["model"],config["strength"],config["temperature"],config["seed"],
            thermostat=config["stage"]!="nve",equilibration_steps=config.get("equilibration_steps",4000),
            timestep_ps=config.get("timestep_ps",DT_PS))
        if config["stage"]=="nve":
            simulation.drive_enabled=False
            simulation.lmp.command("run 0")
            simulation._check()
        description=simulation.description()
        initial=simulation.sample()
        # An arbitrary molecular perturbation checks the converted external
        # UFF+surface/driver force, independently of the native pair models.
        coordinates=simulation.coordinates()
        direction=np.random.default_rng(87123).standard_normal(coordinates.shape)
        direction/=np.linalg.norm(direction)
        energy,force=simulation.external(coordinates)
        epsilon=1e-5
        numeric=(simulation.external(coordinates+epsilon*direction)[0]-simulation.external(coordinates-epsilon*direction)[0])/(2*epsilon)
        analytic=-float(np.sum(force*direction))
        force_error=abs(numeric-analytic)
        if force_error>1e-6: raise AssertionError(f"Converted external force gradient failed: {force_error} eV/A")
        native_pair=float(simulation.lmp.get_thermo("epair"))
        accounting_error=abs(initial["total_energy_eV"]-(energy+native_pair+initial["kinetic_energy_eV"]))
        if accounting_error>1e-7: raise AssertionError("Native plus external energy accounting failed.")
        samples=[initial]
        wall,cpu,sampling_wall=0.0,0.0,0.0
        start=time.perf_counter()
        sample_steps=config.get("sample_steps",200)
        remaining=config["steps"]
        while remaining:
            count=min(sample_steps,remaining)
            timing=simulation.advance(count)
            wall+=timing["wall_seconds"];cpu+=timing["cpu_seconds"]
            sampling_start=time.perf_counter()
            samples.append(simulation.sample())
            sampling_wall+=time.perf_counter()-sampling_start
            remaining-=count
        end_to_end=time.perf_counter()-start
        csv_start=time.perf_counter()
        trajectory=directory/"trajectory.csv"
        with trajectory.open("w",newline="",encoding="utf-8") as handle:
            writer=csv.writer(handle)
            header,_=columns_and_row(initial,initial)
            writer.writerow(header)
            for row in samples: writer.writerow(columns_and_row(row,initial)[1])
        csv_wall=time.perf_counter()-csv_start
        final=samples[-1]
        energies=np.array([row["total_energy_eV"] for row in samples])
        drift=float(np.max(np.abs(energies-energies[0])))
        if config["stage"]=="nve" and drift>0.01:
            raise AssertionError(f"Zero-field NVE absolute energy drift exceeded 0.01 eV: {drift}")
        displacement=np.asarray(final["center_angstrom"])-initial["center_angstrom"]
        duration=config["steps"]*simulation.timestep_ps
        velocities=np.array([row["velocity_angstrom_per_ps"] for row in samples])
        summary.update(status="passed",**description,steps=config["steps"],duration_ps=duration,
            sample_steps=sample_steps,samples=len(samples),trajectory_csv=str(trajectory),csv_columns=header,
            initial=initial,final=final,displacement_angstrom=displacement.tolist(),
            average_forward_velocity_angstrom_per_ps=float(displacement[0]/duration),
            sampled_mean_in_plane_speed_angstrom_per_ps=float(np.linalg.norm(velocities[:,:2],axis=1).mean()),
            maximum_abs_lateral_excursion_angstrom=max(abs(row["center_angstrom"][1]-initial["center_angstrom"][1]) for row in samples),
            maximum_carbon_plane_rms_angstrom=max(row["carbon_plane_rms_angstrom"] for row in samples),
            maximum_carbon_tilt_deg=max(row["carbon_plane_tilt_deg"] for row in samples),
            production_monitor=simulation.monitor,
            performance={"integration_wall_seconds":wall,"integration_cpu_seconds":cpu,
                "integration_seconds_per_step":wall/config["steps"],"steps_per_wall_second":config["steps"]/wall,
                "ps_per_wall_second":duration/wall,"ns_per_day":duration/wall*86.4,
                "sampling_wall_seconds":sampling_wall,"CSV_write_wall_seconds":csv_wall,
                "propagation_plus_sampling_wall_seconds":end_to_end},
            energy_validation={"external_force_gradient_absolute_error_eV_per_angstrom":force_error,
                "native_plus_external_energy_accounting_absolute_error_eV":accounting_error,
                "maximum_absolute_energy_change_eV":drift,
                "maximum_absolute_energy_change_eV_per_car_atom":drift/simulation.n_car,
                "relative_to_initial_absolute_total_energy":drift/max(abs(energies[0]),1e-12),
                "interpretation":"NVE drift" if config["stage"]=="nve" else "Thermostatted/driven energy change is not conservation error"},
            runtime_BLAS=actual_blas_threads())
        if simulation.n_gold:
            np.save(directory/"Au-tail.npy",simulation._all_property("x")[:simulation.n_gold])
        print(f"PASS {name}: dx={displacement[0]:.4f} A; dy={displacement[1]:.4f} A; {config['steps']/wall:.1f} steps/s; setup={simulation.preparation_wall_seconds:.2f}s",flush=True)
    except BaseException as error:
        summary.update(status="failed",error=str(error),traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        if simulation is not None:simulation.close()
        save(summary_path,summary)
    return 0 if summary["status"]=="passed" else 1


def protocol():
    return {"temperature_K":600,"timestep_ps":.0005,"Langevin_damping_ps":.35,"strengths":[.75,.25],
        "seeds":list(SEEDS),"race_target_duration_ps":20,"model_labels":{
            "continuous":"original effective Morse + triangular surface without soft track boundaries",
            "frozen_au":"frozen 2880-atom Au(111) + SAIP; no static Au-Au EAM calculation per step",
            "mobile_au":"bottom 2 layers frozen; 1920 mobile Au with EAM + SAIP and separate Au Langevin bath"},
        "units":{"length":"angstrom","time":"ps","energy":"eV","velocity":"angstrom/ps"},
        "energy_kcal_mol_per_eV":23.0605478306,
        "driver":"same fixed broad gradient: phi=-(0.5*S+0.0015*S*X)*z; X=unwrapped atom x-initial molecular COM x; U=q*phi-0.5*500*alpha*|E|^2",
        "initialization":"same angle-0 planar UFF benzene; fixed x/y, zero-field z relaxation on each own surface; identical NumPy PCG64 car velocity array for paired seed, zero COM, target kinetic temperature",
        "Au_equilibration":"mobile Au only: 2 ps with car held, target 600 K, independent bath seed; not included in production clock",
        "noise":"car Langevin zero no, no temp/com bias; mobile Au uses independent bath with zero yes; matched inputs do not claim bitwise-identical dynamical noise",
        "unwrapping":"native image masks/flags for all car coordinates; continuous original 2.88 A period evaluated on unwrapped positions; explicit Au period 2.884995667 A; no field reset at box crossing",
        "measurement":"MD run calls include callback/neighbour/integration costs; setup, sampling and CSV I/O separately timed inside each child process; 1 LAMMPS/BLAS thread",
        "limitations":["Benzene UFF has fixed topology; no reactions, electron/charge transfer or STM excitation.",
            "The shared driving field is a game effective response, not an experimentally calibrated electric field.",
            "SAIP was fitted to graphene/benzene on Au; molecular tilt/bending are monitored because H normals default to z.",
            "Small periodic unreconstructed Au slab; repeated images exclude steps/defects/finite-track physics.",
            "Only three random seeds and a finite window; no converged diffusion coefficient or race-failure probability.",
            "Mobile Au includes extra surface motion and dissipation and is an extended comparison, not solely a corrugation substitution."]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",choices=("pilot","nve","races","all"),default="pilot")
    parser.add_argument("--window-ps",type=float,default=20)
    parser.add_argument("--mobile-ps",type=float,default=20)
    parser.add_argument("--worker",type=Path)
    args=parser.parse_args()
    if args.worker:return worker(json.loads(args.worker.read_text(encoding="utf-8")))
    OUTPUT.mkdir(parents=True,exist_ok=True)
    report_path=OUTPUT/"report.json"
    if report_path.exists():report=json.loads(report_path.read_text(encoding="utf-8"))
    else:report={"protocol":protocol(),"trials":[],"status":"running"}
    report["updated_UTC"]=datetime.now(timezone.utc).isoformat()
    report["protocol"]=protocol()
    report["source_sha256"]=source_fingerprint()
    if args.stage in ("races","all"):
        report["protocol"]["requested_window_ps_by_model"]={model:args.mobile_ps if model=="mobile_au" else args.window_ps for model in MODELS}
    report["status"]="running"
    save(report_path,report)
    stages=("pilot","nve","races") if args.stage=="all" else (args.stage,)
    configurations=[]
    for stage in stages:
        seeds=SEEDS if stage=="races" else SEEDS[:1]
        strengths=(.75,.25) if stage=="races" else ((0,) if stage=="nve" else (.75,))
        for seed_index,seed in enumerate(seeds):
            for strength_index,strength in enumerate(strengths):
                offset=(seed_index+strength_index)%len(MODELS)
                order=MODELS[offset:]+MODELS[:offset] if stage=="races" else MODELS
                for model in order:
                    duration=args.mobile_ps if model=="mobile_au" else args.window_ps
                    steps=1000 if stage=="pilot" else 2000 if stage=="nve" else int(round(duration/.0005))
                    configurations.append({"stage":stage,"model":model,"strength":strength,"temperature":10 if stage=="nve" else 600,
                        "seed":seed,"steps":steps,"sample_steps":100 if stage=="pilot" else 200,
                        "equilibration_steps":0 if stage=="nve" else 4000})
    for index,config in enumerate(configurations,1):
        name=trial_name(config)
        summary_path=OUTPUT/"trials"/name/"summary.json"
        previous=json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else None
        if previous and previous["status"]=="passed" and previous.get("configuration")==config and previous.get("source_sha256")==source_fingerprint():
            result=previous
            print(f"RESUME {index}/{len(configurations)} {name}",flush=True)
        else:
            config_path=OUTPUT/"trials"/name/"configuration.json"
            save(config_path,config)
            print(f"START {index}/{len(configurations)} {name}",flush=True)
            process=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--worker",str(config_path)],cwd=ROOT,check=False)
            result=json.loads(summary_path.read_text(encoding="utf-8"))
            if process.returncode:report["status"]="failed"
        report["trials"]=[trial for trial in report["trials"] if trial["trial_id"]!=name]+[result]
        save(report_path,report)
        if result["status"]!="passed":return 1
    report["status"]="passed"
    report["last_completed_stage"]=args.stage
    report["updated_UTC"]=datetime.now(timezone.utc).isoformat()
    save(report_path,report)
    print(f"Report: {report_path}",flush=True)
    return 0


if __name__=="__main__":raise SystemExit(main())
