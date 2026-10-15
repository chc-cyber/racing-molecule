"""Exercise the running local game's drive API and restore its default scene.

Run after restarting the current server on http://127.0.0.1:8765.  This uses
only the public HTTP API; it does not start a server or import MD internals.
The final scene is the paused 52-atom default car at step zero.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
import time
import traceback
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import ProxyHandler, Request, build_opener


ROOT = Path(__file__).resolve().parents[1]
REPORT_PATH = ROOT / "outputs" / "drive-api-check" / "api-report.json"
DEFAULT_DRIVE = {"mode": "gradient", "strength": 0.75}
DEFAULT_SETTINGS = {"temperature": 10.0, "damping_fs": 350.0, "corrugation": 0.08}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


class LocalAPI:
    def __init__(self, base: str):
        parsed = urlsplit(base)
        require(parsed.scheme == "http" and parsed.hostname in {"127.0.0.1", "localhost"}
                and parsed.path in {"", "/"} and not parsed.query and not parsed.fragment
                and parsed.username is None, "Use a local HTTP game URL")
        self.base = base.rstrip("/")
        # Do not send local game requests through an inherited system proxy.
        self.opener = build_opener(ProxyHandler({}))

    def request(self, path: str, payload=None, *, expected=200, json_response=True):
        data = None if payload is None else json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
        headers = {"Origin": self.base}
        if data is not None:
            headers["Content-Type"] = "application/json"
        request = Request(self.base + path, data=data, headers=headers,
                          method="GET" if data is None else "POST")
        try:
            response = self.opener.open(request, timeout=10)
        except HTTPError as error:
            response = error
        with response:
            status = response.status
            raw = response.read().decode("utf-8")
        require(status == expected, f"{path}: expected HTTP {expected}, got {status}: {raw[:300]}")
        return json.loads(raw) if json_response else raw

    def state(self):
        return self.request("/api/state")

    def wait(self, predicate, *, timeout=30):
        deadline = time.monotonic() + timeout
        latest = None
        while time.monotonic() < deadline:
            latest = self.state()
            if latest["status"] == "error":
                raise AssertionError(f"Simulation error: {latest.get('error')}")
            if predicate(latest):
                return latest
            time.sleep(0.04)
        raise AssertionError(f"Timed out waiting for game state: status={latest.get('status') if latest else None}")

    def paused(self):
        current = self.state()
        if current["status"] == "running":
            self.request("/api/action", {"action": "pause"})
        return self.wait(lambda state: state["status"] == "paused" and state["frame"] is not None)

    def change(self, payload, predicate):
        before = self.paused()
        self.request("/api/control", payload)
        after = self.wait(lambda state: state["status"] == "paused"
                         and state["generation"] == before["generation"]
                         and state["frame"]["seq"] > before["frame"]["seq"]
                         and predicate(state))
        require(after["frame"]["step"] == before["frame"]["step"], "Paused control change advanced time")
        require(after["frame"]["positions"] == before["frame"]["positions"], "Paused control change altered coordinates")
        return before, after


def compact(state):
    return {
        "generation": state["generation"],
        "status": state["status"],
        "atoms": len(state["molecule"]["atoms"]) if state.get("molecule") else None,
        "step": state["frame"]["step"] if state.get("frame") else None,
        "control": state["control"],
        "settings": state["settings"],
    }


def synchronized(state, strength):
    return (state["control"]["drive"]["strength"] == strength
            and state["control"]["tip"]["strength"] == strength
            and state["frame"]["drive"]["strength"] == strength
            and state["frame"]["tip"]["strength"] == strength)


def restore_default(api: LocalAPI):
    state = api.state()
    if state["status"] == "running":
        api.request("/api/action", {"action": "pause"})
    elif state["status"] in {"error", "finished", "idle"}:
        api.request("/api/load", {"example": "tetraphenylbenzene.mol"})
    api.wait(lambda state: state["status"] == "paused" and state["frame"] is not None)
    api.request("/api/control", {"drive": dict(DEFAULT_DRIVE), **DEFAULT_SETTINGS, "autopilot": True})
    # Load uses independently retained accepted controls, even before worker ACK.
    response = api.request("/api/load", {"example": "tetraphenylbenzene.mol"})
    generation = response["generation"]
    final = api.wait(lambda state: state["generation"] == generation and state["status"] == "paused" and state["frame"] is not None)
    require(final["frame"]["step"] == 0 and len(final["molecule"]["atoms"]) == 52, "Default scene was not restored")
    require(final["control"]["drive"] == DEFAULT_DRIVE and synchronized(final, 0.75), "Default drive was not restored")
    require(final["settings"] == DEFAULT_SETTINGS, "Default bath/surface settings were not restored")
    return compact(final)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default="http://127.0.0.1:8765")
    args = parser.parse_args()
    api = LocalAPI(args.base_url)
    started = time.perf_counter()
    report = {"passed": False, "created_utc": datetime.now(timezone.utc).isoformat(),
              "base_url": api.base, "checks": []}

    def record(name, **details):
        report["checks"].append({"name": name, "passed": True, **details})
        print(f"Passed: {name}", flush=True)

    try:
        initial = api.paused()
        require(initial["track"]["finish_x"] == 1000.0 and initial["track"]["half_width"] == 250.0,
                "Track is not 100 nm by 50 nm")
        require(initial["control"]["drive"] == DEFAULT_DRIVE and synchronized(initial, 0.75), "Restarted default drive is not gradient 0.75")
        require(len(initial["molecule"]["atoms"]) == 52, "Default car is not the 52-atom fixture")
        record("default gradient and 1000 by 500 angstrom track", state=compact(initial))

        before, after = api.change({"drive": {"mode": "tip"}}, lambda state: state["frame"]["drive"]["mode"] == "tip")
        record("paused mode switch preserves step and coordinates", step=after["frame"]["step"], generation=after["generation"])
        _, after = api.change({"drive": {"strength": 0.4}}, lambda state: synchronized(state, 0.4))
        record("drive strength synchronizes tip strength", strength=0.4)
        _, after = api.change({"tip": {"strength": 0.63}}, lambda state: synchronized(state, 0.63))
        record("tip strength synchronizes drive strength", strength=0.63)

        generation = after["generation"]
        wanted_drive = {"mode": "tip", "strength": 0.42}
        wanted_settings = {"temperature": 37.0, "damping_fs": 650.0, "corrugation": 0.12}
        api.request("/api/control", {"drive": wanted_drive, **wanted_settings})
        # Deliberately do not poll for acknowledgment between control and reset.
        api.request("/api/action", {"action": "reset"})
        reset = api.wait(lambda state: state["generation"] > generation and state["status"] == "paused" and state["frame"] is not None)
        require(reset["control"]["drive"] == wanted_drive and synchronized(reset, 0.42), "Immediate reset lost accepted drive controls")
        require(reset["settings"] == wanted_settings, "Immediate reset lost accepted temperature/damping/corrugation")
        require(reset["frame"]["step"] == 0, "Reset did not reset the integration step")
        record("immediate reset retains accepted controls before worker acknowledgment", state=compact(reset))

        rejected = []
        invalid = (
            {"drive": {"mode": "invalid"}},
            {"drive": {"strength": -0.01}},
            {"drive": {"strength": 1.01}},
            {"drive": {"strength": True}},
            {"tip": {"x": 1006.1}},
            {"tip": {"y": 250.1}},
            {"tip": {"z": 3.9}},
            {"tip": {"strength": 0.6}, "drive": {"strength": 0.5}},
        )
        for payload in invalid:
            result = api.request("/api/control", payload, expected=400)
            require(bool(result.get("error")), "Rejected controls lack a readable error")
            rejected.append({"request": payload, "error": result["error"], "status": 400})
        unchanged = api.state()
        require(compact(unchanged) == compact(reset) and unchanged["frame"]["positions"] == reset["frame"]["positions"],
                "Invalid controls changed the simulation")
        record("invalid modes strengths and tip bounds rejected without state changes", rejected=rejected)

        # Start the pulse from a fresh default scene; strength 0.75 clamps to 1.0.
        restore_default(api)
        pulse_start = api.paused()
        api.request("/api/action", {"action": "pulse"})
        pulse_end = api.wait(lambda state: state["generation"] == pulse_start["generation"]
                            and state["status"] == "paused"
                            and state["frame"]["step"] >= pulse_start["frame"]["step"]+800)
        require(pulse_end["frame"]["step"] == pulse_start["frame"]["step"]+800, "Single gradient pulse did not stop at 800 steps")
        replay = api.request("/api/replay")
        pulse_frames = [frame for frame in replay["frames"] if frame["step"] > pulse_start["frame"]["step"]]
        require(len(pulse_frames) == 8, "Pulse replay does not contain the eight 100-step frames")
        for frame in pulse_frames:
            require(frame["drive"] == DEFAULT_DRIVE, "Pulse changed nominal driver controls")
            require(all(math.isfinite(number) for number in frame["positions"]) and math.isfinite(frame["energy"]), "Pulse produced nonfinite trajectory data")
            active = frame["step"] < pulse_start["frame"]["step"]+800
            require(frame["pulse_active"] == active, "Pulse-active metadata is incorrect")
            require(frame["effective_drive"] == {"mode": "gradient", "strength": 1.0 if active else 0.75},
                    "Effective pulse strength is not recorded")
            require(all(frame["tip"][axis] == pulse_start["frame"]["tip"][axis] for axis in ("x", "y", "z")),
                    "Gradient pulse moved the tip")
        time.sleep(0.12)
        require(api.state()["frame"]["step"] == pulse_end["frame"]["step"], "Pulse continued integrating after pause")
        record("gradient pulse advances 800 steps then pauses with effective drive metadata",
               frames=[{key: frame[key] for key in ("step", "drive", "effective_drive", "pulse_active", "pulse_until_step")} for frame in pulse_frames])

        runtime = replay["runtime"]
        require(runtime["gradient"] == {"base_field": 0.5, "field_slope": 0.0015, "polarization_game_gain": 500.0}, "Replay lacks exact gradient model parameters")
        require(runtime["langevin_seed"] and runtime["velocity_seed"] and runtime["timestep_fs"] == 0.5, "Replay lacks reproducibility metadata")
        require(replay["control"]["drive"] == DEFAULT_DRIVE and all("drive" in frame and "effective_drive" in frame for frame in replay["frames"]),
                "Replay omits nominal or effective drive data")
        record("replay contains runtime parameters and nominal/effective drive", runtime=runtime, frames=len(replay["frames"]))

        editor_before = api.state()
        html = api.request("/editor/index.html", json_response=False)
        require("<html" in html.lower(), "Editor page is not served as HTML")
        content = (ROOT / "examples" / "benzene.mol").read_text(encoding="utf-8")
        imported = api.request("/api/editor/import", {"filename": "benzene.mol", "content": content})
        require(len(imported["project"]["atoms"]) == 12 and len(imported["project"]["bonds"]) == 12, "Editor benzene import failed")
        require({atom["element"] for atom in imported["project"]["atoms"]} == {"C", "H"}, "Editor benzene identity changed")
        editor_after = api.state()
        require(editor_after["generation"] == editor_before["generation"]
                and editor_after["frame"]["step"] == editor_before["frame"]["step"]
                and editor_after["frame"]["positions"] == editor_before["frame"]["positions"],
                "Editor import disturbed the racing simulation")
        record("editor page and benzene import remain available without changing race generation",
               generation=editor_after["generation"], editor_atoms=12, editor_bonds=12)
        report["passed"] = True
    except Exception as error:
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
        print(f"Drive API check failed: {error}", file=sys.stderr, flush=True)
    finally:
        try:
            report["restored_default_scene"] = restore_default(api)
            print("Restored: paused 52-atom car, gradient 0.75, 10 K, 350 fs, corrugation 0.08.", flush=True)
        except Exception as error:
            report["passed"] = False
            report["cleanup_error"] = str(error)
            print(f"Default scene restoration failed: {error}", file=sys.stderr, flush=True)
        report["wall_seconds"] = time.perf_counter()-started
        REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
        REPORT_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
        print(f"Report: {REPORT_PATH}", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
