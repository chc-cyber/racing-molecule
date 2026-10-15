"""Actual default-surface candidate lifecycle on an isolated 8787 instance."""
import json
import hashlib
from pathlib import Path
import time
import urllib.request
import uuid

BASE = "http://127.0.0.1:8787"
ROOT = Path(__file__).resolve().parents[1]


def request(path, data=None):
    body = None if data is None else json.dumps(data).encode()
    req = urllib.request.Request(BASE + path, body, {"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=20) as response:
        return json.load(response)


def until(get, statuses, seconds=180):
    deadline = time.monotonic() + seconds
    seen = []
    while time.monotonic() < deadline:
        result = get()
        if not seen or seen[-1] != result["status"]:
            seen.append(result["status"])
        if result["status"] in statuses:
            return result, seen
        time.sleep(.2)
    raise AssertionError("Preparation timeout")


def main():
    out = ROOT / "outputs/editor-handoff-default"
    out.mkdir(exist_ok=True)
    files = ["scripts/verify_editor_handoff.py", "molecular_car/server.py",
             "molecular_car/editor/race_handoff.py", "molecular_car/engine.py", "molecular_car/fields.py"]
    report = {"checks": [], "failed": 0, "sourceHashes": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in files}}
    try:
        active, stages = until(lambda: request("/api/state"), {"paused", "error"})
        report["startup"] = {"status": active["status"], "error": active["error"], "stages": stages}
        assert active["status"] == "paused", active["error"]
        assert active["settings"]["surface"] == "fourier-au"
        request("/api/action", {"action": "start"})
        time.sleep(.4)
        before = request("/api/state")
        old_replay = request("/api/replay")
        job_id = str(uuid.uuid4())
        candidate = {"id": job_id, "filename": "tetraphenylbenzene.mol",
                     "content": (ROOT / "examples/tetraphenylbenzene.mol").read_text(encoding="utf-8")}
        request("/api/editor/race/prepare", candidate)
        result, stages = until(lambda: request("/api/editor/race/status?id=" + job_id), {"ready", "error"})
        after = request("/api/state")
        replay = request("/api/replay")
        assert after["generation"] == before["generation"]
        assert after["molecule"] == before["molecule"]
        assert replay["frames"][:len(old_replay["frames"])] == old_replay["frames"]
        assert replay["frames"][-1]["seq"] >= old_replay["frames"][-1]["seq"]
        report["checks"].append({"name": "Default-surface preparation does not touch active race",
            "candidate": result, "stages": stages, "seqBefore": before["frame"]["seq"], "seqAfter": after["frame"]["seq"]})
        if result["status"] == "ready":
            request("/api/editor/race/cancel", {"id": job_id})
        else:
            assert "preparing" in stages, "Expected a full surface-preparation failure after chemistry checking"
            report["checks"].append({"name": "Actual joint-relaxation failure retains old race and frames",
                                      "error": result["error"]})
        # Prepare and promote a known small candidate, then actually advance it.
        good = {"id": str(uuid.uuid4()), "filename": "benzene.mol",
                "content": (ROOT / "examples/benzene.mol").read_text(encoding="utf-8")}
        request("/api/editor/race/prepare", good)
        ready, stages = until(lambda: request("/api/editor/race/status?id=" + good["id"]), {"ready", "error"})
        assert ready["status"] == "ready", ready["error"]
        committed = request("/api/editor/race/commit", {"id": good["id"]})
        same = request("/api/editor/race/commit", {"id": good["id"]})
        assert same == committed
        current = request("/api/state")
        assert current["status"] == "paused" and current["frame"]["step"] == 0
        assert current["generation"] == before["generation"] + 1
        assert current["runtime"]["preparation"]["accepted"]
        assert current["runtime"]["preparation"]["final_coordinate_force_verified"]
        previous = request("/api/replay/previous")
        assert previous["molecule"] == before["molecule"]
        request("/api/action", {"action": "start"})
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            current = request("/api/state")
            if current["frame"]["step"] > 0:
                break
            time.sleep(.1)
        assert current["frame"]["step"] > 0
        request("/api/action", {"action": "pause"})
        report["checks"].append({"name": "Default-surface candidate commits once, retains previous replay and runs in LAMMPS",
                                  "stages": stages, "generation": current["generation"], "step": current["frame"]["step"],
                                  "preparation": current["runtime"]["preparation"]})
        report["passed"] = len(report["checks"])
    except Exception as error:
        report.update(failed=1, error=str(error))
        raise
    finally:
        report["sourcesChangedDuringRun"] = [name for name in files if
            report["sourceHashes"][name] != hashlib.sha256((ROOT / name).read_bytes()).hexdigest()]
        (out / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"passed": report["passed"], "candidateStatus": result["status"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
