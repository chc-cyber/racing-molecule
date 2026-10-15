"""Prepare a candidate in isolation, then promote its live, paused worker.

All state changes use the controller's lifecycle -> state lock order. The
existing /api/load path remains deliberately unchanged for older clients.
"""
import copy
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import queue
import threading
import time
import uuid

from ..start_pose import validate_start_pose
from ..competition import validate_competition, engine_parameters


def candidate_worker(content, filename, settings, commands, events, drive, pose):
    # This hook is local to the spawned process. It reports progress around the
    # existing preparation function; its result is used by Simulation itself.
    # No second UFF/embedding pass and no change to the MD preparation algorithm.
    from .. import engine
    original = engine.prepare_molecule

    def checked(*args, **kwargs):
        events.put({"kind": "stage", "status": "checking"})
        prepared = original(*args, **kwargs)
        events.put({"kind": "stage", "status": "preparing"})
        return prepared

    engine.prepare_molecule = checked
    try:
        engine.worker_main(content, filename, settings, commands, events, drive, pose)
    finally:
        engine.prepare_molecule = original


def stop_worker(process, commands, events):
    """Stop captured worker references, never the controller's new worker."""
    if process is None:
        return
    if process.is_alive():
        try:
            commands.put({"action": "stop"})
        except (OSError, ValueError):
            pass
        process.join(timeout=.3)
    if process.is_alive():
        process.terminate()
        process.join(timeout=2)
    for channel in (commands, events):
        if channel is not None:
            channel.cancel_join_thread()
            channel.close()


class RaceHandoff:
    PREPARATION_SECONDS = 180
    READY_SECONDS = 600

    def __init__(self, game):
        self.game = game
        self.jobs = {}

    def _job(self, job_id):
        job = self.jobs.get(job_id)
        if job is None:
            raise ValueError("这次准备已过期，请重新准备。")
        return job

    @staticmethod
    def _public(job):
        return {key: copy.deepcopy(job[key]) for key in
                ("id", "status", "error", "baseGeneration", "generation", "atomCount")}

    def status(self, job_id):
        with self.game.lock:
            return self._public(self._job(job_id))

    def begin(self, data):
        job_id = data.get("id")
        try:
            if str(uuid.UUID(job_id)) != job_id:
                raise ValueError()
        except (ValueError, TypeError, AttributeError):
            raise ValueError("准备编号无效。") from None
        content, filename = data.get("content"), data.get("filename")
        if not isinstance(content, str) or not content.strip():
            raise ValueError("分子文件不能为空。")
        if not isinstance(filename, str) or Path(filename).suffix.lower() not in (".mol", ".sdf"):
            raise ValueError("请使用 MOL 或单分子 SDF 文件。")
        if len(content.encode("utf-8")) > 2_000_000:
            raise ValueError("文件超过 2 MB。")
        competition = validate_competition(data.get("competition"))
        settings, drive = engine_parameters(competition)
        fingerprint = hashlib.sha256((filename + "\0" + content + "\0" + json.dumps(competition, sort_keys=True)).encode("utf-8")).hexdigest()
        game = self.game
        with game.lifecycle:
            with game.lock:
                if job_id in self.jobs:
                    job = self._job(job_id)
                    if job["fingerprint"] is not None and job["fingerprint"] != fingerprint:
                        raise ValueError("准备编号对应的分子已改变，请重新准备。")
                    return self._public(job)
                if any(j["status"] in ("checking", "preparing", "ready") for j in self.jobs.values()):
                    raise ValueError("另一辆赛车正在准备，请稍后重试。")
                while len(self.jobs) >= 8:
                    del self.jobs[next(iter(self.jobs))]
                job = {"id": job_id, "status": "checking", "error": None,
                       "baseGeneration": game.generation, "generation": None, "atomCount": None,
                       "fingerprint": fingerprint, "input": (content, Path(filename).name),
                       "settings": settings, "drive": drive, "competition": competition,
                       "baseline_settings": dict(game.requested_settings), "baseline_drive": dict(game.requested_drive),
                       "pose": copy.deepcopy(competition["start_pose"]), "process": None, "commands": None,
                       "events": None, "ready": None, "created": time.monotonic()}
                self.jobs[job_id] = job
            context = mp.get_context("spawn")
            job["commands"], job["events"] = context.Queue(), context.Queue(maxsize=32)
            job["process"] = context.Process(target=candidate_worker, args=(content, filename,
                job["settings"], job["commands"], job["events"], job["drive"], job["pose"]), daemon=True)
            try:
                job["process"].start()
            except Exception:
                self._dispose(job, "error", "无法启动准备进程，请重试。")
                return self._public(job)
            threading.Thread(target=self._receive, args=(job,), daemon=True).start()
            timer = threading.Timer(self.PREPARATION_SECONDS, self._expire, args=(job_id, False))
            timer.daemon = True
            timer.start()
            return self.status(job_id)

    def _dispose(self, job, status, error=None):
        with self.game.lock:
            job.update(status=status, error=error, ready=None)
            refs = job["process"], job["commands"], job["events"]
            job.update(process=None, commands=None, events=None)
        stop_worker(*refs)

    def _receive(self, job):
        process, events = job["process"], job["events"]
        while True:
            with self.game.lock:
                if job["status"] not in ("checking", "preparing"):
                    return
            try:
                event = events.get(timeout=.2)
            except queue.Empty:
                if process.is_alive():
                    continue
                event = {"kind": "error", "error": f"准备进程意外退出（{process.exitcode}），请重试。"}
            except (OSError, ValueError, EOFError):
                return
            with self.game.lock:
                if job["status"] not in ("checking", "preparing"):
                    return
                if event["kind"] == "stage":
                    job["status"] = event["status"]
                elif event["kind"] == "ready":
                    job.update(status="ready", ready=event, atomCount=len(event["molecule"]["atoms"]))
                    timer = threading.Timer(self.READY_SECONDS, self._expire, args=(job["id"], True))
                    timer.daemon = True
                    timer.start()
                    return
                elif event["kind"] == "error":
                    job.update(status="error", error=event["error"])
                else:
                    continue
            if event["kind"] == "error":
                with self.game.lifecycle:
                    self._dispose(job, "error", event["error"])
                return

    def _expire(self, job_id, ready):
        with self.game.lifecycle:
            with self.game.lock:
                job = self.jobs.get(job_id)
                if not job or job["status"] not in (("ready",) if ready else ("checking", "preparing")):
                    return
            self._dispose(job, "error", "准备已过期，请重新准备。" if ready else "准备超时，请调整结构后重试。")

    def cancel(self, job_id):
        with self.game.lifecycle:
            if job_id not in self.jobs:
                try:
                    if str(uuid.UUID(job_id)) != job_id:
                        raise ValueError()
                except (ValueError, TypeError, AttributeError):
                    raise ValueError("准备编号无效。") from None
                # A cancel may arrive before a retried/slow prepare request.
                # Keep a tombstone so that same request cannot start afterwards.
                with self.game.lock:
                    while len(self.jobs) >= 8:
                        terminal = next((key for key, value in self.jobs.items()
                                         if value["status"] not in ("checking", "preparing", "ready")), None)
                        if terminal is None:
                            break
                        del self.jobs[terminal]
                    self.jobs[job_id] = {"id": job_id, "status": "cancelled", "error": None,
                        "baseGeneration": self.game.generation, "generation": None, "atomCount": None,
                        "fingerprint": None}
            job = self._job(job_id)
            if job["status"] in ("checking", "preparing", "ready"):
                self._dispose(job, "cancelled")
            return self.status(job_id)

    def commit(self, job_id):
        game = self.game
        with game.lifecycle:
            with game.lock:
                job = self._job(job_id)
                if job["status"] == "committed":
                    if job["generation"] != game.generation:
                        raise ValueError("赛道已切换到另一辆赛车，请返回编辑。")
                    return self._public(job)
                if job["status"] != "ready":
                    raise ValueError(job["error"] or "赛车尚未准备完成。")
                reason = None
                if job["baseGeneration"] != game.generation:
                    reason = "赛道已切换，请重新准备。"
                elif job["baseline_settings"] != game.requested_settings or job["baseline_drive"] != game.requested_drive:
                    reason = "赛道设置已更新，请重新准备。"
                elif not job["process"].is_alive():
                    reason = "准备进程已退出，请重新准备。"
                if reason is None:
                    old_refs = game.process, game.commands, game.events
                    if game.state["molecule"] is not None and game.frames:
                        game.previous_replay = copy.deepcopy(game.replay())
                    game.generation += 1
                    game.input = job["input"]
                    game.competition = copy.deepcopy(job["competition"])
                    game.requested_settings = dict(job["settings"])
                    game.requested_drive = dict(job["drive"])
                    game.requested_start_pose = copy.deepcopy(job["pose"])
                    game.process, game.commands, game.events = job["process"], job["commands"], job["events"]
                    game.frames.clear()
                    event = job["ready"]
                    game.state.update({key: copy.deepcopy(event[key]) for key in
                                       ("molecule", "runtime", "frame", "control", "settings")})
                    game.state.update(status="paused", error=None, stats=None, finish_reason=None,
                                      generation=game.generation, start_pose=copy.deepcopy(job["pose"]),
                                      competition=copy.deepcopy(job["competition"]))
                    game.frames.append(game.state["frame"])
                    job.update(status="committed", generation=game.generation, ready=None,
                               process=None, commands=None, events=None)
                    threading.Thread(target=game._receive, args=(game.generation, game.process,
                                                               game.events), daemon=True).start()
                    result = self._public(job)
            if reason:
                self._dispose(job, "error", reason)
                raise ValueError(reason)
            stop_worker(*old_refs)
            return result

    def close_locked(self):
        for job in list(self.jobs.values()):
            if job["status"] in ("checking", "preparing", "ready"):
                self._dispose(job, "cancelled")
