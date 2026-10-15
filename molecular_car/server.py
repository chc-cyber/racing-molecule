"""Same-origin local HTTP UI and one isolated simulation worker."""

import argparse
from collections import deque
import copy
from functools import partial
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import math
import mimetypes
import multiprocessing as mp
import os
from pathlib import Path
import queue
import socket
import threading
from urllib.parse import parse_qs, unquote, urlsplit

from .engine import (DEFAULT_SETTINGS, DEFAULT_DRIVE, TRACK, TIMESTEP_FS,
                     worker_main, validate_surface_mode)
from .start_pose import validate_start_pose
from .competition import validate_competition, engine_parameters

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = {"tetraphenylbenzene.mol", "benzene.mol"}


class LocalHTTPServer(ThreadingHTTPServer):
    # Windows SO_REUSEADDR can let a second instance hijack an occupied port.
    # Exclusive binding makes accidental duplicate starts fail immediately.
    allow_reuse_address = os.name != "nt"

    def server_bind(self):
        if os.name == "nt":
            self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        super().server_bind()


class GameController:
    def __init__(self):
        self.lock = threading.RLock()
        self.lifecycle = threading.Lock()
        self.process = None
        self.commands = None
        self.events = None
        self.generation = 0
        self.input = None
        # Accepted controls survive an immediate reset, before worker acknowledgement.
        self.competition = validate_competition()
        self.requested_drive = dict(DEFAULT_DRIVE)
        self.requested_settings = dict(DEFAULT_SETTINGS)
        self.requested_start_pose = validate_start_pose()
        self.frames = deque(maxlen=5010)
        self.previous_replay = None
        from .editor.race_handoff import RaceHandoff
        self.handoff = RaceHandoff(self)
        self.state = {"status": "idle", "error": None, "molecule": None, "frame": None,
                      "track": dict(TRACK), "settings": dict(DEFAULT_SETTINGS),
                      "control": {"drive": dict(DEFAULT_DRIVE)}, "competition": copy.deepcopy(self.competition),
                      "runtime": None, "stats": None, "generation": 0,
                      "start_pose": copy.deepcopy(self.requested_start_pose)}

    def snapshot(self):
        with self.lock:
            result = copy.deepcopy(self.state)
            previous = self.previous_replay
            result["previousReplay"] = ({"name": previous["molecule"]["name"],
                "frames": len(previous["frames"]), "truncated": previous["truncated"]} if previous else None)
            return result

    def _stop(self):
        if self.process is not None:
            if self.process.is_alive():
                self.commands.put({"action": "stop"})
                self.process.join(timeout=2)
            if self.process.is_alive():
                self.process.terminate()
                self.process.join(timeout=2)
            self.commands.close()
            self.events.close()
            self.process = None

    def close(self):
        with self.lifecycle:
            self.handoff.close_locked()
            self.generation += 1
            self._stop()

    def load(self, content, filename, *, start_pose=None, competition=None):
        if not isinstance(content, str) or not content.strip():
            raise ValueError("分子文件不能为空。")
        if not isinstance(filename, str) or Path(filename).suffix.lower() not in (".mol", ".sdf"):
            raise ValueError("首版支持 MOL 或单分子 SDF 文件。")
        if len(content.encode("utf-8")) > 2_000_000:
            raise ValueError("文件超过 2 MB。")
        rules = validate_competition(competition)
        if start_pose is not None:
            rules = validate_competition(dict(rules, start_pose=start_pose))
        with self.lifecycle:
            self.competition = rules
            self.requested_settings, self.requested_drive = engine_parameters(rules)
            self.requested_start_pose = copy.deepcopy(rules["start_pose"])
            self._load_locked(content, filename)

    def _load_locked(self, content, filename):
        """Start a fresh worker while the caller owns the lifecycle lock."""
        self.generation += 1
        generation = self.generation
        self._stop()
        self.input = (content, Path(filename).name)
        with self.lock:
            settings = dict(self.requested_settings)
            drive = dict(self.requested_drive)
            pose = copy.deepcopy(self.requested_start_pose)
            self.frames.clear()
            self.state.update(status="preparing", error=None, molecule=None, frame=None,
                              runtime=None, stats=None, generation=generation,
                              finish_reason=None, settings=settings, start_pose=pose,
                              control={"drive": drive}, competition=copy.deepcopy(self.competition))
        context = mp.get_context("spawn")
        self.commands = context.Queue()
        self.events = context.Queue(maxsize=32)
        self.process = context.Process(target=worker_main,
            args=(content, filename, settings, self.commands, self.events, drive, pose), daemon=True)
        self.process.start()
        threading.Thread(target=self._receive, args=(generation, self.process, self.events), daemon=True).start()

    def _receive(self, generation, process, events):
        while generation == self.generation:
            try:
                event = events.get(timeout=0.2)
            except queue.Empty:
                if not process.is_alive():
                    with self.lock:
                        if generation == self.generation and self.state["status"] != "error":
                            self.state.update(status="error", error=f"模拟进程意外退出（{process.exitcode}），请重置。")
                    return
                continue
            except (OSError, ValueError, EOFError):
                return
            with self.lock:
                if generation != self.generation:
                    return
                kind = event["kind"]
                if kind == "ready":
                    for key in ("molecule", "runtime", "frame", "control", "settings"):
                        self.state[key] = event[key]
                    self.state["status"] = "paused"
                    self.frames.append(event["frame"])
                elif kind == "frame":
                    self.state["frame"] = event["frame"]
                    self.state["stats"] = event["stats"]
                    self.state["control"]["drive"] = event["frame"]["drive"]
                    self.frames.append(event["frame"])
                elif kind == "control":
                    self.state["control"] = event["control"]
                    self.state["settings"] = event["settings"]
                elif kind == "status":
                    self.state["status"] = event["status"]
                    if "reason" in event:
                        self.state["finish_reason"] = event["reason"]
                elif kind == "error":
                    self.state.update(status="error", error=event["error"])

    def action(self, action):
        if action == "reset":
            with self.lifecycle:
                if not self.input:
                    raise ValueError("请先加载分子。")
                self._load_locked(*self.input)
            return
        if action not in ("start", "pause"):
            raise ValueError("未知的比赛操作。")
        with self.lock:
            if self.state["status"] not in ("paused", "running"):
                raise ValueError("请等待分子准备完成；比赛结束后请重置。")
            self.commands.put({"action": action})
            self.state["status"] = "running" if action == "start" else "paused"

    def start_pose(self, data):
        pose = validate_start_pose(data)
        with self.lifecycle:
            with self.lock:
                if not self.input or self.state["status"] in ("preparing", "running", "idle"):
                    raise ValueError("请等待准备完成或暂停比赛，再修改起跑姿态。")
                if self.state.get("molecule"):
                    validate_start_pose(pose, atom_count=len(self.state["molecule"]["atoms"]))
                rules = validate_competition(dict(self.competition, start_pose=pose))
                self.competition = rules
                self.requested_start_pose = pose
                content, filename = self.input
            self._load_locked(content, filename)

    def control(self, data):
        raise ValueError("赛事参数在准备时固定，请选择赛事配置后重新准备。")

    def stream(self, after, playback_time_fs=None, generation=None):
        with self.lock:
            if playback_time_fs is not None:
                if not math.isfinite(playback_time_fs) or playback_time_fs < 0:
                    raise ValueError("Invalid playback time")
                if generation == self.generation and self.commands is not None:
                    # At most 10 ps ahead; old API clients can still run uncapped.
                    self.commands.put({"action": "playback", "until_fs": playback_time_fs + 10000})
            latest = self.frames[-1]["seq"] if self.frames else -1
            earliest = self.frames[0]["seq"] if self.frames else 0
            dropped = after > latest or after < earliest - 1
            frames = list(self.frames)[:250] if dropped else [f for f in self.frames if f["seq"] > after][:250]
            return {"frames": frames, "latest_seq": latest, "status": self.state["status"],
                    "dropped": dropped, "generation": self.generation}

    def replay(self):
        with self.lock:
            return {"format": "molecular-car-replay-v1", "molecule": self.state["molecule"],
                    "track": self.state["track"], "settings": self.state["settings"],
                    "competition": copy.deepcopy(self.competition),
                    "start_pose": copy.deepcopy(self.state["start_pose"]),
                    "runtime": self.state["runtime"], "timestep_fs": TIMESTEP_FS,
                    "surface": (self.state["runtime"] or {}).get("surface"),
                    "control": self.state["control"],
                    "model": (self.state["runtime"] or {}).get("model", "uff-gradient-v3"),
                    "finish_reason": self.state.get("finish_reason"),
                    "truncated": bool(self.frames and self.frames[0]["seq"] > 0),
                    "frames": list(self.frames)}


class Handler(BaseHTTPRequestHandler):
    def __init__(self, *args, game, **kwargs):
        self.game = game
        super().__init__(*args, **kwargs)

    def log_message(self, format, *args):
        # Continuous polling should not fill the user's terminal.
        if args and str(args[1] if len(args) > 1 else "") not in ("200", "304"):
            super().log_message(format, *args)

    def _json(self, payload, status=200, attachment=None):
        content = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store")
        if attachment:
            self.send_header("Content-Disposition", f'attachment; filename="{attachment}"')
        self.end_headers()
        self.wfile.write(content)

    def do_GET(self):
        url = urlsplit(self.path)
        if url.path == "/api/state":
            self._json(self.game.snapshot())
        elif url.path == "/api/frames":
            try:
                query = parse_qs(url.query)
                after = int(query.get("after", ["-1"])[0])
                played = float(query["played"][0]) if "played" in query else None
                generation = int(query["generation"][0]) if "generation" in query else None
                self._json(self.game.stream(after, played, generation))
            except ValueError:
                self._json({"error": "帧编号格式不正确。"}, 400)
        elif url.path == "/api/replay":
            download = parse_qs(url.query).get("download", ["0"])[0] == "1"
            self._json(self.game.replay(), attachment="molecular-car-replay.json" if download else None)
        elif url.path == "/api/replay/previous":
            with self.game.lock:
                previous = copy.deepcopy(self.game.previous_replay)
            if previous is None:
                self._json({"error": "尚无上一场轨迹。"}, 404)
            else:
                self._json(previous, attachment="molecular-car-previous-replay.json")
        elif url.path == "/api/editor/race/status":
            try:
                self._json(self.game.handoff.status(parse_qs(url.query).get("id", [None])[0]))
            except ValueError as error:
                self._json({"error": str(error)}, 400)
        elif url.path.startswith("/api/"):
            self._json({"error": "接口不存在。"}, 404)
        else:
            frontend = (ROOT / "frontend").resolve()
            relative = unquote(url.path).lstrip("/") or "index.html"
            selected = (frontend / relative).resolve()
            if not selected.is_relative_to(frontend) or not selected.is_file():
                self.send_error(404)
                return
            content = selected.read_bytes()
            # Windows registry associations may label .mjs as text/plain.
            # ES-module imports require a JavaScript MIME type for both suffixes.
            mime = "text/javascript" if selected.suffix.lower() in (".js", ".mjs") else mimetypes.guess_type(str(selected))[0] or "application/octet-stream"
            self.send_response(200)
            self.send_header("Content-Type", mime + ("; charset=utf-8" if mime.startswith("text/") else ""))
            self.send_header("Content-Length", str(len(content)))
            self.send_header("Cache-Control", "no-cache")
            self.end_headers()
            self.wfile.write(content)

    def do_POST(self):
        try:
            origin = self.headers.get("Origin")
            port = self.server.server_port
            if origin and origin not in (f"http://127.0.0.1:{port}", f"http://localhost:{port}"):
                raise ValueError("请从本机游戏页面发起操作。")
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 3_000_000:
                raise ValueError("请求为空或过大。")
            data = json.loads(self.rfile.read(length))
            if not isinstance(data, dict):
                raise ValueError("请求格式不正确。")
            path = urlsplit(self.path).path
            if path.startswith("/api/editor/race/"):
                if path == "/api/editor/race/prepare":
                    result = self.game.handoff.begin(data)
                elif path == "/api/editor/race/cancel":
                    result = self.game.handoff.cancel(data.get("id"))
                elif path == "/api/editor/race/commit":
                    result = self.game.handoff.commit(data.get("id"))
                else:
                    raise ValueError("未知的赛车准备操作。")
                self._json(result)
                return
            if path.startswith("/api/editor/"):
                from .editor import dispatch
                self._json(dispatch(path, data))
                return
            if path == "/api/load":
                if data.get("reuse"):
                    if not self.game.input:
                        raise ValueError("请先选择分子。")
                    self.game.load(*self.game.input, competition=data.get("competition"))
                elif "example" in data:
                    if data["example"] not in EXAMPLES:
                        raise ValueError("未知的示例分子。")
                    selected = ROOT / "examples" / data["example"]
                    self.game.load(selected.read_text(encoding="utf-8"), selected.name, start_pose=data.get("start_pose"), competition=data.get("competition"))
                else:
                    self.game.load(data.get("content"), data.get("filename"), start_pose=data.get("start_pose"), competition=data.get("competition"))
            elif path == "/api/start-pose":
                self.game.start_pose(data)
            elif path == "/api/action":
                self.game.action(data.get("action"))
            elif path == "/api/control":
                self.game.control(data)
            elif path == "/api/shutdown":
                self._json({"status": "stopping"})
                # shutdown must run off the HTTP serving thread; request handlers
                # already have their own threads in ThreadingHTTPServer.
                threading.Thread(target=self.server.shutdown, daemon=True).start()
                return
            else:
                self._json({"error": "接口不存在。"}, 404)
                return
            self._json(self.game.snapshot())
        except (ValueError, TypeError, KeyError, UnicodeError) as error:
            self._json({"error": str(error)}, 400)


def main():
    parser = argparse.ArgumentParser(description="Local LAMMPS molecular racing prototype")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--molecule", type=Path, help="Developer supplied MOL/SDF file")
    args = parser.parse_args()
    if not (1 <= args.port <= 65535):
        parser.error("port must be in 1..65535")
    if not (ROOT / "frontend" / "vendor" / "three.module.js").is_file():
        parser.error("Run scripts/fetch_frontend.py to install the browser assets.")
    game = GameController()
    if args.molecule:
        game.input = (args.molecule.read_text(encoding="utf-8-sig"), args.molecule.name)
        game.state["selected_input"] = args.molecule.name
    try:
        server = LocalHTTPServer(("127.0.0.1", args.port), partial(Handler, game=game))
    except OSError as error:
        parser.error(f"Port {args.port} is unavailable. Open an existing local game or choose --port. ({error})")
    print(f"Racing Molecule: http://127.0.0.1:{args.port}", flush=True)
    print("Choose race rules and prepare a molecule in the browser; Ctrl+C stops the local server.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        game.close()


if __name__ == "__main__":
    mp.freeze_support()
    main()
