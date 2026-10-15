"""Transactional preparation lifecycle, independently of slow physical MD."""
import copy
import queue
import unittest
from unittest.mock import Mock, patch
import uuid

from molecular_car.server import GameController
from molecular_car.editor.race_handoff import candidate_worker


class HandoffTests(unittest.TestCase):
    def setUp(self):
        self.game = GameController()
        self.game.process, self.game.commands, self.game.events = Mock(), Mock(), Mock()
        self.game.generation = 7
        self.game.input = ("original", "old.mol")
        self.game.state.update(status="running", generation=7, molecule={"name": "old", "atoms": [{"id": 1}]},
                               frame={"seq": 9, "step": 900}, runtime={"model": "test"})
        self.game.frames.extend([{"seq": 8, "step": 800}, {"seq": 9, "step": 900}])
        self.payload = {"id": str(uuid.uuid4()), "filename": "new.mol", "content": "candidate"}

    def begin(self, payload=None):
        self.context = Mock()
        self.context.Process.return_value.is_alive.return_value = True
        with patch("molecular_car.editor.race_handoff.mp.get_context", return_value=self.context), \
                patch("molecular_car.editor.race_handoff.threading.Thread"), \
                patch("molecular_car.editor.race_handoff.threading.Timer"):
            return self.game.handoff.begin(payload or self.payload)

    def ready(self):
        self.begin()
        self.job = self.game.handoff.jobs[self.payload["id"]]
        self.job.update(status="ready", atomCount=2, ready={"kind": "ready",
            "molecule": {"name": "new", "atoms": [{"id": 1}, {"id": 2}]},
            "frame": {"seq": 0, "step": 0}, "runtime": {"model": "new"},
            "settings": dict(self.game.requested_settings), "control": copy.deepcopy(self.game.state["control"])})
        return self.job

    def test_begin_keeps_live_race_and_replay(self):
        state, replay = self.game.snapshot(), copy.deepcopy(self.game.replay())
        old_worker = self.game.process
        result = self.begin()
        self.assertEqual(result["status"], "checking")
        self.assertEqual(self.game.snapshot(), state)
        self.assertEqual(self.game.replay(), replay)
        self.assertIs(self.game.process, old_worker)
        self.game.commands.put.assert_not_called()

    def test_prepare_retry_is_idempotent_and_conflicting_input_is_rejected(self):
        first = self.begin()
        self.assertEqual(self.game.handoff.begin(self.payload), first)
        self.context.Process.assert_called_once()
        with self.assertRaisesRegex(ValueError, "已改变"):
            self.game.handoff.begin({**self.payload, "content": "different"})

    def test_second_candidate_does_not_replace_first(self):
        self.begin()
        with self.assertRaisesRegex(ValueError, "另一辆"):
            self.game.handoff.begin({**self.payload, "id": str(uuid.uuid4())})
        self.assertEqual(len(self.game.handoff.jobs), 1)

    def test_cancel_stops_only_candidate_and_keeps_old_race(self):
        old = self.game.snapshot()
        self.ready()
        with patch("molecular_car.editor.race_handoff.stop_worker") as stop:
            cancelled = self.game.handoff.cancel(self.payload["id"])
        self.assertEqual(cancelled["status"], "cancelled")
        self.assertEqual(self.game.snapshot(), old)
        self.assertIs(stop.call_args.args[0], self.context.Process.return_value)
        self.assertIsNot(stop.call_args.args[0], self.game.process)
        self.assertEqual(self.game.handoff.cancel(self.payload["id"]), cancelled)

    def test_cancel_before_prepare_blocks_delayed_request(self):
        result = self.game.handoff.cancel(self.payload["id"])
        with patch("molecular_car.editor.race_handoff.mp.get_context") as context:
            self.assertEqual(self.game.handoff.begin(self.payload), result)
            context.assert_not_called()

    def test_success_promotes_same_worker_paused_and_archives_deep_copy(self):
        job = self.ready()
        candidate = job["process"]
        previous, old_worker = copy.deepcopy(self.game.replay()), self.game.process
        with patch("molecular_car.editor.race_handoff.stop_worker") as stop, \
                patch("molecular_car.editor.race_handoff.threading.Thread") as receiver:
            result = self.game.handoff.commit(job["id"])
        self.assertEqual(result["generation"], 8)
        self.assertEqual(self.game.state["status"], "paused")
        self.assertEqual(self.game.state["frame"]["step"], 0)
        self.assertEqual(list(self.game.frames), [{"seq": 0, "step": 0}])
        self.assertIs(self.game.process, candidate)
        self.assertIs(stop.call_args.args[0], old_worker)
        self.assertIs(receiver.call_args.kwargs["args"][1], candidate)
        self.assertEqual(self.game.previous_replay, previous)
        self.game.previous_replay["frames"][0]["step"] = -1
        self.assertEqual(previous["frames"][0]["step"], 800)
        self.assertEqual(self.game.snapshot()["previousReplay"]["name"], "old")

    def test_commit_retry_does_not_restart_or_archive_again(self):
        job = self.ready()
        with patch("molecular_car.editor.race_handoff.stop_worker"), patch("molecular_car.editor.race_handoff.threading.Thread"):
            result = self.game.handoff.commit(job["id"])
        with patch("molecular_car.editor.race_handoff.stop_worker") as stop:
            self.assertEqual(self.game.handoff.commit(job["id"]), result)
            self.assertEqual(self.game.handoff.cancel(job["id"]), result)
            stop.assert_not_called()
        self.assertEqual(self.game.generation, 8)

    def test_stale_generation_settings_and_drive_never_promote(self):
        for field in ("generation", "settings", "drive"):
            with self.subTest(field=field):
                self.setUp()
                job = self.ready()
                if field == "generation":
                    self.game.generation += 1
                elif field == "settings":
                    self.game.requested_settings["temperature"] += 1
                else:
                    self.game.requested_drive["strength"] += .01
                old = self.game.snapshot()
                with patch("molecular_car.editor.race_handoff.stop_worker"), self.assertRaises(ValueError):
                    self.game.handoff.commit(job["id"])
                self.assertEqual(self.game.snapshot(), old)
                self.assertEqual(job["status"], "error")
                self.assertIsNone(self.game.previous_replay)

    def test_dead_candidate_and_not_ready_are_rejected(self):
        self.begin()
        with self.assertRaisesRegex(ValueError, "尚未"):
            self.game.handoff.commit(self.payload["id"])
        self.game.handoff.jobs[self.payload["id"]]["status"] = "cancelled"
        self.game.handoff.jobs.clear()
        job = self.ready()
        job["process"].is_alive.return_value = False
        with patch("molecular_car.editor.race_handoff.stop_worker"), self.assertRaisesRegex(ValueError, "退出"):
            self.game.handoff.commit(job["id"])
        self.assertEqual(self.game.state["molecule"]["name"], "old")

    def test_worker_failure_is_separate_from_current_race(self):
        self.begin()
        job = self.game.handoff.jobs[self.payload["id"]]
        events = queue.Queue()
        events.put({"kind": "stage", "status": "preparing"})
        events.put({"kind": "error", "error": "joint relaxation failed"})
        job["events"] = events
        old = self.game.snapshot()
        with patch("molecular_car.editor.race_handoff.stop_worker"):
            self.game.handoff._receive(job)
        self.assertEqual(self.game.handoff.status(job["id"])["error"], "joint relaxation failed")
        self.assertEqual(self.game.snapshot(), old)

    def test_timeout_and_shutdown_dispose_candidates_only(self):
        job = self.ready()
        old = self.game.snapshot()
        with patch("molecular_car.editor.race_handoff.stop_worker") as stop:
            self.game.handoff._expire(job["id"], False)
            stop.assert_not_called()
            self.game.handoff._expire(job["id"], True)
            self.assertEqual(job["status"], "error")
            self.assertEqual(self.game.snapshot(), old)
        self.game.handoff.jobs.clear()
        job = self.ready()
        with patch("molecular_car.editor.race_handoff.stop_worker"), patch.object(self.game, "_stop"):
            self.game.close()
        self.assertEqual(job["status"], "cancelled")

    def test_private_progress_hook_reuses_exact_prepared_object_once(self):
        from molecular_car import engine
        events, prepared = Mock(), object()
        prepare = Mock(return_value=prepared)
        def worker(*args):
            self.assertIs(engine.prepare_molecule("m", "f"), prepared)
        with patch.object(engine, "prepare_molecule", prepare), patch.object(engine, "worker_main", worker):
            candidate_worker("m", "f", {}, Mock(), events, {}, {})
            self.assertIs(engine.prepare_molecule, prepare)
        prepare.assert_called_once()
        self.assertEqual([call.args[0]["status"] for call in events.put.call_args_list], ["checking", "preparing"])


if __name__ == "__main__":
    unittest.main()
