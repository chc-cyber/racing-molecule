"""Playback scheduling must preserve complete, ordered MD samples."""
import queue
import threading
import unittest
from unittest.mock import patch
from molecular_car.engine import worker_main
from molecular_car.server import GameController


class FakeSimulation:
    molecule = runtime = settings = drive = {}
    def __init__(self, *args, **kwargs): self.step = 0
    def frame(self): return dict(seq=self.step // 50, step=self.step, time_fs=self.step, progress=self.step / 1000)
    def advance(self): self.step += 50; return self.frame()
    def close(self): pass


class PlaybackBufferTest(unittest.TestCase):
    def test_worker_window_waits_without_skipping_steps(self):
        commands, events = queue.Queue(), queue.Queue()
        commands.put(dict(action='playback', until_fs=100))
        commands.put(dict(action='start'))
        with patch('molecular_car.engine.Simulation', FakeSimulation):
            thread = threading.Thread(target=worker_main, args=('', '', {}, commands, events), daemon=True)
            thread.start()
            try:
                frames=[]
                while len(frames)<2:
                    event=events.get(timeout=2)
                    if event['kind']=='frame': frames.append(event['frame'])
                with self.assertRaises(queue.Empty): events.get(timeout=.1)
                commands.put(dict(action='playback', until_fs=1000))
                while True:
                    event=events.get(timeout=2)
                    if event['kind']=='frame': frames.append(event['frame'])
                    if event.get('status')=='finished': break
                self.assertEqual([f['step'] for f in frames], list(range(50,1001,50)))
            finally:
                commands.put(dict(action='stop'));thread.join(timeout=2)
            self.assertFalse(thread.is_alive())

    def test_stream_pages_from_first_unread_and_ignores_stale_credit(self):
        game=GameController();game.commands=queue.Queue()
        game.frames.extend(dict(seq=n) for n in range(600))
        first=game.stream(-1, 50, game.generation)
        self.assertEqual([f['seq'] for f in first['frames']],list(range(250)))
        self.assertFalse(first['dropped'])
        self.assertEqual(game.commands.get_nowait()['until_fs'],10050)
        second=game.stream(249,100,game.generation-1)
        self.assertEqual(second['frames'][0]['seq'],250)
        self.assertTrue(game.commands.empty())


if __name__ == '__main__': unittest.main()
