"""Local static preview only. This server has no molecular calculation API."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class StaticHandler(SimpleHTTPRequestHandler):
    # Windows registry associations can report .mjs as text/plain. ES modules
    # and streaming WebAssembly require explicit, portable content types.
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map,
                      ".js": "text/javascript", ".mjs": "text/javascript",
                      ".wasm": "application/wasm", ".gz": "application/gzip"}

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8783)
    args = parser.parse_args()
    handler = partial(StaticHandler, directory=str(ROOT / "wasm" / "dist"))
    server = ThreadingHTTPServer(("127.0.0.1", args.port), handler)
    print(f"Racing Molecule browser preview: http://127.0.0.1:{server.server_port}/", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
