"""Serve the integrated browser site; no simulation API is provided."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

class Handler(SimpleHTTPRequestHandler):
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map,
                      '.js': 'text/javascript', '.mjs': 'text/javascript',
                      '.svg': 'image/svg+xml', '.ttf': 'font/ttf'}

    def end_headers(self):
        self.send_header('Cache-Control', 'no-store')
        super().end_headers()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=8796)
    args = parser.parse_args()
    directory = str(Path(__file__).resolve().parents[1] / 'wasm' / 'dist')
    server = ThreadingHTTPServer(('127.0.0.1', args.port), partial(Handler, directory=directory))
    print(f'Racing Molecule interface preview: http://127.0.0.1:{server.server_port}/', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
