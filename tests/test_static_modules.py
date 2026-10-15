"""Exercise the production file handler with Windows' incorrect MIME mapping."""
import io
import unittest
from unittest.mock import Mock, patch

from molecular_car.server import Handler, ROOT


class StaticModuleTests(unittest.TestCase):
    def test_js_modules_override_platform_plain_text_mapping(self):
        for path in ('/app.js', '/race-course.mjs', '/race-course.mjs?v=reload',
                     '/browser/worker.mjs'):
            with self.subTest(path=path):
                handler = Handler.__new__(Handler)
                handler.path = path
                handler.wfile = io.BytesIO()
                handler.send_response = Mock()
                handler.send_header = Mock()
                handler.end_headers = Mock()
                with patch('molecular_car.server.mimetypes.guess_type', return_value=('text/plain', None)):
                    handler.do_GET()
                handler.send_response.assert_called_once_with(200)
                headers = dict(call.args for call in handler.send_header.call_args_list)
                self.assertEqual(headers['Content-Type'], 'text/javascript; charset=utf-8')
                self.assertEqual(handler.wfile.getvalue(),
                    (ROOT/'frontend'/path.split('?')[0].lstrip('/')).read_bytes())


if __name__ == '__main__':
    unittest.main()
