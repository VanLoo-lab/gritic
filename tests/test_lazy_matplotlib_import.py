import os
from pathlib import Path
import subprocess
import sys
import unittest


class MatplotlibImportTests(unittest.TestCase):
    def test_gritictimer_import_does_not_import_matplotlib(self):
        project_root = Path(__file__).resolve().parents[1]
        environment = os.environ.copy()
        environment['PYTHONPATH'] = os.pathsep.join(
            filter(None, (str(project_root), environment.get('PYTHONPATH')))
        )
        result = subprocess.run(
            [
                sys.executable,
                '-c',
                (
                    'import sys; '
                    'import gritic.gritictimer; '
                    'assert "matplotlib" not in sys.modules'
                ),
            ],
            cwd=project_root,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
