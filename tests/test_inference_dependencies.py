import os
from pathlib import Path
import subprocess
import sys
import unittest


class InferenceDependencyTests(unittest.TestCase):
    def test_import_and_interval_calculation_do_not_import_plotting_libraries(self):
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
                    'from gritic.distributiontools import get_interval_bounds, IntervalSpec; '
                    'get_interval_bounds([0.1, 0.2, 0.8], IntervalSpec(0.9)); '
                    'assert not any(name.split(".")[0] in '
                    '{"matplotlib", "seaborn", "plotly"} for name in sys.modules)'
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
