import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from tests._plot_cache import MPL_CONFIG_DIRECTORY, XDG_CACHE_DIRECTORY


class CliModulePlotIntegrationTest(unittest.TestCase):
    @staticmethod
    def write_input_tables(directory):
        copy_number_path = directory / 'copy-number.tsv'
        mutation_path = directory / 'mutations.tsv'
        pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Segment_Start': [0, 1_000],
            'Segment_End': [1_000, 1_200],
            'Major_CN': [1, 2],
            'Minor_CN': [1, 1],
        }).to_csv(copy_number_path, sep='\t', index=False)
        pd.DataFrame({
            'Chromosome': ['1'] * 24,
            'Position': (
                list(range(10, 22)) + list(range(1_010, 1_022))
            ),
            'Tumor_Ref_Count': [20] * 24,
            'Tumor_Alt_Count': [5] * 24,
        }).to_csv(mutation_path, sep='\t', index=False)
        return copy_number_path, mutation_path

    def test_python_module_guard_runs_real_argv_to_outputs_and_tree_pdfs(self):
        project_root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            copy_number_path, mutation_path = self.write_input_tables(
                temporary_path
            )
            output_root = temporary_path / 'output'
            environment = os.environ.copy()
            environment.update({
                'MPLBACKEND': 'Agg',
                'MPLCONFIGDIR': str(MPL_CONFIG_DIRECTORY),
                'XDG_CACHE_HOME': str(XDG_CACHE_DIRECTORY),
                'PYTHONHASHSEED': '0',
            })

            result = subprocess.run(
                [
                    sys.executable,
                    '-m',
                    'gritic.cli',
                    '--mutation-table',
                    str(mutation_path),
                    '--copy-number-table',
                    str(copy_number_path),
                    '--purity',
                    '0.8',
                    '--sample-id',
                    'CLI_PLOT',
                    '--output',
                    str(output_root),
                    '--sample-sex',
                    'XX',
                    '--wgd-count',
                    '0',
                    '--random-seed',
                    '20260828',
                    '--no-merge-adjacent-segments',
                    '--plot-trees',
                ],
                cwd=project_root,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
                timeout=180,
            )

            self.assertEqual(
                result.returncode,
                0,
                f'stdout:\n{result.stdout}\nstderr:\n{result.stderr}',
            )
            sample_output = output_root / 'CLI_PLOT'
            route_path = sample_output / 'CLI_PLOT_route_table.tsv'
            gain_path = sample_output / 'CLI_PLOT_gain_timing_table.tsv'
            self.assertTrue(route_path.is_file())
            self.assertTrue(gain_path.is_file())

            route_table = pd.read_csv(route_path, sep='\t')
            gain_table = pd.read_csv(gain_path, sep='\t')
            self.assertEqual(
                set(route_table['Segment_ID']),
                {'1-0-1000', '1-1000-1200'},
            )
            self.assertEqual(set(gain_table['Segment_ID']), {'1-1000-1200'})

            plot_directory = (
                sample_output
                / 'CLI_PLOT_tree_plots'
                / '1-1000-1200'
            )
            plot_paths = sorted(plot_directory.glob('route_*.pdf'))
            self.assertTrue(plot_paths)
            for plot_path in plot_paths:
                with self.subTest(plot=plot_path.name):
                    self.assertGreater(plot_path.stat().st_size, 1_000)
                    self.assertEqual(plot_path.read_bytes()[:4], b'%PDF')


if __name__ == '__main__':
    unittest.main()
