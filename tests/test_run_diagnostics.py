import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

from gritic import diagnostics


def usage(user, system, peak):
    return SimpleNamespace(ru_utime=user, ru_stime=system, ru_maxrss=peak)


class RunDiagnosticsTest(unittest.TestCase):
    def test_cpu_is_interval_delta_but_memory_is_lifetime_peak(self):
        initial = {'self': usage(20, 4, 900), 'children': usage(40, 6, 1500)}
        final = {'self': usage(23, 5, 900), 'children': usage(47, 8, 1500)}
        with (
            mock.patch.object(diagnostics, '_usage', side_effect=[(initial, None), (final, None)]),
            mock.patch.object(diagnostics.time, 'perf_counter', side_effect=[100, 105]),
            mock.patch.object(diagnostics.sys, 'platform', 'linux'),
        ):
            report = diagnostics.RunDiagnostics(scope='test.operation').report('gritic')
        execution = report['execution']
        self.assertEqual(execution['elapsed_seconds'], 5)
        self.assertEqual(execution['cpu_seconds']['self'], {'user': 3, 'system': 1})
        self.assertEqual(execution['cpu_seconds']['children'], {'user': 7, 'system': 2})
        self.assertEqual(execution['cpu_seconds']['total'], {'user': 10, 'system': 3})
        self.assertEqual(
            execution['cpu_seconds']['children_scope'],
            'lifetime_of_children_reaped_during_interval',
        )
        peaks = execution['memory_peaks']
        self.assertEqual(peaks['self_peak_rss_bytes'], 900 * 1024)
        self.assertEqual(peaks['children_peak_rss_bytes'], 1500 * 1024)
        self.assertEqual(peaks['scope'], 'process_lifetime')
        self.assertFalse(peaks['includes_concurrent_process_tree'])
        self.assertEqual(peaks['children_aggregation'], 'largest_reaped_child')

    def test_macos_peak_is_already_bytes(self):
        snapshot = {'self': usage(0, 0, 1_000_000), 'children': usage(0, 0, 2_000_000)}
        with (
            mock.patch.object(diagnostics, '_usage', return_value=(snapshot, None)),
            mock.patch.object(diagnostics.sys, 'platform', 'darwin'),
        ):
            report = diagnostics.RunDiagnostics(scope='test').report('mutic')
        self.assertEqual(report['execution']['memory_peaks']['self_peak_rss_bytes'], 1_000_000)
        self.assertEqual(report['execution']['memory_peaks']['children_peak_rss_bytes'], 2_000_000)

    def test_unavailable_resource_counters_are_explicitly_null(self):
        with mock.patch.object(diagnostics, 'resource', None):
            report = diagnostics.RunDiagnostics(scope='test').report('plotic')
        self.assertIsNone(report['execution']['cpu_seconds']['total'])
        self.assertIsNone(report['execution']['memory_peaks']['self_peak_rss_bytes'])
        self.assertIn('unavailable', report['execution']['cpu_seconds']['unavailable_reason'])
        self.assertGreaterEqual(report['execution']['elapsed_seconds'], 0)

    def test_unrecognized_rss_units_do_not_change_cpu_accounting(self):
        snapshot = {'self': usage(0, 0, 100), 'children': usage(0, 0, 0)}
        with (
            mock.patch.object(diagnostics, '_usage', return_value=(snapshot, None)),
            mock.patch.object(diagnostics.sys, 'platform', 'unknown'),
        ):
            report = diagnostics.RunDiagnostics(scope='test').report('gritic')
        self.assertEqual(report['execution']['cpu_seconds']['self'], {'user': 0, 'system': 0})
        self.assertIsNone(report['execution']['memory_peaks']['self_peak_rss_bytes'])
        self.assertIn('units', report['execution']['memory_peaks']['unavailable_reason'])

    @unittest.skipIf(diagnostics.resource is None, 'requires getrusage')
    def test_cpu_includes_child_that_exits_before_report(self):
        tracker = diagnostics.RunDiagnostics(scope='test.child')
        subprocess.run([
            sys.executable, '-c',
            'import time; start = time.process_time(); '
            '\nwhile time.process_time() - start < 0.05: pass',
        ], check=True)
        report = tracker.report('gritic')
        child = report['execution']['cpu_seconds']['children']
        self.assertGreater(child['user'] + child['system'], 0.025)
        self.assertGreaterEqual(report['execution']['elapsed_seconds'], 0.025)

    def test_json_preserves_numpy_values_and_undefined_diagnostics(self):
        report = diagnostics.RunDiagnostics(scope='test').report(
            'gritic', sample_id='SAMPLE',
            inputs={'table': Path('input.tsv')},
            settings={'enabled': np.bool_(True), 'count': np.int64(7)},
            analysis={'values': np.array([1.0, np.nan, np.inf, -np.inf])},
        )
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'diagnostics.json'
            diagnostics.write_diagnostics(target, report)
            saved = json.loads(target.read_text())
            self.assertEqual(saved['settings'], {'enabled': True, 'count': 7})
            self.assertEqual(saved['analysis']['values'], [1.0, 'NaN', 'Infinity', '-Infinity'])
            self.assertEqual(saved['inputs']['table'], 'input.tsv')
            self.assertEqual(list(Path(directory).iterdir()), [target])

    def test_failed_atomic_install_preserves_existing_report(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'diagnostics.json'
            target.write_text('original')
            with mock.patch.object(diagnostics.os, 'replace', side_effect=OSError('failed install')):
                with self.assertRaisesRegex(OSError, 'failed install'):
                    diagnostics.write_diagnostics(target, {'status': 'completed'})
            self.assertEqual(target.read_text(), 'original')
            self.assertEqual(list(Path(directory).iterdir()), [target])

    def test_missing_distribution_version_is_null(self):
        with mock.patch.object(diagnostics.metadata, 'version', side_effect=diagnostics.metadata.PackageNotFoundError):
            versions = diagnostics.software_versions(package_names=('uninstalled',))
        self.assertEqual(versions['packages'], {'uninstalled': None})


if __name__ == '__main__':
    unittest.main()
