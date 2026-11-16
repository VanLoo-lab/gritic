"""Common provenance and resource measurements for completed suite runs.

CPU counters cover the measured interval; RSS counters are process-lifetime
high-water marks. Reaped-child RSS is a largest-process peak, never a sum of
simultaneously resident memory. Repeated calls in one interpreter retain earlier
memory peaks. Keep these distinctions in the saved report as well as the docs.
"""

from collections.abc import Mapping
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from importlib import metadata
import json
import math
from numbers import Integral, Real
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

try:
    import resource
except ImportError:
    resource = None


_PACKAGES = (
    'gritic', 'mutic', 'sigtic', 'plotic', 'hopsy', 'numpy', 'pandas',
    'scipy', 'numba', 'arviz-stats', 'matplotlib', 'av',
)


def software_versions(package_names=None):
    """Describe the interpreter and installed distributions without importing them."""
    packages = {}
    for name in _PACKAGES if package_names is None else package_names:
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            packages[name] = None
    return {
        'python': platform.python_version(),
        'platform': platform.platform(),
        'packages': packages,
    }


def _json_values(value):
    """Retain undefined diagnostics as strings in standards-compliant JSON."""
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return _json_values(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_values(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_values(item) for item in value]
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        value = float(value)
        if math.isfinite(value):
            return value
        return 'NaN' if math.isnan(value) else ('Infinity' if value > 0 else '-Infinity')
    # NumPy arrays/scalars expose tolist(); no numerical imports are required
    # by this module, including when a CLI imports it to provide --help.
    if hasattr(value, 'tolist'):
        return _json_values(value.tolist())
    raise TypeError(f'Unsupported diagnostic value: {type(value).__name__}')


def write_diagnostics(path, report):
    """Atomically publish one JSON file; callers own output collision policy."""
    path = Path(path)
    encoded = json.dumps(_json_values(report), indent=2, allow_nan=False) + '\n'
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode='w', encoding='utf-8', dir=path.parent,
            prefix=f'.{path.name}.', suffix='.tmp', delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(encoded)
        os.replace(temporary_path, path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return path


def _usage():
    if resource is None:
        return None, 'getrusage is unavailable on this platform'
    try:
        return {
            'self': resource.getrusage(resource.RUSAGE_SELF),
            'children': resource.getrusage(resource.RUSAGE_CHILDREN),
        }, None
    except OSError as error:
        return None, f'getrusage failed: {error}'


def _rss_multiplier():
    if sys.platform == 'darwin':
        return 1  # Darwin reports bytes.
    if sys.platform.startswith('linux'):
        return 1024  # Linux reports KiB.
    return None


class RunDiagnostics:
    """Measure one synchronous operation, without polling or a new dependency.

Elapsed/CPU measurements end when report() is called, excluding report creation
and any later output publication. CPU counters include all threads of the
    current process and the lifetime CPU of children reaped during the interval.
    A child started before the interval can therefore contribute earlier CPU time.
    Counters do not cover still-running children or isolate concurrent work in
    one process.
"""

    def __init__(self, *, scope):
        self.scope = scope
        self.started_at = datetime.now(timezone.utc).isoformat()
        self._started = time.perf_counter()
        self._initial, self._initial_error = _usage()

    def report(
        self, package, *, sample_id=None, command=None,
        inputs=None, settings=None, outputs=None, analysis=None,
    ):
        """Build a completion record to publish with the operation's outputs."""
        current, current_error = _usage()
        elapsed = time.perf_counter() - self._started
        cpu = {
            'method': 'getrusage',
            'scope': 'interval',
            'children_scope': 'lifetime_of_children_reaped_during_interval',
            'self': None,
            'children': None,
            'total': None,
            'unavailable_reason': self._initial_error or current_error,
        }
        if self._initial is not None and current is not None:
            for who in ('self', 'children'):
                cpu[who] = {
                    'user': current[who].ru_utime - self._initial[who].ru_utime,
                    'system': current[who].ru_stime - self._initial[who].ru_stime,
                }
            cpu['total'] = {
                mode: cpu['self'][mode] + cpu['children'][mode]
                for mode in ('user', 'system')
            }
        multiplier = _rss_multiplier()
        memory = {
            'method': 'getrusage.ru_maxrss',
            'scope': 'process_lifetime',
            'self_peak_rss_bytes': None,
            'children_peak_rss_bytes': None,
            'children_aggregation': 'largest_reaped_child',
            'includes_concurrent_process_tree': False,
            'unavailable_reason': current_error,
        }
        if current is not None:
            if multiplier is None:
                memory['unavailable_reason'] = 'peak RSS units are unsupported on this platform'
            else:
                for who in ('self', 'children'):
                    memory[f'{who}_peak_rss_bytes'] = int(current[who].ru_maxrss) * multiplier
        report = {
            'format': f'{package}_run_diagnostics',
            'version': 1,
            'status': 'completed',
            'outputs_written': True,
            'inputs': {} if inputs is None else inputs,
            'settings': {} if settings is None else settings,
            'outputs': [] if outputs is None else outputs,
            'software': software_versions(),
            'execution': {
                'scope': self.scope,
                'started_at': self.started_at,
                'recorded_at': datetime.now(timezone.utc).isoformat(),
                'elapsed_seconds': elapsed,
                'cpu_seconds': cpu,
                'memory_peaks': memory,
            },
        }
        if sample_id is not None:
            report['sample_id'] = str(sample_id)
        if command is not None:
            report['command'] = command
        if analysis is not None:
            report['analysis'] = analysis
        return _json_values(report)
