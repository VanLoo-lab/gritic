"""Biological exclusions have a distinct CLI contract from invalid inputs."""
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from gritic import cli, gritictimer
from gritic.outcomes import EXIT_INELIGIBLE, IneligibleSampleError


def arguments(tmp_path, *, alt_count=1):
    mutations = tmp_path / 'mutations.tsv'
    copy_number = tmp_path / 'copy-number.tsv'
    pd.DataFrame({
        'Chromosome': ['1'], 'Position': [10], 'Mutation_ID': ['m1'],
        'Tumor_Ref_Count': [20], 'Tumor_Alt_Count': [alt_count],
    }).to_csv(mutations, sep='\t', index=False)
    pd.DataFrame({
        'Chromosome': ['1'], 'Segment_Start': [0], 'Segment_End': [100],
        'Major_CN': [1], 'Minor_CN': [1],
    }).to_csv(copy_number, sep='\t', index=False)
    return ['--sample-id', 'SA1', '--sample-dir', str(tmp_path / 'output'),
            '--purity', '.8', '--mutation-table', str(mutations), '--copy-number-table', str(copy_number)]


def test_real_cli_returns_reserved_status_for_biological_exclusion(tmp_path):
    args = arguments(tmp_path)
    env = {**os.environ, 'PYTHONPATH': str(Path(cli.__file__).resolve().parents[1])}
    result = subprocess.run([sys.executable, '-m', 'gritic.cli', *args],
                            cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == EXIT_INELIGIBLE == 3
    outcome = json.loads(result.stderr.splitlines()[-1])
    assert outcome == {
        'format': 'gritic-outcome', 'status': 'ineligible', 'sample_id': 'SA1',
        'reason': 'no_mutations_after_read_filters',
        'message': 'No mutations remain after applying the configured minimum alternate-read count and coverage thresholds',
    }
    assert not (tmp_path / 'output').exists()
    assert 'Traceback' not in result.stderr


def test_cli_invalid_input_is_not_reclassified_as_biological_exclusion(tmp_path, capsys):
    args = arguments(tmp_path, alt_count=-1)
    with pytest.raises(ValueError) as error:
        cli.main(args)
    assert not isinstance(error.value, IneligibleSampleError)
    assert 'gritic-outcome' not in capsys.readouterr().err
    assert not (tmp_path / 'output').exists()


def test_cli_operational_and_numerical_errors_propagate(tmp_path, monkeypatch, capsys):
    args = arguments(tmp_path)
    for error in (OSError('storage unavailable'), ValueError('nonfinite timing posterior')):
        def fail(*args):
            raise error
        monkeypatch.setattr(cli, '_run_analysis', fail)
        with pytest.raises(type(error), match=str(error)):
            cli.main(args)
    assert 'gritic-outcome' not in capsys.readouterr().err


def test_unsupported_modal_copy_number_is_structured_before_output_creation(tmp_path, monkeypatch):
    monkeypatch.setattr(gritictimer, 'get_major_cn_mode', lambda sample: 3)
    with pytest.raises(IneligibleSampleError) as error:
        gritictimer.process_sample(SimpleNamespace(sample_id='SA1'), tmp_path / 'output')
    assert error.value.reason == 'unsupported_major_copy_number_mode'
    assert not (tmp_path / 'output').exists()


def test_missing_wgd_segments_is_ineligible(tmp_path, monkeypatch):
    sample = SimpleNamespace(sample_id='SA1')
    monkeypatch.setattr(gritictimer, 'get_potential_wgd_segments', lambda sample: [])
    with pytest.raises(IneligibleSampleError) as error:
        gritictimer.time_wgd_major_cn_2(sample, tmp_path, tmp_path / 'timing', rng=np.random.default_rng(1))
    assert error.value.reason == 'no_eligible_wgd_segments'


def test_nonfinite_wgd_posterior_is_not_a_biological_exclusion(tmp_path, monkeypatch):
    sample, segment = SimpleNamespace(sample_id='SA1'), SimpleNamespace(major_cn=2)
    monkeypatch.setattr(gritictimer, 'get_potential_wgd_segments', lambda sample: [segment])
    monkeypatch.setattr(gritictimer, '_get_wgd_timing_model', lambda segment: (0, None))
    monkeypatch.setattr(gritictimer, 'RouteClassifier', lambda *args, **kwargs: SimpleNamespace(routes={}))
    monkeypatch.setattr(gritictimer, 'fit_route_classifiers', lambda *args: None)
    monkeypatch.setattr(gritictimer, '_get_wgd_segment_result', lambda *args: (pd.DataFrame(), np.array([np.nan])))
    with pytest.raises(ValueError, match='finite WGD') as error:
        gritictimer.time_wgd_major_cn_2(sample, tmp_path, tmp_path / 'timing', rng=np.random.default_rng(1))
    assert not isinstance(error.value, IneligibleSampleError)


def test_segment_mutation_threshold_uses_structured_exclusion(tmp_path, capsys):
    args = arguments(tmp_path, alt_count=5)
    assert cli.main(args) == EXIT_INELIGIBLE
    outcome = json.loads(capsys.readouterr().err.splitlines()[-1])
    assert outcome['reason'] == 'no_segments_meeting_mutation_count'
