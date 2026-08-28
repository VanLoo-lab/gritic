import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pandas as pd

from gritic import gritictimer, sampletools


class ProcessSampleGuardTest(unittest.TestCase):
    @staticmethod
    def lightweight_sample(sample_id):
        return SimpleNamespace(
            sample_id=sample_id,
            sex=None,
            get_mutation_table=lambda: pd.DataFrame({
                'Mutation_ID': ['mutation'],
                'Segment_ID': ['segment'],
                'Segment_Mutation_Index': [0],
                'Phase_Group_ID': [0],
                'Chromosome': ['1'],
                'Major_CN': [1],
                'Minor_CN': [0],
                'Tumor_Ref_Count': [10],
                'Tumor_Alt_Count': [5],
                'Phasing': [pd.NA],
            }),
            get_count_group_table=lambda: pd.DataFrame({
                'Count_Group_ID': [0],
                'Tumor_Ref_Count': [10],
                'Tumor_Alt_Count': [5],
            }),
            get_phase_group_table=lambda: pd.DataFrame({
                'Phase_Group_ID': [0],
                'Count_Group_ID': [0],
                'Phasing': ['non_phased'],
            }),
            get_likelihood_context_table=lambda: pd.DataFrame({
                'Likelihood_Context_ID': [0],
                'Major_CN': [1],
                'Minor_CN': [0],
                'Normal_Total_CN': [2],
            }),
            get_segment_context_table=lambda: pd.DataFrame({
                'Segment_ID': ['segment'],
                'Likelihood_Context_ID': [0],
            }),
            get_count_group_likelihood_table=lambda: pd.DataFrame({
                'Likelihood_Context_ID': [0],
                'Count_Group_ID': [0],
                'Prob_Mult_1': [1.0],
            }),
            get_segment_group_table=lambda: pd.DataFrame({
                'Segment_ID': ['segment'],
                'Phase_Group_ID': [0],
                'N_Mutations': [1],
            }),
            get_subclone_table=lambda: None,
        )

    def test_interval_config_type_is_checked_before_sample_or_filesystem_access(self):
        with self.assertRaisesRegex(
            TypeError,
            'interval_config must be a TimingIntervalConfig',
        ):
            gritictimer.process_sample(
                object(),
                'unused-output',
                interval_config=object(),
            )

    def test_unsupported_modal_major_copy_number_fails_before_output_creation(self):
        sample = self.lightweight_sample('UNSUPPORTED')
        with tempfile.TemporaryDirectory() as directory:
            output_root = Path(directory)
            with mock.patch.object(
                gritictimer,
                'get_major_cn_mode',
                return_value=3,
            ):
                with self.assertRaisesRegex(
                    ValueError,
                    'supports only modal major copy numbers 1 and 2',
                ):
                    gritictimer.process_sample(sample, output_root)

            self.assertFalse((output_root / sample.sample_id).exists())

    def test_sample_output_path_that_is_a_file_is_rejected(self):
        sample = self.lightweight_sample('FILE')
        with tempfile.TemporaryDirectory() as directory:
            sample_output_path = Path(directory) / sample.sample_id
            sample_output_path.write_text('occupied', encoding='utf-8')

            with mock.patch.object(
                gritictimer,
                'get_major_cn_mode',
                return_value=1,
            ):
                with self.assertRaisesRegex(
                    FileExistsError,
                    'Sample output path is not a directory',
                ):
                    gritictimer.process_sample(sample, directory)

            self.assertEqual(
                sample_output_path.read_text(encoding='utf-8'),
                'occupied',
            )

    def test_preexisting_empty_sample_output_directory_is_accepted(self):
        sample = self.lightweight_sample('EMPTY')
        with tempfile.TemporaryDirectory() as directory:
            sample_output_path = Path(directory) / sample.sample_id
            sample_output_path.mkdir()

            with mock.patch.object(
                gritictimer,
                'get_major_cn_mode',
                return_value=1,
            ), mock.patch.object(gritictimer, '_run_sample') as run_sample:
                gritictimer.process_sample(
                    sample,
                    directory,
                    wgd_count=0,
                )

            self.assertTrue(
                (sample_output_path / 'EMPTY_mutation_table.tsv').is_file()
            )
            self.assertTrue(
                (
                    sample_output_path / 'EMPTY_count_group_table.tsv'
                ).is_file()
            )
            self.assertTrue(
                (
                    sample_output_path / 'EMPTY_phase_group_table.tsv'
                ).is_file()
            )
            self.assertTrue(
                (
                    sample_output_path
                    / 'EMPTY_likelihood_context_table.tsv'
                ).is_file()
            )
            self.assertTrue(
                (
                    sample_output_path / 'EMPTY_segment_context_table.tsv'
                ).is_file()
            )
            self.assertTrue(
                (
                    sample_output_path
                    / 'EMPTY_count_group_likelihood_table.tsv'
                ).is_file()
            )
            self.assertTrue(
                (
                    sample_output_path / 'EMPTY_segment_group_table.tsv'
                ).is_file()
            )
            self.assertTrue(
                (sample_output_path / 'EMPTY_subclone_table.tsv').is_file()
            )
            timing_directory = sample_output_path / 'EMPTY_timing_dicts'
            self.assertTrue(timing_directory.is_dir())
            run_sample.assert_called_once()
            self.assertEqual(run_sample.call_args.args[1], sample_output_path)
            self.assertEqual(
                Path(run_sample.call_args.args[2]),
                timing_directory,
            )


class SexChromosomeProcessIntegrationTest(unittest.TestCase):
    @staticmethod
    def make_xy_sample():
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1', 'X'],
            'Segment_Start': [0, 0],
            'Segment_End': [1_000, 1_000],
            'Major_CN': [1, 2],
            'Minor_CN': [1, 0],
        })
        mutation_table = pd.DataFrame({
            'Chromosome': ['1'] * 12 + ['X'] * 12,
            'Position': list(range(10, 22)) * 2,
            'Tumor_Ref_Count': [20] * 24,
            'Tumor_Alt_Count': [5] * 24,
        })
        return sampletools.Sample(
            mutation_table,
            copy_number_table,
            None,
            sample_id='XY_SAMPLE',
            purity=0.8,
            sex='XY',
            merge_cn=False,
        )

    def test_xy_gain_is_timed_and_written_through_process_sample(self):
        sample = self.make_xy_sample()
        with tempfile.TemporaryDirectory() as directory:
            gritictimer.process_sample(
                sample,
                directory,
                wgd_count=0,
                random_seed=20260831,
            )
            output_directory = Path(directory) / sample.sample_id

            route_table = pd.read_csv(
                output_directory / 'XY_SAMPLE_route_table.tsv',
                sep='\t',
            )
            gain_table = pd.read_csv(
                output_directory / 'XY_SAMPLE_gain_timing_table.tsv',
                sep='\t',
            )
            with (
                output_directory / 'XY_SAMPLE_wgd_calling_info.json'
            ).open(encoding='utf-8') as calling_info_file:
                calling_info = json.load(calling_info_file)

        x_routes = route_table.loc[route_table['Segment_ID'].eq('X-0-1000')]
        self.assertFalse(x_routes.empty)
        self.assertEqual(set(x_routes['Chromosome']), {'X'})
        x_gains = gain_table.loc[gain_table['Segment_ID'].eq('X-0-1000')]
        self.assertFalse(x_gains.empty)
        self.assertEqual(set(x_gains['Node_Phasing']), {'Major'})
        self.assertTrue(x_gains['Timing'].between(0, 1).all())
        self.assertEqual(calling_info['Major_CN_Mode'], 1)
        self.assertIs(calling_info['WGD_Status'], False)


if __name__ == '__main__':
    unittest.main()
