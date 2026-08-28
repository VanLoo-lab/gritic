import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from gritic import gritictimer, sampletools, tableschemas, timingio
from gritic.tableschemas import (
    ROUTE_PARTICLES_REPRESENTATION,
    TIMING_REPRESENTATION_COLUMN,
    UNIFORM_NO_GAIN_REPRESENTATION,
)


class ProcessSampleIntegrationTest(unittest.TestCase):
    @staticmethod
    def make_sample():
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1'],
            'Segment_Start': [0],
            'Segment_End': [1_000],
            'Major_CN': [1],
            'Minor_CN': [1],
        })
        mutation_table = pd.DataFrame({
            'Chromosome': ['1'] * 12,
            'Position': list(range(10, 22)),
            'Tumor_Ref_Count': [20] * 12,
            'Tumor_Alt_Count': [5] * 12,
        })
        return sampletools.Sample(
            mutation_table,
            copy_number_table,
            None,
            sample_id='SMOKE',
            purity=0.8,
        )

    @staticmethod
    def make_gained_sample():
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Segment_Start': [0, 1_000],
            'Segment_End': [1_000, 1_200],
            'Major_CN': [1, 2],
            'Minor_CN': [1, 1],
        })
        mutation_table = pd.DataFrame({
            'Chromosome': ['1'] * 24,
            'Position': (
                list(range(10, 22)) + list(range(1_010, 1_022))
            ),
            'Tumor_Ref_Count': [20] * 24,
            'Tumor_Alt_Count': [5] * 24,
        })
        return sampletools.Sample(
            mutation_table,
            copy_number_table,
            None,
            sample_id='GAIN',
            purity=0.8,
            merge_cn=False,
        )

    @staticmethod
    def make_wgd_sample():
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1'],
            'Segment_Start': [0],
            'Segment_End': [1_000],
            'Major_CN': [2],
            'Minor_CN': [1],
        })
        mutation_table = pd.DataFrame({
            'Chromosome': ['1'] * 12,
            'Position': list(range(10, 22)),
            'Tumor_Ref_Count': [20] * 12,
            'Tumor_Alt_Count': [5] * 12,
        })
        return sampletools.Sample(
            mutation_table,
            copy_number_table,
            None,
            sample_id='WGD',
            purity=0.8,
            merge_cn=False,
        )

    def test_clonal_diploid_sample_runs_through_every_output_boundary(self):
        sample = self.make_sample()

        with tempfile.TemporaryDirectory() as temporary_directory:
            gritictimer.process_sample(
                sample,
                temporary_directory,
                wgd_count=0,
                random_seed=20260828,
            )
            output_directory = Path(temporary_directory) / sample.sample_id

            expected_files = {
                'SMOKE_gain_timing_table.tsv',
                'SMOKE_mutation_table.tsv',
                'SMOKE_posterior_timing_table_summary_penalty_False.tsv',
                'SMOKE_posterior_timing_table_summary_penalty_True.tsv',
                'SMOKE_route_table.tsv',
                'SMOKE_subclone_table.tsv',
                'SMOKE_wgd_calling_info.json',
            }
            self.assertTrue(
                expected_files.issubset(
                    path.name
                    for path in output_directory.iterdir()
                    if path.is_file()
                )
            )

            route_table = pd.read_csv(
                output_directory / 'SMOKE_route_table.tsv',
                sep='\t',
            )
            self.assertEqual(len(route_table), 1)
            self.assertEqual(
                route_table.columns.tolist(),
                tableschemas.ROUTE_TABLE_COLUMNS,
            )
            self.assertEqual(route_table.loc[0, 'Sample_ID'], 'SMOKE')
            self.assertEqual(route_table.loc[0, 'Segment_ID'], '1-0-1000')
            self.assertEqual(route_table.loc[0, 'Probability'], 1.0)
            self.assertEqual(
                route_table.loc[0, 'Penalized_Probability'],
                1.0,
            )
            self.assertEqual(
                route_table.loc[0, TIMING_REPRESENTATION_COLUMN],
                UNIFORM_NO_GAIN_REPRESENTATION,
            )

            gain_table = pd.read_csv(
                output_directory / 'SMOKE_gain_timing_table.tsv',
                sep='\t',
            )
            self.assertTrue(gain_table.empty)
            self.assertEqual(
                gain_table.columns.tolist(),
                tableschemas.GAIN_TIMING_TABLE_COLUMNS,
            )

            mutation_table = pd.read_csv(
                output_directory / 'SMOKE_mutation_table.tsv',
                sep='\t',
            )
            self.assertEqual(len(mutation_table), 12)
            self.assertTrue(mutation_table['GRITIC_Mutation_ID'].is_unique)
            self.assertEqual(
                sorted(mutation_table['Segment_Mutation_Index']),
                list(range(12)),
            )

            subclone_table = pd.read_csv(
                output_directory / 'SMOKE_subclone_table.tsv',
                sep='\t',
            )
            self.assertTrue(subclone_table.empty)

            for apply_penalty in (False, True):
                summary = pd.read_csv(
                    output_directory
                    / (
                        'SMOKE_posterior_timing_table_summary_'
                        f'penalty_{apply_penalty}.tsv'
                    ),
                    sep='\t',
                )
                self.assertTrue(summary.empty)
                self.assertEqual(
                    summary.columns.tolist(),
                    tableschemas.SUMMARY_COLUMNS,
                )

            with (
                output_directory / 'SMOKE_wgd_calling_info.json'
            ).open(encoding='utf-8') as calling_info_file:
                calling_info = json.load(calling_info_file)
            self.assertEqual(
                calling_info,
                {
                    'WGD_Timing': None,
                    'WGD_Timing_CI_Low': None,
                    'WGD_Timing_CI_High': None,
                    'Major_CN_Mode': 1,
                    'Overlap_Proportion': None,
                    'WGD_Status': False,
                    'Best_Overlap_Timing': None,
                },
            )

            archive_path, manifest_path = timingio.get_timing_archive_paths(
                output_directory / 'SMOKE_timing_dicts',
                '1-0-1000',
            )
            self.assertFalse(archive_path.exists())
            self.assertFalse(manifest_path.exists())

    def test_existing_sample_output_is_not_overwritten(self):
        sample = self.make_sample()

        with tempfile.TemporaryDirectory() as temporary_directory:
            sample_output = Path(temporary_directory) / sample.sample_id
            sample_output.mkdir()
            (sample_output / 'existing-output').write_text(
                'keep',
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                FileExistsError,
                'must be absent or empty',
            ):
                gritictimer.process_sample(
                    sample,
                    temporary_directory,
                    wgd_count=0,
                )

    def test_gained_segment_produces_aligned_route_and_posterior_outputs(self):
        sample = self.make_gained_sample()

        with tempfile.TemporaryDirectory() as temporary_directory:
            gritictimer.process_sample(
                sample,
                temporary_directory,
                wgd_count=0,
                random_seed=20260829,
            )
            output_directory = Path(temporary_directory) / sample.sample_id

            route_table = pd.read_csv(
                output_directory / 'GAIN_route_table.tsv',
                sep='\t',
            )
            probability_sums = route_table.groupby('Segment_ID')[
                ['Probability', 'Penalized_Probability']
            ].sum()
            np.testing.assert_allclose(probability_sums, 1.0)

            gained_routes = route_table.loc[
                route_table['Segment_ID'].eq('1-1000-1200')
            ]
            self.assertEqual(len(gained_routes), 1)
            self.assertEqual(
                gained_routes.iloc[0][TIMING_REPRESENTATION_COLUMN],
                ROUTE_PARTICLES_REPRESENTATION,
            )
            route_id = gained_routes.iloc[0]['Route']

            gain_table = pd.read_csv(
                output_directory / 'GAIN_gain_timing_table.tsv',
                sep='\t',
            )
            self.assertEqual(len(gain_table), 1)
            self.assertEqual(gain_table.loc[0, 'Route'], route_id)
            self.assertEqual(gain_table.loc[0, 'Node_Phasing'], 'Major')
            self.assertTrue(
                gain_table.loc[0, [
                    'Timing',
                    'Timing_CI_Low',
                    'Timing_CI_High',
                ]].between(0, 1).all()
            )
            self.assertLessEqual(
                gain_table.loc[0, 'Timing_CI_Low'],
                gain_table.loc[0, 'Timing_CI_High'],
            )

            for apply_penalty in (False, True):
                summary = pd.read_csv(
                    output_directory
                    / (
                        'GAIN_posterior_timing_table_summary_'
                        f'penalty_{apply_penalty}.tsv'
                    ),
                    sep='\t',
                )
                self.assertEqual(len(summary), 1)
                self.assertEqual(summary.loc[0, 'Gain_Index'], 1)
                self.assertEqual(summary.loc[0, 'Proportion'], 1.0)
                self.assertTrue(
                    summary.loc[0, [
                        'Timing_Median',
                        'Timing_Low_CI',
                        'Timing_High_CI',
                    ]].between(0, 1).all()
                )

            archive_path, manifest_path = timingio.get_timing_archive_paths(
                output_directory / 'GAIN_timing_dicts',
                '1-1000-1200',
            )
            timing_hierarchy = timingio.load_timing_archive(
                archive_path,
                manifest_path,
            )
            self.assertEqual(set(timing_hierarchy), {route_id})
            route_timing = timing_hierarchy[route_id]
            node = int(gain_table.loc[0, 'Node'])
            expected_shape = (
                gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT,
            )
            self.assertEqual(route_timing['Timing'][node].shape, expected_shape)
            self.assertEqual(
                route_timing['Timing']['WGD'].shape,
                expected_shape,
            )
            self.assertTrue(np.isnan(route_timing['Timing']['WGD']).all())
            self.assertEqual(
                route_timing['Mult'].shape[0],
                gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT,
            )

    def test_forced_wgd_run_aligns_overlap_call_json_and_pooled_archive(self):
        sample = self.make_wgd_sample()

        with tempfile.TemporaryDirectory() as temporary_directory:
            gritictimer.process_sample(
                sample,
                temporary_directory,
                wgd_count=1,
                random_seed=20260830,
            )
            output_directory = Path(temporary_directory) / sample.sample_id

            with (
                output_directory / 'WGD_wgd_calling_info.json'
            ).open(encoding='utf-8') as calling_info_file:
                calling_info = json.load(calling_info_file)
            self.assertEqual(calling_info['Major_CN_Mode'], 2)
            self.assertIs(calling_info['WGD_Status'], True)
            self.assertEqual(calling_info['Overlap_Proportion'], 1.0)
            self.assertTrue(
                all(
                    0 <= calling_info[key] <= 1
                    for key in (
                        'WGD_Timing',
                        'WGD_Timing_CI_Low',
                        'WGD_Timing_CI_High',
                        'Best_Overlap_Timing',
                    )
                )
            )
            self.assertLessEqual(
                calling_info['WGD_Timing_CI_Low'],
                calling_info['WGD_Timing_CI_High'],
            )

            candidate_table = pd.read_csv(
                output_directory / 'WGD_gain_timing_table_wgd_segments.tsv',
                sep='\t',
            )
            self.assertEqual(len(candidate_table), 1)
            self.assertEqual(candidate_table.loc[0, 'Segment_ID'], '1-0-1000')
            self.assertEqual(candidate_table.loc[0, 'Node_Phasing'], 'Major')
            self.assertTrue(bool(candidate_table.loc[0, 'Intersecting']))
            self.assertAlmostEqual(
                candidate_table.loc[0, 'Overlap_Proportion'],
                calling_info['Overlap_Proportion'],
            )
            self.assertAlmostEqual(
                candidate_table.loc[0, 'Best_Overlap_Timing'],
                calling_info['Best_Overlap_Timing'],
            )

            route_table = pd.read_csv(
                output_directory / 'WGD_route_table.tsv',
                sep='\t',
            )
            self.assertTrue(route_table.empty)
            self.assertEqual(
                route_table.columns.tolist(),
                tableschemas.ROUTE_TABLE_COLUMNS,
            )

            archive_path, manifest_path = timingio.get_timing_archive_paths(
                output_directory / 'WGD_timing_dicts',
                'WGD_minor_cn_1',
            )
            timing_hierarchy = timingio.load_timing_archive(
                archive_path,
                manifest_path,
            )
            self.assertEqual(len(timing_hierarchy), 1)
            route_timing = next(iter(timing_hierarchy.values()))
            expected_shape = (
                gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT,
            )
            self.assertEqual(
                route_timing['Timing']['WGD'].shape,
                expected_shape,
            )
            self.assertTrue(np.isnan(route_timing['Timing']['WGD']).all())
            node_keys = [
                key
                for key in route_timing['Timing']
                if isinstance(key, int)
            ]
            self.assertEqual(len(node_keys), 1)
            pooled_wgd_node_timing = route_timing['Timing'][node_keys[0]]
            self.assertEqual(pooled_wgd_node_timing.shape, expected_shape)
            self.assertTrue(np.isfinite(pooled_wgd_node_timing).all())
            self.assertTrue(
                (
                    (pooled_wgd_node_timing >= 0)
                    & (pooled_wgd_node_timing <= 1)
                ).all()
            )


if __name__ == '__main__':
    unittest.main()
