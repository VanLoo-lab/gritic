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
                'SMOKE_count_group_table.tsv',
                'SMOKE_count_group_likelihood_table.tsv',
                'SMOKE_gain_timing_table.tsv',
                'SMOKE_likelihood_context_table.tsv',
                'SMOKE_mutation_table.tsv',
                'SMOKE_phase_group_table.tsv',
                'SMOKE_posterior_timing_table_summary_penalty_False.tsv',
                'SMOKE_posterior_timing_table_summary_penalty_True.tsv',
                'SMOKE_route_table.tsv',
                'SMOKE_segment_context_table.tsv',
                'SMOKE_segment_group_table.tsv',
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
            self.assertIn('Phase_Group_ID', mutation_table.columns)
            self.assertFalse(any(
                column.startswith('Prob_')
                or column.startswith('Alt_Count_Correction_')
                for column in mutation_table.columns
            ))
            self.assertTrue({
                'Tumor_Ref_Count',
                'Tumor_Alt_Count',
                'Phasing',
                'Major_CN',
                'Minor_CN',
                'Total_CN',
                'Gain_Type',
            }.isdisjoint(mutation_table.columns))

            count_group_table = pd.read_csv(
                output_directory / 'SMOKE_count_group_table.tsv',
                sep='\t',
            )
            self.assertEqual(
                count_group_table.columns.tolist(),
                tableschemas.COUNT_GROUP_TABLE_COLUMNS,
            )
            self.assertEqual(len(count_group_table), 1)
            self.assertEqual(count_group_table.loc[0, 'Sample_ID'], 'SMOKE')
            self.assertEqual(count_group_table.loc[0, 'Count_Group_ID'], 0)
            self.assertEqual(count_group_table.loc[0, 'Tumor_Ref_Count'], 20)
            self.assertEqual(count_group_table.loc[0, 'Tumor_Alt_Count'], 5)

            phase_group_table = pd.read_csv(
                output_directory / 'SMOKE_phase_group_table.tsv',
                sep='\t',
            )
            self.assertEqual(
                phase_group_table.columns.tolist(),
                tableschemas.PHASE_GROUP_TABLE_COLUMNS,
            )
            self.assertEqual(len(phase_group_table), 1)
            self.assertEqual(phase_group_table.loc[0, 'Sample_ID'], 'SMOKE')
            self.assertEqual(phase_group_table.loc[0, 'Phase_Group_ID'], 0)
            self.assertEqual(phase_group_table.loc[0, 'Count_Group_ID'], 0)
            self.assertEqual(
                phase_group_table.loc[0, 'Phasing'],
                'non_phased',
            )

            likelihood_context_table = pd.read_csv(
                output_directory / 'SMOKE_likelihood_context_table.tsv',
                sep='\t',
            )
            self.assertEqual(
                likelihood_context_table.columns.tolist(),
                tableschemas.LIKELIHOOD_CONTEXT_TABLE_COLUMNS,
            )
            self.assertEqual(len(likelihood_context_table), 1)
            self.assertEqual(
                likelihood_context_table.loc[0, [
                    'Major_CN',
                    'Minor_CN',
                    'Normal_Total_CN',
                ]].tolist(),
                [1, 1, 2],
            )

            segment_context_table = pd.read_csv(
                output_directory / 'SMOKE_segment_context_table.tsv',
                sep='\t',
            )
            self.assertEqual(
                segment_context_table.columns.tolist(),
                tableschemas.SEGMENT_CONTEXT_TABLE_COLUMNS,
            )
            self.assertEqual(len(segment_context_table), 1)
            self.assertEqual(
                segment_context_table.loc[0, 'Segment_ID'],
                '1-0-1000',
            )
            self.assertEqual(
                segment_context_table.loc[0, 'Likelihood_Context_ID'],
                0,
            )

            count_group_likelihood_table = pd.read_csv(
                output_directory
                / 'SMOKE_count_group_likelihood_table.tsv',
                sep='\t',
            )
            self.assertEqual(
                count_group_likelihood_table.columns.tolist(),
                tableschemas.COUNT_GROUP_LIKELIHOOD_BASE_COLUMNS
                + ['Prob_Mult_1'],
            )
            self.assertEqual(len(count_group_likelihood_table), 1)
            self.assertEqual(
                count_group_likelihood_table.loc[0, 'Prob_Mult_1'],
                1.0,
            )

            segment_group_table = pd.read_csv(
                output_directory / 'SMOKE_segment_group_table.tsv',
                sep='\t',
            )
            self.assertEqual(
                segment_group_table.columns.tolist(),
                tableschemas.SEGMENT_GROUP_TABLE_COLUMNS,
            )
            self.assertEqual(len(segment_group_table), 1)
            self.assertEqual(segment_group_table.loc[0, 'N_Mutations'], 12)

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

            count_group_table = pd.read_csv(
                output_directory / 'GAIN_count_group_table.tsv',
                sep='\t',
            )
            phase_group_table = pd.read_csv(
                output_directory / 'GAIN_phase_group_table.tsv',
                sep='\t',
            )
            count_group_likelihood_table = pd.read_csv(
                output_directory
                / 'GAIN_count_group_likelihood_table.tsv',
                sep='\t',
            )
            segment_group_table = pd.read_csv(
                output_directory / 'GAIN_segment_group_table.tsv',
                sep='\t',
            )
            mutation_table = pd.read_csv(
                output_directory / 'GAIN_mutation_table.tsv',
                sep='\t',
            )
            self.assertEqual(len(count_group_table), 1)
            self.assertEqual(len(phase_group_table), 1)
            self.assertEqual(len(count_group_likelihood_table), 2)
            self.assertEqual(segment_group_table['N_Mutations'].sum(), 24)
            self.assertEqual(len(mutation_table), 24)
            self.assertEqual(
                segment_group_table.groupby('Segment_ID')[
                    'N_Mutations'
                ].sum().sort_index().tolist(),
                [12, 12],
            )

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
            self.assertEqual(route_timing['Timing'].shape, (expected_shape[0], 1))
            np.testing.assert_array_equal(route_timing['Timing_Node_ID'], [node])
            self.assertEqual(
                route_timing['WGD_Timing'].shape,
                expected_shape,
            )
            self.assertTrue(np.isnan(route_timing['WGD_Timing']).all())
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
                route_timing['WGD_Timing'].shape,
                expected_shape,
            )
            self.assertTrue(np.isnan(route_timing['WGD_Timing']).all())
            self.assertEqual(route_timing['Timing_Node_ID'].shape, (1,))
            pooled_wgd_node_timing = route_timing['Timing'][:, 0]
            self.assertEqual(pooled_wgd_node_timing.shape, expected_shape)
            self.assertTrue(np.isfinite(pooled_wgd_node_timing).all())
            self.assertTrue(
                (
                    (pooled_wgd_node_timing >= 0)
                    & (pooled_wgd_node_timing <= 1)
                ).all()
            )
            np.testing.assert_array_equal(route_timing['Archive_Kind'], [1])
            np.testing.assert_array_equal(route_timing['Target_Major_CN'], [2])
            np.testing.assert_array_equal(route_timing['Target_Minor_CN'], [1])
            np.testing.assert_array_equal(route_timing['Target_WGD_Status'], [1])
            np.testing.assert_array_equal(route_timing['Model_Major_CN'], [2])
            np.testing.assert_array_equal(route_timing['Model_Minor_CN'], [1])
            np.testing.assert_array_equal(route_timing['Model_WGD_Status'], [0])


if __name__ == '__main__':
    unittest.main()
