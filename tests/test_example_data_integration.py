import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from gritic import dataloader, gritictimer, sampletools, timingio


class ExampleDataIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        examples = Path(__file__).resolve().parents[1] / 'examples'
        cls.copy_number_table, cls.mutation_table = (
            dataloader.load_input_tables(
                examples / 'cn_table_example.tsv',
                examples / 'snv_table_example.tsv',
            )
        )
        cls.input_subclone_table = pd.read_csv(
            examples / 'subclone_table_example.tsv',
            sep='\t',
        )
        cls.sample = sampletools.Sample._from_validated_input_tables(
            cls.mutation_table,
            cls.copy_number_table,
            cls.input_subclone_table,
            sample_id='EXAMPLE_INTEGRATION',
            purity=0.5,
        )

    def test_documented_example_builds_a_consistent_sample(self):
        sample = self.sample

        self.assertGreater(len(self.mutation_table), 5_000)
        self.assertGreater(len(sample.copy_number_table), 50)
        self.assertGreater(len(sample.segments), 50)
        self.assertEqual(
            sum(segment.n_mutations for segment in sample.segments),
            len(sample.mutation_table),
        )
        self.assertTrue(sample.copy_number_table['Segment_ID'].is_unique)
        self.assertTrue(sample.mutation_table['GRITIC_Mutation_ID'].is_unique)
        self.assertFalse(sample.mutation_table['Segment_ID'].isna().any())

        mutation_rows = sample.mutation_table
        self.assertTrue(
            (mutation_rows['Position'] >= mutation_rows['Segment_Start']).all()
        )
        self.assertTrue(
            (mutation_rows['Position'] < mutation_rows['Segment_End']).all()
        )
        self.assertTrue(
            mutation_rows['Total_CN'].eq(
                mutation_rows['Major_CN'] + mutation_rows['Minor_CN']
            ).all()
        )

    def test_segment_mutation_indexes_are_canonical_and_consecutive(self):
        mutation_table = self.sample.get_mutation_table()

        for _, segment_rows in mutation_table.groupby(
            'Segment_ID',
            observed=True,
            sort=False,
        ):
            canonical_rows = segment_rows.sort_values(
                'GRITIC_Mutation_ID',
                kind='mergesort',
            )
            self.assertEqual(
                canonical_rows['Segment_Mutation_Index'].tolist(),
                list(range(len(canonical_rows))),
            )

    def test_every_segment_has_normalized_finite_multiplicity_models(self):
        for segment in self.sample.segments:
            probabilities = segment.get_multiplicity_probabilities_array(
                'all'
            )
            corrections = segment.get_reads_correction_array('all')

            self.assertEqual(probabilities.shape[0], segment.n_mutations)
            self.assertEqual(
                probabilities.shape[1],
                segment.major_cn + segment.n_subclones,
            )
            np.testing.assert_allclose(probabilities.sum(axis=1), 1.0)
            self.assertTrue(np.isfinite(probabilities).all())
            self.assertTrue((probabilities >= 0).all())
            self.assertTrue(np.isfinite(corrections).all())
            self.assertTrue(((corrections >= 0) & (corrections <= 1)).all())

    def test_example_subclone_and_wgd_inputs_survive_processing(self):
        subclones = self.sample.get_subclone_table()

        self.assertEqual(
            subclones.columns.tolist(),
            list(dataloader.SUBCLONE_OUTPUT_COLUMNS),
        )
        np.testing.assert_array_equal(
            subclones['N_SNVs'],
            np.round(
                len(self.sample.mutation_table)
                * subclones['Subclone_Fraction']
            ).astype(int),
        )
        self.assertEqual(sampletools.get_major_cn_mode(self.sample), 2)

    def test_downsampled_real_data_runs_forced_non_wgd_smoke(self):
        selected_segments = (
            ('1', 21_089_441, 247_489_947),
            ('3', 0, 116_026_931),
            ('2', 21_726_957, 186_101_780),
        )
        copy_number_parts = []
        mutation_parts = []
        for chromosome, segment_start, segment_end in selected_segments:
            copy_number_parts.append(self.copy_number_table.loc[
                self.copy_number_table['Chromosome'].eq(chromosome)
                & self.copy_number_table['Segment_Start'].eq(segment_start)
                & self.copy_number_table['Segment_End'].eq(segment_end)
            ])
            mutation_parts.append(self.mutation_table.loc[
                self.mutation_table['Chromosome'].eq(chromosome)
                & self.mutation_table['Position'].ge(segment_start)
                & self.mutation_table['Position'].lt(segment_end)
            ].head(12))

        sample = sampletools.Sample._from_validated_input_tables(
            pd.concat(mutation_parts, ignore_index=True),
            pd.concat(copy_number_parts, ignore_index=True),
            self.input_subclone_table,
            sample_id='EXAMPLE_DOWNSAMPLED',
            purity=0.5,
            merge_cn=False,
        )
        self.assertEqual(sampletools.get_major_cn_mode(sample), 2)

        with tempfile.TemporaryDirectory() as directory:
            with self.assertWarnsRegex(
                UserWarning,
                'WGD count 0 but major CN mode is 2',
            ):
                gritictimer.process_sample(
                    sample,
                    directory,
                    wgd_count=0,
                    random_seed=20260901,
                )
            output_directory = Path(directory) / sample.sample_id

            subclone_table = pd.read_csv(
                output_directory / 'EXAMPLE_DOWNSAMPLED_subclone_table.tsv',
                sep='\t',
            )
            self.assertEqual(len(subclone_table), 1)
            self.assertEqual(subclone_table.loc[0, 'N_SNVs'], 7)

            route_table = pd.read_csv(
                output_directory / 'EXAMPLE_DOWNSAMPLED_route_table.tsv',
                sep='\t',
            )
            self.assertEqual(
                set(zip(route_table['Major_CN'], route_table['Minor_CN'])),
                {(2, 0), (2, 1), (3, 1)},
            )
            np.testing.assert_allclose(
                route_table.groupby('Segment_ID')['Probability'].sum(),
                1.0,
            )

            archive_path, manifest_path = timingio.get_timing_archive_paths(
                output_directory / 'EXAMPLE_DOWNSAMPLED_timing_dicts',
                '2-21726957-186101780',
            )
            timing_hierarchy = timingio.load_timing_archive(
                archive_path,
                manifest_path,
            )
            self.assertTrue(timing_hierarchy)
            for route_timing in timing_hierarchy.values():
                self.assertEqual(route_timing['Mult'].shape[1], 8)
                self.assertTrue(np.isfinite(route_timing['Mult']).all())
                self.assertTrue(
                    (
                        (route_timing['Mult'][:, -1] >= 0)
                        & (route_timing['Mult'][:, -1] <= 1)
                    ).all()
                )


if __name__ == '__main__':
    unittest.main()
