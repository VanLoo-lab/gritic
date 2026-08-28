import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from gritic import cli, dataloader, sampletools


class PhasingLabelTest(unittest.TestCase):
    @staticmethod
    def copy_number_table():
        return pd.DataFrame({
            'Chromosome': ['1'],
            'Segment_Start': [0],
            'Segment_End': [100],
            'Major_CN': [2],
            'Minor_CN': [1],
            'Segment_ID': ['segment'],
        })

    @staticmethod
    def mutation_table(phasing):
        row_count = len(phasing)
        return pd.DataFrame({
            'Chromosome': ['1'] * row_count,
            'Tumor_Ref_Count': [12] * row_count,
            'Tumor_Alt_Count': [5] * row_count,
            'Mutation_ID': [f'mutation-{index}' for index in range(row_count)],
            'Segment_ID': ['segment'] * row_count,
            'Phasing': phasing,
        })

    @staticmethod
    def required_cli_arguments():
        return [
            '--mutation-table',
            'mutations.tsv',
            '--copy-number-table',
            'copy_number.tsv',
            '--purity',
            '0.8',
            '--sample-id',
            'TEST',
            '--output',
            'output',
        ]

    def test_canonicalizes_case_and_whitespace_and_preserves_missing(self):
        table = self.mutation_table([
            'Major ',
            'MAJOR',
            ' major ',
            'Minor ',
            'MINOR',
            ' minor ',
            None,
            np.nan,
        ])

        result = dataloader.validate_or_drop_phasing_labels(table)

        self.assertEqual(
            result['Phasing'].iloc[:6].tolist(),
            ['major', 'major', 'major', 'minor', 'minor', 'minor'],
        )
        self.assertTrue(result['Phasing'].iloc[6:].isna().all())
        self.assertEqual(table['Phasing'].iloc[0], 'Major ')

    def test_absent_phasing_column_remains_valid_and_unphased(self):
        mutation_table = self.mutation_table([None]).drop(columns='Phasing')

        validated = dataloader.validate_or_drop_phasing_labels(mutation_table)
        sample = sampletools.Sample(
            mutation_table,
            self.copy_number_table(),
            None,
            'TEST',
            0.8,
            sex='XX',
        )

        self.assertNotIn('Phasing', validated.columns)
        self.assertTrue(sample.mutation_table['Phasing'].isna().all())

    def test_rejects_every_nonmissing_label_outside_major_and_minor(self):
        for invalid_value in ('unknown', 'major allele', '', '   ', 1, True):
            with self.subTest(invalid_value=invalid_value):
                with self.assertRaisesRegex(
                    ValueError,
                    'Phasing values must be major, minor, or missing',
                ):
                    dataloader.validate_or_drop_phasing_labels(
                        self.mutation_table([invalid_value])
                    )

    def test_drop_policy_warns_once_and_drops_affected_mutations(self):
        table = self.mutation_table([
            'Major ',
            'unknown',
            None,
            ' MINOR ',
            'other',
        ])

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = dataloader.validate_or_drop_phasing_labels(
                table,
                drop_unrecognized_phasing=True,
            )

        self.assertEqual(
            result['Mutation_ID'].tolist(),
            ['mutation-0', 'mutation-2', 'mutation-3'],
        )
        self.assertEqual(result['Phasing'].iloc[[0, 2]].tolist(), [
            'major',
            'minor',
        ])
        self.assertTrue(pd.isna(result['Phasing'].iloc[1]))
        user_warnings = [
            str(item.message)
            for item in caught
            if issubclass(item.category, UserWarning)
        ]
        self.assertEqual(len(user_warnings), 1)
        self.assertIn('Dropping 2 mutations', user_warnings[0])
        self.assertIn("'unknown'", user_warnings[0])
        self.assertIn("'other'", user_warnings[0])

    def test_file_loader_canonicalizes_phasing(self):
        mutation_table = self.mutation_table([' MAJOR', 'minor ', None])
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            copy_number_path = directory / 'copy-number.tsv'
            mutation_path = directory / 'mutations.tsv'
            self.copy_number_table().to_csv(
                copy_number_path,
                sep='\t',
                index=False,
            )
            mutation_table.to_csv(mutation_path, sep='\t', index=False)

            _, loaded_mutations = dataloader.load_input_tables(
                copy_number_path,
                mutation_path,
                sex='XX',
            )

        self.assertEqual(
            loaded_mutations['Phasing'].iloc[:2].tolist(),
            ['major', 'minor'],
        )
        self.assertTrue(pd.isna(loaded_mutations['Phasing'].iloc[2]))

    def test_direct_sample_uses_canonical_labels_in_phased_likelihoods(self):
        sample = sampletools.Sample(
            self.mutation_table([' Major ', 'MINOR', None]),
            self.copy_number_table(),
            None,
            'TEST',
            0.8,
            sex='XX',
        )

        self.assertEqual(
            sample.mutation_table['Phasing'].iloc[:2].tolist(),
            ['major', 'minor'],
        )
        self.assertTrue(pd.isna(sample.mutation_table['Phasing'].iloc[2]))
        probability_store = sample.segments[0].multiplicity_probabilities
        self.assertTrue(probability_store.use_major)
        self.assertTrue(probability_store.use_minor)
        self.assertTrue(probability_store.use_non_phased)
        self.assertEqual(probability_store.major_array.shape[0], 1)
        self.assertEqual(probability_store.minor_array.shape[0], 1)
        self.assertEqual(probability_store.non_phased_array.shape[0], 1)

    def test_direct_sample_drop_policy_removes_invalid_rows(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            sample = sampletools.Sample(
                self.mutation_table(['major', 'unknown', 'minor']),
                self.copy_number_table(),
                None,
                'TEST',
                0.8,
                sex='XX',
                drop_unrecognized_phasing=True,
            )

        self.assertEqual(
            sample.mutation_table['Mutation_ID'].tolist(),
            ['mutation-0', 'mutation-2'],
        )
        self.assertTrue(any(
            'Dropping 1 mutation' in str(item.message)
            for item in caught
            if issubclass(item.category, UserWarning)
        ))

    def test_all_invalid_rows_warn_then_fail_when_building_sample(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            with self.assertRaisesRegex(
                ValueError,
                'No mutations remain after dropping mutations with '
                'unrecognized Phasing values',
            ):
                sampletools.Sample(
                    self.mutation_table(['unknown', 'other']),
                    self.copy_number_table(),
                    None,
                    'TEST',
                    0.8,
                    sex='XX',
                    drop_unrecognized_phasing=True,
                )

        self.assertTrue(any(
            'Dropping 2 mutations' in str(item.message)
            for item in caught
            if issubclass(item.category, UserWarning)
        ))

    def test_cli_flag_defaults_off_and_can_be_enabled(self):
        parser = cli.build_parser()
        default_args = parser.parse_args(self.required_cli_arguments())
        self.assertFalse(default_args.drop_unrecognized_phasing)

        enabled_args = parser.parse_args([
            *self.required_cli_arguments(),
            '--drop-unrecognized-phasing',
        ])
        self.assertTrue(enabled_args.drop_unrecognized_phasing)

if __name__ == '__main__':
    unittest.main()
