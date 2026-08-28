import tempfile
import unittest
from collections.abc import Mapping
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.linalg import null_space

from gritic import gritictimer, hitandrun, sampletools, timingio


def make_subclonal_gain_sample():
    copy_number_table = pd.DataFrame({
        'Chromosome': ['1', '1'],
        'Segment_Start': [0, 1_000],
        'Segment_End': [1_000, 1_200],
        'Major_CN': [1, 2],
        'Minor_CN': [1, 1],
    })
    mutation_table = pd.DataFrame({
        'Chromosome': ['1'] * 24,
        'Position': list(range(10, 22)) + list(range(1_010, 1_022)),
        'Tumor_Ref_Count': [20] * 24,
        'Tumor_Alt_Count': [5] * 24,
    })
    subclone_table = pd.DataFrame({
        'Cluster': ['subclone-a'],
        'Subclone_CCF': [0.4],
        'Subclone_Fraction': [0.2],
    })
    return sampletools.Sample(
        mutation_table,
        copy_number_table,
        subclone_table,
        sample_id='SEEDED_SUBCLONE',
        purity=0.8,
        merge_cn=False,
    )


class RandomSeedValidationTest(unittest.TestCase):
    def test_validator_accepts_full_numpy_seed_range(self):
        for seed in (0, 1, np.int64(17), 2**32 - 1, None):
            with self.subTest(seed=seed):
                expected = None if seed is None else int(seed)
                self.assertEqual(
                    gritictimer._validate_random_seed(seed),
                    expected,
                )

    def test_process_sample_rejects_invalid_seed_before_creating_output(self):
        invalid_seeds = (
            True,
            np.bool_(False),
            -1,
            2**32,
            1.5,
            '17',
        )
        sample = SimpleNamespace(sample_id='INVALID_SEED')
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / sample.sample_id
            for seed in invalid_seeds:
                with self.subTest(seed=seed):
                    with self.assertRaisesRegex(
                        ValueError,
                        r'random_seed.*between 0 and 2\*\*32 - 1',
                    ):
                        gritictimer.process_sample(
                            sample,
                            temporary_directory,
                            random_seed=seed,
                        )
                    self.assertFalse(output_path.exists())

    def test_numba_hit_and_run_stream_can_be_reseeded_exactly(self):
        state = np.array([0.2, 0.3, 0.5])
        null_basis = null_space(np.ones((1, state.size)))

        hitandrun.seed_random(90210)
        first = hitandrun.hit_and_run(
            null_basis,
            state,
            n_samples=12,
            burn_in=4,
            skips=2,
        )
        hitandrun.seed_random(90210)
        repeated = hitandrun.hit_and_run(
            null_basis,
            state,
            n_samples=12,
            burn_in=4,
            skips=2,
        )
        hitandrun.seed_random(90211)
        different = hitandrun.hit_and_run(
            null_basis,
            state,
            n_samples=12,
            burn_in=4,
            skips=2,
        )

        np.testing.assert_array_equal(first, repeated)
        self.assertFalse(np.array_equal(first, different))


class SeededSubcloneProcessSampleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directories = [
            tempfile.TemporaryDirectory(),
            tempfile.TemporaryDirectory(),
        ]
        for temporary_directory in cls.temporary_directories:
            cls.addClassCleanup(temporary_directory.cleanup)

        sample = make_subclonal_gain_sample()
        cls.output_directories = []
        for temporary_directory in cls.temporary_directories:
            gritictimer.process_sample(
                sample,
                temporary_directory.name,
                wgd_count=0,
                random_seed=731_921,
            )
            cls.output_directories.append(
                Path(temporary_directory.name) / sample.sample_id
            )

    def assert_nested_exactly_equal(self, first, second):
        if isinstance(first, Mapping):
            self.assertIsInstance(second, Mapping)
            self.assertEqual(list(first), list(second))
            for key in first:
                self.assert_nested_exactly_equal(first[key], second[key])
            return
        if isinstance(first, np.ndarray):
            np.testing.assert_array_equal(first, second)
            return
        np.testing.assert_equal(first, second)

    def read_table(self, output_directory, suffix):
        return pd.read_csv(
            output_directory / f'SEEDED_SUBCLONE_{suffix}.tsv',
            sep='\t',
        )

    def load_gain_archive(self, output_directory):
        archive_path, manifest_path = timingio.get_timing_archive_paths(
            output_directory / 'SEEDED_SUBCLONE_timing_dicts',
            '1-1000-1200',
        )
        return timingio.load_timing_archive(archive_path, manifest_path)

    def test_subclonal_run_writes_nonempty_inferential_outputs(self):
        output_directory = self.output_directories[0]
        subclones = self.read_table(output_directory, 'subclone_table')
        routes = self.read_table(output_directory, 'route_table')
        gains = self.read_table(output_directory, 'gain_timing_table')
        posterior = self.read_table(
            output_directory,
            'posterior_timing_table_summary_penalty_False',
        )

        self.assertEqual(subclones['Cluster'].tolist(), ['subclone-a'])
        self.assertEqual(subclones['N_SNVs'].tolist(), [5])
        gained_routes = routes.loc[
            routes['Segment_ID'].eq('1-1000-1200')
        ]
        self.assertFalse(gained_routes.empty)
        self.assertFalse(gains.empty)
        self.assertFalse(posterior.empty)
        np.testing.assert_allclose(gained_routes['Probability'].sum(), 1.0)

        timing_hierarchy = self.load_gain_archive(output_directory)
        self.assertTrue(timing_hierarchy)
        for route_samples in timing_hierarchy.values():
            # 2 * major CN + minor CN + one subclone coordinate.
            self.assertEqual(
                route_samples['Mult'].shape,
                (gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT, 6),
            )
            subclone_shares = route_samples['Mult'][:, -1]
            self.assertTrue(np.isfinite(subclone_shares).all())
            self.assertTrue(((subclone_shares > 0) & (subclone_shares < 1)).all())

    def test_same_seed_reproduces_tables_and_full_timing_archive(self):
        first_output, second_output = self.output_directories
        table_suffixes = (
            'route_table',
            'gain_timing_table',
            'posterior_timing_table_summary_penalty_False',
            'posterior_timing_table_summary_penalty_True',
        )
        for suffix in table_suffixes:
            with self.subTest(table=suffix):
                first = self.read_table(first_output, suffix)
                second = self.read_table(second_output, suffix)
                # Wall-clock runtime is intentionally not reproducible.
                first = first.drop(columns=['Time'], errors='ignore')
                second = second.drop(columns=['Time'], errors='ignore')
                pd.testing.assert_frame_equal(
                    first,
                    second,
                    check_exact=True,
                )

        self.assert_nested_exactly_equal(
            self.load_gain_archive(first_output),
            self.load_gain_archive(second_output),
        )


if __name__ == '__main__':
    unittest.main()
