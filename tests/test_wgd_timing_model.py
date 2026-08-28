import unittest

import numpy as np
import pandas as pd

from gritic import gritictimer, sampletools


class WgdTimingModelTest(unittest.TestCase):
    @staticmethod
    def make_segment(minor_cn, phasing):
        row_count = len(phasing)
        mutation_table = pd.DataFrame({
            'Segment_ID': ['segment'] * row_count,
            'Chromosome': ['1'] * row_count,
            'Segment_Start': [0] * row_count,
            'Segment_End': [100] * row_count,
            'Major_CN': [2] * row_count,
            'Minor_CN': [minor_cn] * row_count,
            'Tumor_Ref_Count': np.arange(20, 20 + row_count),
            'Tumor_Alt_Count': np.arange(4, 4 + row_count),
            'Phasing': phasing,
        })
        return sampletools.Segment(
            mutation_table,
            None,
            purity=0.8,
            sex='XX',
        )

    def test_balanced_segment_uses_every_phase_group_once(self):
        segment = self.make_segment(2, [np.nan, 'major', 'minor'])
        source_store = segment.multiplicity_probabilities
        source_phasing = segment.mutation_table['Phasing'].copy()

        pseudo_minor_cn, wgd_store = gritictimer._get_wgd_timing_model(
            segment
        )

        self.assertEqual(pseudo_minor_cn, 0)
        self.assertIsNot(wgd_store, source_store)
        self.assertEqual(wgd_store.major_cn, 2)
        self.assertEqual(wgd_store.minor_cn, 0)
        self.assertEqual(wgd_store.n_mutations, 3)
        self.assertTrue(wgd_store.use_non_phased)
        self.assertFalse(wgd_store.use_major)
        self.assertFalse(wgd_store.use_minor)
        np.testing.assert_allclose(
            wgd_store.non_phased_array,
            source_store.combined_array,
        )
        np.testing.assert_allclose(
            wgd_store.combined_array,
            source_store.combined_array,
        )
        np.testing.assert_allclose(
            wgd_store.reads_correction_non_phased_array,
            source_store.reads_correction_combined_array,
        )

        self.assertEqual(source_store.major_cn, 2)
        self.assertEqual(source_store.minor_cn, 2)
        self.assertTrue(source_store.use_non_phased)
        self.assertTrue(source_store.use_major)
        self.assertTrue(source_store.use_minor)
        pd.testing.assert_series_equal(
            segment.mutation_table['Phasing'],
            source_phasing,
        )

    def test_minor_phased_row_is_evaluated_without_broadcasting(self):
        segment = self.make_segment(2, [np.nan, 'major', 'minor'])
        _, wgd_store = gritictimer._get_wgd_timing_model(segment)
        full_states = np.array([[0.6, 0.4, 0.6, 0.4]])

        likelihood = wgd_store.evaluate_likelihood_array(full_states)

        corrected_state = (
            full_states[:, :2]
            * wgd_store.reads_correction_non_phased_array
        )
        corrected_state /= corrected_state.sum(axis=1, keepdims=True)
        per_mutation_likelihood = (
            wgd_store.non_phased_array @ corrected_state[0]
        )
        expected = np.log(per_mutation_likelihood + 2.2e-300).sum()
        np.testing.assert_allclose(likelihood, [expected])
        self.assertTrue(np.isfinite(likelihood).all())

    def test_non_balanced_segment_reuses_original_model(self):
        segment = self.make_segment(1, [np.nan, 'major', 'minor'])
        source_store = segment.multiplicity_probabilities

        pseudo_minor_cn, wgd_store = gritictimer._get_wgd_timing_model(
            segment
        )

        self.assertEqual(pseudo_minor_cn, 1)
        self.assertIs(wgd_store, source_store)
        self.assertEqual(source_store.major_cn, 2)
        self.assertEqual(source_store.minor_cn, 1)
        self.assertTrue(source_store.use_non_phased)
        self.assertTrue(source_store.use_major)
        self.assertTrue(source_store.use_minor)
        likelihood = wgd_store.evaluate_likelihood_array(
            np.array([[0.6, 0.4, 0.6, 0.4, 1.0]])
        )
        self.assertTrue(np.isfinite(likelihood).all())


if __name__ == '__main__':
    unittest.main()
