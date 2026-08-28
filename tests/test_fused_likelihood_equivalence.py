import unittest

import numpy as np

from gritic import sampletools


def reference_log_likelihood(mult_array, mult_states):
    likelihoods = np.zeros(mult_states.shape[0])
    for state_index in range(mult_states.shape[0]):
        mult_state = np.clip(mult_states[state_index], 0.0, 1.0)
        mutation_likelihoods = np.multiply(mult_state, mult_array)
        mutation_likelihoods = np.sum(mutation_likelihoods, axis=1)
        likelihoods[state_index] = np.sum(
            np.log(mutation_likelihoods + 2.2e-300)
        )
    return likelihoods


def reference_array_likelihood(
    full_states,
    probability_array,
    reads_correction_array,
    tolerance=1e-8,
):
    corrected_states = np.multiply(full_states, reads_correction_array)
    corrected_states = corrected_states / np.sum(
        corrected_states,
        axis=1,
        keepdims=True,
    )
    probability_sums = np.sum(probability_array, axis=1)
    valid_rows = (
        (probability_sums > 1.0 - tolerance)
        & (probability_sums < 1.0 + tolerance)
    )
    probability_array = probability_array[valid_rows, :]
    if probability_array.shape[0] == 0:
        return np.full(full_states.shape[0], np.nan)
    return reference_log_likelihood(probability_array, corrected_states)


class FusedLikelihoodEquivalenceTest(unittest.TestCase):
    def test_fused_kernel_matches_allocation_based_reference(self):
        rng = np.random.default_rng(2197)
        probability_array = rng.dirichlet(np.ones(7), size=83)
        full_states = rng.normal(0.4, 0.7, size=(127, 7))
        full_states[0] = [-np.inf, -0.0, 0.0, 0.5, 1.0, 1.5, np.inf]

        expected = reference_log_likelihood(
            probability_array,
            full_states,
        )
        actual = sampletools.log_likelihood_numba_parallel(
            probability_array,
            full_states,
        )

        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)

    def test_read_correction_and_probability_filter_match_reference(self):
        probability_array = np.array([
            [0.65, 0.20, 0.10, 0.05],
            [0.10, 0.20, 0.30, 0.40 + 5e-9],
            [0.30, 0.20, 0.10, 0.39],
        ])
        full_states = np.array([
            [0.30, 0.25, 0.20, 0.25],
            [1.20, -0.10, 0.40, 0.25],
            [0.00, 0.25, 0.50, 0.25],
        ])
        reads_correction_array = np.array([0.25, 0.75, 0.50, 1.00])

        expected = reference_array_likelihood(
            full_states,
            probability_array,
            reads_correction_array,
        )
        actual = sampletools.evaluate_likelihood_array_numba(
            full_states,
            probability_array,
            reads_correction_array,
        )

        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)

    def test_mult_probability_store_preserves_phase_specific_inputs(self):
        non_phased_array = np.array([
            [0.70, 0.15, 0.10, 0.05],
            [0.10, 0.20, 0.30, 0.40],
        ])
        major_array = np.array([
            [0.50, 0.20, 0.10, 0.20],
            [0.15, 0.35, 0.25, 0.25],
        ])
        minor_array = np.array([
            [0.60, 0.25, 0.15],
            [0.20, 0.30, 0.50],
        ])
        reads_correction_store = {
            'Non_Phased': np.array([0.30, 0.50, 0.70, 0.90]),
            'Major': np.array([0.20, 0.45, 0.80, 0.65]),
            'Minor': np.array([0.35, 0.75, 0.60]),
            'All': np.array([0.30, 0.50, 0.70, 0.90]),
        }
        store = sampletools.MultProbabilityStore(
            {
                'Non_Phased': non_phased_array,
                'Major': major_array,
                'Minor': minor_array,
                'All': non_phased_array,
            },
            reads_correction_store,
            major_cn=3,
            minor_cn=2,
            n_subclones=1,
        )
        full_states = np.array([
            [0.15, 0.25, 0.60, 0.40, 0.35, 0.25, 0.80, 0.20, 0.45],
            [1.10, -0.20, 0.25, 0.55, 0.20, 0.25, 0.30, 0.70, 0.35],
        ])

        expected = reference_array_likelihood(
            full_states[:, [0, 1, 2, 8]],
            non_phased_array,
            reads_correction_store['Non_Phased'],
        )
        expected += reference_array_likelihood(
            full_states[:, [3, 4, 5, 8]],
            major_array,
            reads_correction_store['Major'],
        )
        expected += reference_array_likelihood(
            full_states[:, [6, 7, 8]],
            minor_array,
            reads_correction_store['Minor'],
        )

        actual = store.evaluate_likelihood_array(full_states)

        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)

    def test_all_filtered_probability_rows_still_return_nan(self):
        full_states = np.array([[0.4, 0.6], [0.8, 0.2]])
        probability_array = np.array([[0.2, 0.3], [0.1, 0.2]])

        actual = sampletools.evaluate_likelihood_array_numba(
            full_states,
            probability_array,
            np.ones(2),
        )

        self.assertTrue(np.isnan(actual).all())


if __name__ == '__main__':
    unittest.main()
