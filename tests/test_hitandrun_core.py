import unittest

import numpy as np
from scipy.linalg import null_space

from gritic import hitandrun


class DirectionGeometryTest(unittest.TestCase):
    def test_random_direction_has_requested_dimension_and_unit_norm(self):
        for dimension in (1, 2, 7):
            with self.subTest(dimension=dimension):
                direction = hitandrun.get_random_direction(dimension)
                self.assertEqual(direction.shape, (dimension,))
                self.assertTrue(np.isfinite(direction).all())
                self.assertAlmostEqual(np.linalg.norm(direction), 1.0, 13)

    def test_positive_direction_limit_stops_at_first_zero_coordinate(self):
        direction = np.array([-0.5, -0.25, 0.75])
        current_position = np.array([0.2, 0.8, 0.0])

        limit = hitandrun.get_direction_limit(direction, current_position)

        # Coordinate zero reaches zero after 0.4 units; coordinate one would
        # permit 3.2 units.
        self.assertAlmostEqual(limit, 0.4)

    def test_direction_range_finds_both_simplex_boundaries(self):
        direction = np.array([-1.0, 1.0]) / np.sqrt(2.0)
        current_position = np.array([0.25, 0.75])

        lower, upper = hitandrun.get_direction_range(
            direction,
            current_position,
        )

        self.assertAlmostEqual(lower, -0.75 * np.sqrt(2.0))
        self.assertAlmostEqual(upper, 0.25 * np.sqrt(2.0))
        np.testing.assert_allclose(
            current_position + lower * direction,
            [1.0, 0.0],
            atol=1e-14,
        )
        np.testing.assert_allclose(
            current_position + upper * direction,
            [0.0, 1.0],
            atol=1e-14,
        )

    def test_new_position_stays_on_simplex_and_does_not_mutate_input(self):
        constraint = np.ones((1, 4))
        null_basis = null_space(constraint)
        start = np.array([0.1, 0.2, 0.3, 0.4])
        original = start.copy()

        for _ in range(12):
            sampled = hitandrun.get_new_position(
                start,
                null_basis.shape[1],
                null_basis,
            )
            self.assertTrue((sampled >= 0.0).all())
            self.assertTrue((sampled <= 1.0).all())
            self.assertAlmostEqual(sampled.sum(), 1.0, 12)

        np.testing.assert_array_equal(start, original)


class HitAndRunChainTest(unittest.TestCase):
    def test_zero_dimensional_polytope_repeats_the_only_state(self):
        state = np.array([0.2, 0.3, 0.5])

        samples, process_limit_exceeded = hitandrun._hit_and_run(
            np.empty((3, 0)),
            state,
            n_samples=6,
            burn_in=100,
            skips=20,
        )

        self.assertFalse(process_limit_exceeded)
        np.testing.assert_array_equal(samples, np.tile(state, (6, 1)))

    def test_one_dimensional_polytope_is_sampled_at_evenly_spaced_points(self):
        state = np.array([0.25, 0.75])
        null_basis = np.array([[-1.0], [1.0]]) / np.sqrt(2.0)

        samples, process_limit_exceeded = hitandrun._hit_and_run(
            null_basis,
            state,
            n_samples=5,
        )

        self.assertFalse(process_limit_exceeded)
        ordered = samples[np.argsort(samples[:, 0])]
        np.testing.assert_allclose(
            ordered[:, 0],
            np.linspace(0.0, 1.0, 5),
            atol=1e-14,
        )
        np.testing.assert_allclose(ordered.sum(axis=1), 1.0, atol=1e-14)
        np.testing.assert_allclose(ordered[:, 1], 1.0 - ordered[:, 0])

    def test_run_chain_honours_burn_in_skip_and_sample_count(self):
        state = np.array([0.2, 0.3, 0.5])
        null_basis = null_space(np.ones((1, 3)))
        burn_in = 4
        skips = 3
        n_samples = 7
        dense_sample_count = burn_in + skips * n_samples

        hitandrun.seed_random(44_021)
        dense_samples = hitandrun.run_chain(
            state,
            null_basis.shape[1],
            null_basis,
            state,
            burn_in=0,
            skips=1,
            n_samples=dense_sample_count,
        )
        hitandrun.seed_random(44_021)
        thinned_samples = hitandrun.run_chain(
            state,
            null_basis.shape[1],
            null_basis,
            state,
            burn_in=burn_in,
            skips=skips,
            n_samples=n_samples,
        )

        np.testing.assert_array_equal(
            thinned_samples,
            dense_samples[burn_in::skips],
        )
        self.assertEqual(thinned_samples.shape, (n_samples, 3))

    def test_multidimensional_sampler_preserves_linear_constraints(self):
        constraint = np.array([
            [1.0, 1.0, 1.0, 1.0],
            [1.0, -1.0, 0.0, 0.0],
        ])
        null_basis = null_space(constraint)
        state = np.array([0.2, 0.2, 0.25, 0.35])
        expected_constraint_values = constraint @ state

        samples, process_limit_exceeded = hitandrun._hit_and_run(
            null_basis,
            state,
            n_samples=20,
            burn_in=5,
            skips=2,
        )

        self.assertFalse(process_limit_exceeded)
        self.assertEqual(samples.shape, (20, 4))
        self.assertTrue((samples >= 0.0).all())
        self.assertTrue((samples <= 1.0).all())
        expected = np.tile(expected_constraint_values, (20, 1))
        np.testing.assert_allclose(samples @ constraint.T, expected, atol=1e-12)
        self.assertTrue(np.any(np.linalg.norm(samples - state, axis=1) > 1e-3))

    def test_public_sampler_returns_same_fixed_polytope_contract(self):
        state = np.array([0.4, 0.6])

        samples = hitandrun.hit_and_run(
            np.empty((2, 0)),
            state,
            n_samples=4,
        )

        np.testing.assert_array_equal(samples, np.tile(state, (4, 1)))


if __name__ == '__main__':
    unittest.main()
