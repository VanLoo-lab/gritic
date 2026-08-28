import types
import unittest
from unittest import mock

import numpy as np

from gritic import multiplicityoptimiser


class FakeMultiplicityProbabilities:
    def __init__(self, major_cn, likelihood):
        self.major_cn = major_cn
        self._likelihood = likelihood

    def evaluate_likelihood(self, state):
        return self._likelihood(np.asarray(state))


class SimplexPointTest(unittest.TestCase):
    def test_generated_point_is_strictly_inside_requested_simplex(self):
        np.random.seed(8128)

        for dimension in (1, 2, 8):
            with self.subTest(dimension=dimension):
                point = multiplicityoptimiser.get_point_on_simplex(dimension)
                self.assertEqual(point.shape, (dimension,))
                self.assertTrue((point > 0.0).all())
                self.assertAlmostEqual(point.sum(), 1.0, 14)

    def test_simplex_sampling_obeys_numpy_seed(self):
        np.random.seed(91)
        first = multiplicityoptimiser.get_point_on_simplex(5)
        np.random.seed(91)
        second = multiplicityoptimiser.get_point_on_simplex(5)

        np.testing.assert_array_equal(first, second)


class MultiplicityOptimisationTest(unittest.TestCase):
    def test_objective_is_negative_log_likelihood(self):
        probabilities = FakeMultiplicityProbabilities(
            2,
            lambda state: 3.5 * state[0] - 2.0 * state[1],
        )

        objective = multiplicityoptimiser.unconstrained_mult_likelihood(
            np.array([0.25, 0.75]),
            probabilities,
        )

        self.assertEqual(objective, -(3.5 * 0.25 - 2.0 * 0.75))

    def test_optimizer_builds_unit_simplex_bounds_and_constraint(self):
        probabilities = FakeMultiplicityProbabilities(3, lambda state: 0.0)
        start = np.array([0.1, 0.2, 0.3, 0.4])
        solution = np.array([0.4, 0.3, 0.2, 0.1])

        with mock.patch.object(
            multiplicityoptimiser,
            'get_point_on_simplex',
            return_value=start,
        ) as point_mock, mock.patch.object(
            multiplicityoptimiser,
            'minimize',
            return_value=types.SimpleNamespace(success=True, x=solution),
        ) as minimize_mock:
            actual = multiplicityoptimiser.unconstrained_mult_optimisation(
                probabilities,
                n_subclones=1,
            )

        point_mock.assert_called_once_with(4)
        self.assertIs(actual, solution)
        positional, keyword = minimize_mock.call_args
        self.assertIs(
            positional[0],
            multiplicityoptimiser.unconstrained_mult_likelihood,
        )
        self.assertIs(positional[1], start)
        self.assertEqual(keyword['bounds'], [(0, 1)] * 4)
        self.assertEqual(keyword['args'], (probabilities,))
        constraint = keyword['constraints']
        np.testing.assert_array_equal(constraint.A, np.ones((1, 4)))
        np.testing.assert_array_equal(constraint.lb, [1])
        np.testing.assert_array_equal(constraint.ub, [1])

    def test_unsuccessful_solver_returns_none(self):
        probabilities = FakeMultiplicityProbabilities(2, lambda state: 0.0)

        with mock.patch.object(
            multiplicityoptimiser,
            'minimize',
            return_value=types.SimpleNamespace(
                success=False,
                x=np.array([0.5, 0.5]),
            ),
        ):
            result = multiplicityoptimiser.unconstrained_mult_optimisation(
                probabilities,
                n_subclones=0,
            )

        self.assertIsNone(result)

    def test_real_solver_recovers_known_concave_mixture_optimum(self):
        weights = np.array([1.0, 2.0, 3.0])
        probabilities = FakeMultiplicityProbabilities(
            3,
            lambda state: np.sum(
                weights * np.log(np.clip(state, 1e-300, None))
            ),
        )
        np.random.seed(1301)

        solution = multiplicityoptimiser.unconstrained_mult_optimisation(
            probabilities,
            n_subclones=0,
        )

        self.assertIsNotNone(solution)
        np.testing.assert_allclose(
            solution,
            weights / weights.sum(),
            rtol=2e-4,
            atol=2e-4,
        )
        self.assertTrue((solution >= 0.0).all())
        self.assertAlmostEqual(solution.sum(), 1.0, 8)


if __name__ == '__main__':
    unittest.main()
