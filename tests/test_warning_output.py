import contextlib
import io
import unittest
import warnings
from unittest import mock

import numpy as np

from gritic import hitandrun, sampletools


class WarningOutputTest(unittest.TestCase):
    def test_hit_and_run_emits_user_warning_without_printing(self):
        expected_samples = np.array([[0.25, 0.75]])
        stdout = io.StringIO()

        with mock.patch.object(
            hitandrun,
            '_hit_and_run',
            return_value=(expected_samples, True),
        ):
            with contextlib.redirect_stdout(stdout):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    samples = hitandrun.hit_and_run(
                        np.zeros((2, 2)),
                        np.array([0.25, 0.75]),
                    )

        self.assertIs(samples, expected_samples)
        self.assertEqual(stdout.getvalue(), '')
        self.assertEqual(len(caught), 1)
        self.assertIs(caught[0].category, UserWarning)
        self.assertIn(
            'process count exceeded 50000',
            str(caught[0].message),
        )

    def test_hit_and_run_preserves_array_return_without_warning(self):
        timing_state = np.array([0.25, 0.75])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            samples = hitandrun.hit_and_run(
                np.empty((timing_state.size, 0)),
                timing_state,
                n_samples=3,
            )

        self.assertEqual(caught, [])
        np.testing.assert_array_equal(
            samples,
            np.tile(timing_state, (3, 1)),
        )

    def test_empty_probability_store_raises_without_printing(self):
        array_store = {
            'Non_Phased': None,
            'Major': None,
            'Minor': None,
            'All': None,
        }
        stdout = io.StringIO()

        with contextlib.redirect_stdout(stdout):
            with self.assertRaisesRegex(ValueError, 'No arrays provided'):
                sampletools.MultProbabilityStore(
                    array_store,
                    array_store.copy(),
                    major_cn=1,
                    minor_cn=0,
                    n_subclones=0,
                )

        self.assertEqual(stdout.getvalue(), '')


if __name__ == '__main__':
    unittest.main()
