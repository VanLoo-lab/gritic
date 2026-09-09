import dataclasses
import unittest

import numpy as np

from gritic import distributiontools


class IntervalWidthValidationTest(unittest.TestCase):
    def test_valid_widths_are_returned_as_builtin_floats(self):
        for width in (1, 0.5, np.float64(0.125)):
            with self.subTest(width=width):
                validated = distributiontools.IntervalSpec(width).width
                self.assertIs(type(validated), float)
                self.assertEqual(validated, float(width))

    def test_invalid_widths_are_rejected(self):
        invalid_widths = (
            0,
            -0.1,
            1.00001,
            True,
            False,
            np.bool_(True),
            np.nan,
            np.inf,
            -np.inf,
            '0.9',
            None,
        )
        for width in invalid_widths:
            with self.subTest(width=width):
                with self.assertRaisesRegex(ValueError, 'interval width'):
                    distributiontools.IntervalSpec(width).width


class IntervalConfigurationTest(unittest.TestCase):
    def test_interval_spec_defaults_to_hpd_and_is_frozen(self):
        spec = distributiontools.IntervalSpec(1)

        self.assertEqual(spec.width, 1.0)
        self.assertEqual(spec.method, distributiontools.HPD)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            spec.width = 0.5

    def test_interval_spec_rejects_unknown_method(self):
        with self.assertRaisesRegex(ValueError, 'interval method'):
            distributiontools.IntervalSpec(0.9, 'central-ish')

    def test_timing_config_has_documented_independent_defaults(self):
        config = distributiontools.TimingIntervalConfig()

        self.assertEqual(config.route_gain, distributiontools.IntervalSpec(0.95))
        self.assertEqual(config.tree_gain, distributiontools.IntervalSpec(0.9))
        self.assertEqual(config.wgd_overlap, distributiontools.IntervalSpec(0.9))
        self.assertEqual(config.sample_wgd, distributiontools.IntervalSpec(0.9))
        self.assertEqual(
            config.posterior_summary,
            distributiontools.IntervalSpec(0.95),
        )
        self.assertIsNot(config.route_gain, config.posterior_summary)

    def test_each_timing_config_field_requires_an_interval_spec(self):
        valid = distributiontools.IntervalSpec(0.8)
        field_names = (
            'route_gain',
            'tree_gain',
            'wgd_overlap',
            'sample_wgd',
            'posterior_summary',
        )
        for field_name in field_names:
            values = {name: valid for name in field_names}
            values[field_name] = 0.8
            with self.subTest(field=field_name):
                with self.assertRaisesRegex(
                    TypeError,
                    f'{field_name} must be an IntervalSpec',
                ):
                    distributiontools.TimingIntervalConfig(**values)


class EmpiricalIntervalBoundsTest(unittest.TestCase):
    def test_equal_tailed_interval_uses_symmetric_empirical_quantiles(self):
        bounds = distributiontools.get_interval_bounds(
            [0.0, 10.0],
            distributiontools.IntervalSpec(0.5, distributiontools.EQUAL_TAILED),
        )

        self.assertEqual(bounds, (2.5, 7.5))

    def test_equal_tailed_full_width_spans_the_sample_range(self):
        bounds = distributiontools.get_interval_bounds(
            [4.0, -2.0, 8.0, 1.0],
            distributiontools.IntervalSpec(1.0, distributiontools.EQUAL_TAILED),
        )

        self.assertEqual(bounds, (-2.0, 8.0))

    def test_hpd_selects_the_shortest_contiguous_window(self):
        bounds = distributiontools.get_interval_bounds(
            [101.0, 1.0, 100.0, 0.0, 2.0],
            distributiontools.IntervalSpec(0.5),
        )

        # ArviZ's nearest HDI uses endpoints floor(n * width) indices apart.
        # Here [0, 2] is the shortest three-observation window.
        self.assertEqual(bounds, (0.0, 2.0))

    def test_hpd_tie_selects_the_lowest_sorted_window(self):
        bounds = distributiontools.get_interval_bounds(
            [3.0, 1.0, 0.0, 2.0],
            distributiontools.IntervalSpec(0.5),
        )

        self.assertEqual(bounds, (0.0, 2.0))

    def test_hpd_full_width_and_singleton_return_extrema(self):
        cases = (
            ([3.0, -1.0, 7.0], 1.0, (-1.0, 7.0)),
            ([9.25], 0.01, (9.25, 9.25)),
        )
        for values, width, expected in cases:
            with self.subTest(values=values, width=width):
                self.assertEqual(
                    distributiontools.get_interval_bounds(
                        values,
                        distributiontools.IntervalSpec(width),
                    ),
                    expected,
                )

    def test_duplicate_samples_can_form_zero_width_hpd(self):
        bounds = distributiontools.get_interval_bounds(
            [0.0, 4.0, 4.0, 4.0, 4.0, 10.0],
            distributiontools.IntervalSpec(0.5),
        )

        self.assertEqual(bounds, (4.0, 4.0))

    def test_nonfinite_samples_produce_nan_bounds(self):
        for nonfinite in (np.nan, np.inf, -np.inf):
            with self.subTest(nonfinite=nonfinite):
                low, high = distributiontools.get_interval_bounds(
                    [0.1, nonfinite, 0.9],
                    distributiontools.IntervalSpec(0.95),
                )
                self.assertTrue(np.isnan(low))
                self.assertTrue(np.isnan(high))

    def test_sample_shape_and_emptiness_are_validated(self):
        spec = distributiontools.IntervalSpec(0.9)
        with self.assertRaisesRegex(ValueError, 'one-dimensional'):
            distributiontools.get_interval_bounds([[0.0, 1.0]], spec)
        with self.assertRaisesRegex(ValueError, 'must not be empty'):
            distributiontools.get_interval_bounds([], spec)

    def test_interval_argument_must_be_an_interval_spec(self):
        for interval in (0.9, {'width': 0.9}, None):
            with self.subTest(interval=interval):
                with self.assertRaisesRegex(TypeError, 'IntervalSpec'):
                    distributiontools.get_interval_bounds([0.0, 1.0], interval)


if __name__ == '__main__':
    unittest.main()
