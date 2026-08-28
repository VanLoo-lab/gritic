import json
import os
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd

from gritic import gritictimer
from gritic.intervaltools import IntervalSpec


class CloneFractionAndRouteBasicsTest(unittest.TestCase):
    def test_sample_clone_fractions_include_residual_clonal_fraction(self):
        np.testing.assert_array_equal(
            gritictimer.get_sample_clone_fractions(None),
            [1.0],
        )
        table = pd.DataFrame({'Subclone_Fraction': [0.2, 0.3]})
        np.testing.assert_allclose(
            gritictimer.get_sample_clone_fractions(table),
            [0.5, 0.2, 0.3],
        )

    def test_route_initialization_and_reset_clear_all_fit_state(self):
        classifier = gritictimer.RouteClassifier(2, 1, False, 'No_WGD')
        route = next(iter(classifier.routes.values()))
        self.assertEqual(route.short_id, route.route_id[:9])
        self.assertEqual(route.total_cn, 3)
        self.assertIsNone(route.mirror_route_id)
        for attribute in (
            'log_evidence',
            'node_timing',
            'wgd_timing_store',
            'n_events_store',
            'mult_store',
            'raw_samples',
            'unphased_mirror_source',
            'density',
            'density_high',
        ):
            self.assertIsNone(getattr(route, attribute), attribute)
        self.assertTrue(np.isnan(route.run_time))

        route.log_evidence = 5
        route.node_timing = np.ones((1, 1))
        route.run_time = 1
        route.reset_fit()
        self.assertIsNone(route.log_evidence)
        self.assertIsNone(route.node_timing)
        self.assertTrue(np.isnan(route.run_time))

    def test_average_events_and_node_timing_before_and_after_fit(self):
        route = next(iter(
            gritictimer.RouteClassifier(2, 1, False, 'No_WGD').routes.values()
        ))
        self.assertTrue(np.isnan(route.get_average_events('N_Events')))
        self.assertTrue(np.isnan(route.get_node_timing(0)))
        route.n_events_store = {'N_Events': [1, 2, 3]}
        route.node_timing = np.array([
            [0.1, 0.2],
            [0.3, 0.4],
            [0.5, 0.6],
            [0.7, 0.8],
        ])
        self.assertEqual(route.get_average_events('N_Events'), 2)
        np.testing.assert_array_equal(route.get_node_timing(2), [0.5, 0.6])

    def test_cumulative_timing_sums_periods_along_each_component_path(self):
        route = next(iter(
            gritictimer.RouteClassifier(2, 1, False, 'No_WGD').routes.values()
        ))
        periods = np.array([
            [0.2, 0.8, 0.8, 1.0],
            [0.3, 0.7, 0.7, 1.0],
        ])
        cumulative = route.get_cumulative_timing(periods)
        np.testing.assert_allclose(cumulative, np.array([
            [0.2, 0.3],
            [1.0, 1.0],
            [1.0, 1.0],
            [1.0, 1.0],
        ]))

    def test_weighted_arrays_resample_aligned_proposals(self):
        route = next(iter(
            gritictimer.RouteClassifier(2, 1, False, 'No_WGD').routes.values()
        ))
        cumulative = np.array([[0.1, 0.2, 0.3], [1.1, 1.2, 1.3]])
        wgd = np.array([0.4, 0.5, 0.6])
        mult = np.array([[1, 10], [2, 20], [3, 30]])
        indexes = np.array([1, 1, 2])
        with mock.patch.object(
            gritictimer.np.random,
            'choice',
            return_value=indexes,
        ) as choice:
            observed = route.get_weighted_arrays(
                cumulative,
                wgd,
                mult,
                np.array([-1e-81, 2.0, 0.0]),
                n_samples=3,
            )
        np.testing.assert_array_equal(observed[0], cumulative[:, indexes])
        np.testing.assert_array_equal(observed[1], wgd[indexes])
        np.testing.assert_array_equal(observed[2], mult[indexes])
        np.testing.assert_array_equal(choice.call_args.kwargs['p'], [0, 1, 0])

    def test_nan_weights_return_aligned_nan_placeholders(self):
        route = next(iter(
            gritictimer.RouteClassifier(2, 1, False, 'No_WGD').routes.values()
        ))
        observed = route.get_weighted_arrays(
            np.ones((2, 4)),
            np.ones(4),
            np.ones((4, 3)),
            np.array([1.0, np.nan, 2.0, 3.0]),
            n_samples=2,
        )
        self.assertEqual(observed[0].shape, (2, 2))
        self.assertEqual(observed[1].shape, (2,))
        self.assertEqual(observed[2].shape, (2, 3))
        for value in observed:
            self.assertTrue(np.isnan(value).all())

    def test_event_estimate_uses_one_shared_random_sample_index(self):
        route = gritictimer.Route.__new__(gritictimer.Route)
        route.route_tree = SimpleNamespace()
        route.route_tree.get_n_events_batch = mock.Mock(return_value=(
            np.array([1, 2, 3]),
            np.array([0, 0, 0]),
            np.array([1, 0, 1]),
        ))
        node_timing = np.arange(12).reshape(3, 4)
        wgd_timing = np.array([0.1, 0.2, 0.3, 0.4])
        indexes = np.array([3, 1, 3])
        with mock.patch.object(
            gritictimer.np.random,
            'choice',
            return_value=indexes,
        ):
            result = route.get_n_events_estimate(
                node_timing,
                wgd_timing,
                n_samples=3,
            )
        route.route_tree.get_n_events_batch.assert_called_once()
        call = route.route_tree.get_n_events_batch.call_args.args
        np.testing.assert_array_equal(call[0], node_timing[:, indexes])
        np.testing.assert_array_equal(call[1], wgd_timing[indexes])
        self.assertEqual(result, {
            'N_Events': [1, 2, 3],
            'Pre_WGD_Losses': [0, 0, 0],
            'Post_WGD_Losses': [1, 0, 1],
        })


class ProposalGeometryValidationTest(unittest.TestCase):
    def setUp(self):
        classifier = gritictimer.RouteClassifier(2, 1, False, 'No_WGD')
        self.route = next(iter(classifier.routes.values()))
        self.geometry = gritictimer.ProposalGeometry(
            mult_store=np.ones((2, 5)),
            timing_store=np.ones((4, 2)),
            wgd_timing_store=np.full(2, np.nan),
            density=np.array([0.9, 0.8]),
        )

    def test_materialization_rejects_wrong_base_geometry_shape(self):
        for value in (np.ones(5), np.ones((2, 4))):
            with self.subTest(shape=value.shape):
                geometry = gritictimer.ProposalGeometry(
                    mult_store=value,
                    timing_store=self.geometry.timing_store,
                    wgd_timing_store=self.geometry.wgd_timing_store,
                    density=self.geometry.density,
                )
                with self.assertRaisesRegex(ValueError, 'unexpected number'):
                    self.route.materialize_mult_store(
                        geometry,
                        alpha=None,
                        n_subclones=0,
                    )

    def test_materialization_rejects_clone_share_for_clonal_sample(self):
        with self.assertRaisesRegex(ValueError, 'must be omitted'):
            self.route.materialize_mult_store(
                self.geometry,
                alpha=None,
                n_subclones=0,
                clone_share=np.ones((2, 1)),
            )

    def test_materialization_validates_prior_and_clone_share_shapes(self):
        with self.assertRaisesRegex(ValueError, 'one clonal'):
            self.route.materialize_mult_store(
                self.geometry,
                alpha=[1],
                n_subclones=1,
            )
        with self.assertRaisesRegex(ValueError, r'shape \(2, 2\)'):
            self.route.materialize_mult_store(
                self.geometry,
                alpha=[1, 1],
                n_subclones=1,
                clone_share=np.ones((2, 3)),
            )

    def test_wgd_geometry_sampling_keeps_selected_wgd_time_aligned(self):
        route = next(iter(
            gritictimer.RouteClassifier(2, 1, True, 'Only_WGD').routes.values()
        ))
        base_mult = np.arange(10, dtype=float).reshape(2, 5)
        timing = np.arange(8, dtype=float).reshape(4, 2)
        with mock.patch.object(
            gritictimer.np.random,
            'choice',
            return_value=0.4,
        ) as choice, mock.patch.object(
            route,
            'sample_mults',
            return_value=(base_mult, timing),
        ) as sample, mock.patch.object(
            route,
            'get_density_estimate',
            return_value=(0.95, 0.5),
        ):
            geometry = route.run_geometry_sampling(
                np.array([0.2, 0.4, 0.6]),
                samples_per_run=2,
                max_samples=2,
                density_cut_off=0.9,
            )
        choice.assert_called_once()
        sample.assert_called_once_with(0.4, 2)
        np.testing.assert_array_equal(geometry.mult_store, base_mult)
        np.testing.assert_array_equal(geometry.timing_store, timing)
        np.testing.assert_array_equal(geometry.wgd_timing_store, [0.4, 0.4])
        np.testing.assert_array_equal(geometry.density, [0.95, 0.5])

    def test_density_estimate_reports_local_neighbor_proportions(self):
        samples = np.array([[0.0], [0.01], [1.0], [2.0]])
        with mock.patch.dict(os.environ, {'LOKY_MAX_CPU_COUNT': '1'}), warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                message='Could not find the number of physical cores',
                category=UserWarning,
            )
            with mock.patch.object(
                gritictimer.np.random,
                'choice',
                return_value=np.arange(4),
            ):
                density, density_high = self.route.get_density_estimate(
                    samples,
                    n_test_points=4,
                    radius=0.05,
                )
        self.assertEqual(density, 0.5)
        self.assertEqual(density_high, 0)


class RouteClassifierOutputTest(unittest.TestCase):
    def fitted_classifier(self):
        classifier = gritictimer.RouteClassifier(2, 1, False, 'No_WGD')
        route = next(iter(classifier.routes.values()))
        route.log_evidence = 0.0
        route.node_timing = np.array([
            [0.1, 0.4, 0.9],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ])
        route.wgd_timing_store = np.full(3, np.nan)
        route.mult_store = np.ones((3, 5))
        route.raw_samples = {}
        route.n_events_store = {
            'N_Events': [1, 2],
            'Pre_WGD_Losses': [np.nan, np.nan],
            'Post_WGD_Losses': [np.nan, np.nan],
        }
        route.run_time = 1.23456
        route.density = 0.8
        route.density_high = 0.6
        classifier._finalize_fit()
        return classifier, route

    def test_output_routes_require_fit(self):
        classifier = gritictimer.RouteClassifier(2, 1, False, 'No_WGD')
        with self.assertRaisesRegex(ValueError, "hasn't been fit"):
            classifier._get_output_routes()

    def test_output_routes_reject_short_identifier_collision(self):
        classifier = gritictimer.RouteClassifier.__new__(
            gritictimer.RouteClassifier
        )
        first = SimpleNamespace(short_id='same')
        second = SimpleNamespace(short_id='same')
        classifier.routes = {'full-a': first, 'full-b': second}
        classifier.route_probabilities = {'full-a': 0.5, 'full-b': 0.5}
        with self.assertRaisesRegex(ValueError, 'collisions: same'):
            classifier._get_output_routes()

    def test_route_table_reports_probability_events_runtime_and_density(self):
        classifier, route = self.fitted_classifier()
        table = classifier.get_route_table()
        self.assertEqual(table['Route'].tolist(), [route.short_id])
        self.assertEqual(table['Probability'].tolist(), [1.0])
        self.assertEqual(table['Average_N_Events'].tolist(), [1.5])
        self.assertTrue(np.isnan(table.loc[0, 'Average_Pre_WGD_Losses']))
        self.assertTrue(np.isnan(table.loc[0, 'Average_Post_WGD_Losses']))
        self.assertEqual(table['Time'].tolist(), [1.235])
        self.assertEqual(table['Density'].tolist(), [0.8])

    def test_timing_table_reports_phasing_median_and_interval(self):
        classifier, route = self.fitted_classifier()
        table = classifier.get_timing_table(
            interval=IntervalSpec(1),
            rounding_digits=2,
        )
        self.assertEqual(len(table), 1)
        self.assertEqual(table.loc[0, 'Route'], route.short_id)
        self.assertEqual(table.loc[0, 'Node'], 0)
        self.assertEqual(table.loc[0, 'Node_Phasing'], 'Major')
        self.assertEqual(table.loc[0, 'Timing'], 0.4)
        self.assertEqual(table.loc[0, 'Timing_CI_Low'], 0.1)
        self.assertEqual(table.loc[0, 'Timing_CI_High'], 0.9)
        np.testing.assert_array_equal(
            classifier.get_best_timing(),
            np.array([[0.1, 0.4, 0.9]]),
        )

    def test_timing_table_rejects_duplicate_route_node_pairs(self):
        classifier, route = self.fitted_classifier()
        route.route_tree.timeable_nodes = [0, 0]
        with self.assertRaisesRegex(ValueError, 'unique Route and Node'):
            classifier.get_timing_table(interval=IntervalSpec(1))

    def test_timing_tree_labels_cover_gain_wgd_and_terminal_nodes(self):
        classifier = gritictimer.RouteClassifier.__new__(
            gritictimer.RouteClassifier
        )
        route = SimpleNamespace(
            route_tree=SimpleNamespace(
                non_phased_node_order=[10, 11, 12],
                timeable_nodes=[10],
                wgd_nodes=[11],
            ),
            get_node_timing=lambda node: np.array([0.2, 0.4, 0.8]),
        )
        labels = classifier.get_timing_tree_labels(
            route,
            {
                'WGD_Timing': 0.5,
                'WGD_Timing_CI_Low': 0.3,
                'WGD_Timing_CI_High': 0.7,
            },
            gain_interval=IntervalSpec(1),
            rounding_digits=2,
        )
        self.assertEqual(labels[10], '0.4 - [0.2,0.8]')
        self.assertEqual(labels[11], '0.5 - [0.3,0.7]')
        self.assertEqual(labels[12], '')

    def test_shared_route_trees_must_match_classifier_state(self):
        source = gritictimer.RouteClassifier(2, 1, False, 'No_WGD')
        route_trees = {
            route_id: route.route_tree
            for route_id, route in source.routes.items()
        }
        with self.assertRaisesRegex(ValueError, 'must match'):
            gritictimer.RouteClassifier(
                3,
                1,
                False,
                'No_WGD',
                route_trees=route_trees,
            )

    def test_fit_route_classifiers_rejects_mismatched_states_before_sampling(self):
        first = gritictimer.RouteClassifier(2, 1, False, 'No_WGD')
        second = gritictimer.RouteClassifier(3, 1, False, 'No_WGD')
        with self.assertRaisesRegex(ValueError, 'matching copy-number state'):
            gritictimer.fit_route_classifiers(
                [(first, None, None), (second, None, None)],
                None,
            )
        gritictimer.fit_route_classifiers([], None)


class TimerHelperTest(unittest.TestCase):
    def test_classifier_output_byte_estimate_uses_all_retained_arrays(self):
        tree = SimpleNamespace(
            major_cn=3,
            minor_cn=1,
            timeable_nodes=[0, 1],
            non_phased_node_order=[0, 1, 2, 3, 4],
        )
        result = gritictimer.estimate_classifier_output_bytes(
            {'route': tree},
            n_subclones=2,
        )
        expected_raw_values = gritictimer.RAW_ROUTE_SAMPLE_COUNT * (2 + 9 + 2)
        expected_conditional_values = (
            gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT * (5 + 9 + 1)
        )
        self.assertEqual(
            result,
            8 * (expected_raw_values + expected_conditional_values + 900),
        )
        self.assertEqual(gritictimer.estimate_classifier_output_bytes({}, 5), 0)

    def test_add_wgd_info_returns_copy_with_all_fields(self):
        source = pd.DataFrame({'Route': ['r']})
        result = gritictimer.add_wgd_info_to_route_table(
            source,
            {
                'WGD_Timing': 0.5,
                'WGD_Timing_CI_Low': 0.2,
                'WGD_Timing_CI_High': 0.8,
            },
            True,
        )
        self.assertEqual(list(source.columns), ['Route'])
        self.assertEqual(result.loc[0, 'WGD_Timing'], 0.5)
        self.assertEqual(result.loc[0, 'WGD_Timing_CI_Low'], 0.2)
        self.assertEqual(result.loc[0, 'WGD_Timing_CI_High'], 0.8)
        self.assertTrue(result.loc[0, 'WGD_Status'])

    def test_wgd_info_handles_absence_and_rounds_distribution_summary(self):
        missing = gritictimer.get_wgd_info(None)
        self.assertTrue(all(np.isnan(value) for value in missing.values()))
        observed = gritictimer.get_wgd_info(
            np.array([0.1111, 0.5555, 0.9999]),
            interval=IntervalSpec(1),
            rounding_digits=3,
        )
        self.assertEqual(observed, {
            'WGD_Timing': 0.556,
            'WGD_Timing_CI_Low': 0.111,
            'WGD_Timing_CI_High': 1.0,
        })

    def test_potential_wgd_segments_filter_cn_mutations_and_autosomes(self):
        segments = [
            SimpleNamespace(major_cn=2, n_mutations=10, chromosome='1'),
            SimpleNamespace(major_cn=2, n_mutations=9, chromosome='1'),
            SimpleNamespace(major_cn=3, n_mutations=20, chromosome='1'),
            SimpleNamespace(major_cn=2, n_mutations=20, chromosome='X'),
        ]
        sample = SimpleNamespace(autosomes=['1', '2'], segments=segments)
        self.assertEqual(
            gritictimer.get_potential_wgd_segments(sample),
            [segments[0]],
        )

    def test_permitted_copy_number_state_boundaries(self):
        cases = (
            ((1, 0, False), True),
            ((1, 1, True), True),
            ((2, 0, False), True),
            ((2, 2, True), False),
            ((3, 3, True), True),
            ((5, 5, True), True),
            ((6, 4, True), True),
            ((6, 5, True), False),
            ((7, 3, False), True),
            ((7, 4, False), False),
            ((8, 1, True), True),
            ((8, 2, True), False),
            ((9, 0, False), False),
            ((2, 3, False), False),
            ((2, -1, False), False),
            ((0, 0, False), False),
        )
        for arguments, expected in cases:
            with self.subTest(arguments=arguments):
                self.assertIs(
                    gritictimer.check_permitted_cn_state(*arguments),
                    expected,
                )

    def test_timeable_segments_groups_and_filters_in_input_order(self):
        segments = [
            SimpleNamespace(
                segment_id='keep-a', major_cn=3, minor_cn=1, n_mutations=5
            ),
            SimpleNamespace(
                segment_id='exclude', major_cn=3, minor_cn=1, n_mutations=9
            ),
            SimpleNamespace(
                segment_id='few', major_cn=4, minor_cn=2, n_mutations=1
            ),
            SimpleNamespace(
                segment_id='keep-b', major_cn=3, minor_cn=1, n_mutations=6
            ),
            SimpleNamespace(
                segment_id='invalid', major_cn=8, minor_cn=2, n_mutations=20
            ),
        ]
        result = gritictimer.get_timeable_segments(
            SimpleNamespace(segments=segments),
            wgd_status=False,
            excluded_segment_ids=['exclude'],
            min_mutations=5,
        )
        self.assertEqual(list(result), [(3, 1)])
        self.assertEqual(result[(3, 1)], [segments[0], segments[3]])

    def test_json_scalar_normalizes_numpy_and_nonfinite_values(self):
        cases = (
            (np.int64(3), 3),
            (np.float64(0.5), 0.5),
            (np.nan, None),
            (np.inf, None),
            (pd.NA, None),
            (None, None),
            ('value', 'value'),
        )
        for value, expected in cases:
            with self.subTest(value=value):
                self.assertEqual(gritictimer._json_scalar(value), expected)

    def test_write_wgd_calling_info_emits_strict_json(self):
        with tempfile.TemporaryDirectory() as directory:
            gritictimer.write_wgd_calling_info(
                {
                    'WGD_Timing': np.float64(0.5),
                    'WGD_Timing_CI_Low': np.nan,
                    'WGD_Timing_CI_High': np.inf,
                },
                overlap_proportion=np.float64(0.7),
                best_overlap_timing=np.nan,
                major_cn_mode=np.int64(2),
                wgd_status=True,
                output_dir=Path(directory),
                sample_id='sample',
            )
            path = Path(directory) / 'sample_wgd_calling_info.json'
            text = path.read_text(encoding='utf-8')
        self.assertTrue(text.endswith('\n'))
        self.assertNotIn('NaN', text)
        self.assertNotIn('Infinity', text)
        self.assertEqual(json.loads(text), {
            'WGD_Timing': 0.5,
            'WGD_Timing_CI_Low': None,
            'WGD_Timing_CI_High': None,
            'Major_CN_Mode': 2,
            'Overlap_Proportion': 0.7,
            'WGD_Status': True,
            'Best_Overlap_Timing': None,
        })

    def test_wgd_count_validation(self):
        self.assertIsNone(gritictimer._validate_wgd_count(None))
        for value in (0, 1, np.int64(1)):
            self.assertEqual(gritictimer._validate_wgd_count(value), int(value))
        for value in (-1, 2, 0.0, True, np.bool_(False), '1'):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, '0 or 1'):
                    gritictimer._validate_wgd_count(value)

    def test_minimum_wgd_overlap_validation(self):
        for value, expected in ((0, 0.0), (0.5, 0.5), (np.float64(1), 1.0)):
            self.assertEqual(
                gritictimer._validate_min_wgd_overlap(value),
                expected,
            )
        for value in (-0.1, 1.1, np.nan, np.inf, True, '0.5'):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
                    gritictimer._validate_min_wgd_overlap(value)


if __name__ == '__main__':
    unittest.main()
