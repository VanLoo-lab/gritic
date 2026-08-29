import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd

from gritic import (
    gritictimer,
    hitandrun,
    posteriortablegen,
    sampletools,
    timingio,
)


class CloneSharePriorValidationTest(unittest.TestCase):
    def test_rejects_invalid_clone_fraction_vectors(self):
        cases = (
            (np.array([]), 'non-empty one-dimensional'),
            (np.array([[0.5, 0.5]]), 'non-empty one-dimensional'),
            (np.array([np.nan, 0.5]), 'finite non-negative'),
            (np.array([np.inf, 0.0]), 'finite non-negative'),
            (np.array([1.1, -0.1]), 'finite non-negative'),
            (np.array([0.4, 0.4]), 'sum to 1'),
        )
        for fractions, message in cases:
            with self.subTest(fractions=fractions):
                with self.assertRaisesRegex(ValueError, message):
                    gritictimer.get_clone_share_prior_alpha(
                        fractions,
                        subclone_fraction_prior='supplied',
                    )

    def test_adjusted_prior_rejects_invalid_correction_arrays(self):
        fractions = np.array([0.7, 0.2, 0.1])
        cases = (
            (None, 'same shape'),
            (np.ones(2), 'same shape'),
            (np.ones((1, 3)), 'same shape'),
            (np.array([1.0, 0.0, 1.0]), 'finite positive'),
            (np.array([1.0, -0.1, 1.0]), 'finite positive'),
            (np.array([1.0, np.nan, 1.0]), 'finite positive'),
            (np.array([1.0, np.inf, 1.0]), 'finite positive'),
        )
        for correction, message in cases:
            with self.subTest(correction=correction):
                with self.assertRaisesRegex(ValueError, message):
                    gritictimer.get_clone_share_prior_alpha(
                        fractions,
                        correction,
                        subclone_fraction_prior='adjusted',
                    )


class DensityRetryTest(unittest.TestCase):
    def setUp(self):
        self.route = next(iter(
            gritictimer.RouteClassifier(2, 1, False, 'No_WGD').routes.values()
        ))

    @staticmethod
    def sample_batch(call_store, clock):
        def sample_mults(_wgd_timing, n_samples):
            if n_samples != 1:
                raise AssertionError('test fixture expects one proposal per batch')
            batch_index = len(call_store)
            call_store.append(batch_index)
            clock.advance(0.01)
            return (
                np.array([[batch_index]], dtype=float),
                np.array([[batch_index]], dtype=float),
            )

        return sample_mults

    def test_density_scheduling_policy(self):
        initial_batch_count = (
            gritictimer._INITIAL_DENSITY_EVALUATION_BATCH_COUNT
        )
        self.assertFalse(gritictimer._density_evaluation_due(
            initial_batch_count - 1,
            sample_count=initial_batch_count - 1,
            max_samples=1_000,
        ))
        self.assertTrue(gritictimer._density_evaluation_due(
            initial_batch_count,
            sample_count=initial_batch_count,
            max_samples=1_000,
        ))
        self.assertFalse(gritictimer._density_evaluation_due(
            initial_batch_count + 1,
            sample_count=initial_batch_count + 1,
            max_samples=1_000,
            current_time=4.0,
            next_evaluation_time=4.0,
        ))
        self.assertTrue(gritictimer._density_evaluation_due(
            initial_batch_count + 2,
            sample_count=initial_batch_count + 2,
            max_samples=1_000,
            current_time=4.01,
            next_evaluation_time=4.0,
        ))
        self.assertTrue(gritictimer._density_evaluation_due(
            batches_sampled=1,
            sample_count=5,
            max_samples=5,
        ))

        self.assertEqual(gritictimer._density_retry_delay(2.0, 0.1), 2.0)
        self.assertEqual(gritictimer._density_retry_delay(0.1, 0.5), 2.5)
        self.assertEqual(
            gritictimer._density_retry_delay(100.0, 100.0),
            gritictimer._MAX_DENSITY_RETRY_DELAY_SECONDS,
        )

    def test_density_miss_is_deferred_then_retried_to_success(self):
        class FakeClock:
            def __init__(self):
                self.now = 0.0

            def __call__(self):
                return self.now

            def advance(self, seconds):
                self.now += seconds

        clock = FakeClock()
        sample_calls = []
        evaluation_batch_counts = []

        def density_estimate(_timing):
            evaluation_batch_counts.append(len(sample_calls))
            clock.advance(0.02)
            if len(evaluation_batch_counts) == 1:
                return 0.2, 0.1
            return 0.95, 0.8

        with mock.patch.object(
            self.route,
            'sample_mults',
            side_effect=self.sample_batch(sample_calls, clock),
        ), mock.patch.object(
            self.route,
            'get_density_estimate',
            side_effect=density_estimate,
        ) as density_estimate, mock.patch.object(
            gritictimer.time,
            'perf_counter',
            side_effect=clock,
        ):
            geometry = self.route.run_geometry_sampling(
                None,
                samples_per_run=1,
                max_samples=250,
                density_cut_off=0.9,
            )

        self.assertEqual(density_estimate.call_count, 2)
        self.assertGreater(
            evaluation_batch_counts[1],
            evaluation_batch_counts[0] + 1,
        )
        self.assertLess(evaluation_batch_counts[1], 250)
        np.testing.assert_array_equal(geometry.density, [0.95, 0.8])
        self.assertEqual(geometry.mult_store.shape[0], len(sample_calls))
        self.assertEqual(geometry.timing_store.shape[1], len(sample_calls))


class HighCopyNumberFitTest(unittest.TestCase):
    def test_real_seven_plus_three_classifier_fit_is_finite_and_normalized(self):
        major_cn = 7
        minor_cn = 3
        probabilities = np.full((5, major_cn), 1.0 / major_cn)
        mult_probabilities = sampletools.MultProbabilityStore(
            {
                'Non_Phased': probabilities,
                'Major': None,
                'Minor': None,
                'All': probabilities,
            },
            {
                'Non_Phased': np.ones(major_cn),
                'Major': None,
                'Minor': None,
                'All': np.ones(major_cn),
            },
            {
                'Non_Phased': np.ones(5, dtype=np.int64),
                'Major': None,
                'Minor': None,
                'All': np.ones(5, dtype=np.int64),
            },
            major_cn=major_cn,
            minor_cn=minor_cn,
            n_subclones=0,
        )
        classifier = gritictimer.RouteClassifier(
            major_cn,
            minor_cn,
            False,
            'No_WGD',
        )
        original_sampling = gritictimer.Route.run_geometry_sampling

        def bounded_real_sampling(route, wgd_timing_distribution):
            return original_sampling(
                route,
                wgd_timing_distribution,
                samples_per_run=40,
                max_samples=40,
                density_cut_off=0.0,
            )

        np.random.seed(1_907)
        hitandrun.seed_random(1_907)
        with mock.patch.object(
            gritictimer.Route,
            'run_geometry_sampling',
            autospec=True,
            side_effect=bounded_real_sampling,
        ), mock.patch.object(
            gritictimer.Route,
            'get_density_estimate',
            autospec=True,
            return_value=(1.0, 1.0),
        ):
            classifier.fit_routes(mult_probabilities, None, None)

        self.assertGreater(len(classifier.routes), 1)
        self.assertEqual(
            set(classifier.route_probabilities),
            set(classifier.routes),
        )
        probabilities = np.asarray(list(classifier.route_probabilities.values()))
        self.assertTrue(np.isfinite(probabilities).all())
        self.assertTrue((probabilities >= 0).all())
        self.assertAlmostEqual(probabilities.sum(), 1.0)
        for route in classifier.routes.values():
            self.assertTrue(np.isfinite(route.log_evidence))
            self.assertEqual(route.node_timing.shape[1], 1_000)
            self.assertTrue(np.isfinite(route.node_timing).all())
            self.assertTrue(
                (
                    (
                        route.node_timing
                        >= -timingio.ROUTE_PARTICLE_TIMING_TOLERANCE
                    )
                    & (
                        route.node_timing
                        <= 1 + timingio.ROUTE_PARTICLE_TIMING_TOLERANCE
                    )
                ).all()
            )

        timing_dict = classifier.get_timing_dict({
            route.short_id: classifier.route_probabilities[route.route_id]
            for route in classifier.routes.values()
        })
        for route_entry in timing_dict.values():
            timing = route_entry[timingio.ROUTE_PARTICLE_TIMING_KEY]
            self.assertTrue((timing >= 0.0).all())
            self.assertTrue((timing <= 1.0).all())


class EmptyProposalValidationTest(unittest.TestCase):
    def test_route_fit_rejects_empty_geometry_before_posterior_sampling(self):
        route = next(iter(
            gritictimer.RouteClassifier(2, 1, False, 'No_WGD').routes.values()
        ))
        geometry = gritictimer.ProposalGeometry(
            mult_store=np.empty((0, 5)),
            timing_store=np.empty((
                len(route.route_tree.non_phased_node_order),
                0,
            )),
            wgd_timing_store=np.empty(0),
            density=np.array([np.nan, np.nan]),
        )
        mult_probabilities = SimpleNamespace(
            use_major=False,
            use_minor=False,
            evaluate_likelihood_array=lambda _states: np.empty(0),
        )

        with self.assertRaisesRegex(
            ValueError,
            'route fit must contain at least one proposal',
        ):
            route.run_sampling(
                mult_probabilities,
                subclone_table=None,
                wgd_timing_distribution=None,
                proposal_geometry=geometry,
            )


class PenalizedProbabilityOverflowTest(unittest.TestCase):
    def test_finite_event_counts_that_overflow_penalty_return_clean_nan(self):
        maximum_finite = np.finfo(float).max
        route_table = pd.DataFrame({
            'Sample_ID': ['sample', 'sample'],
            'Segment_ID': ['segment', 'segment'],
            'Route': ['route-a', 'route-b'],
            'Probability': [0.5, 0.5],
            'Average_N_Events': [maximum_finite, maximum_finite],
        })
        self.assertTrue(np.isfinite(route_table['Average_N_Events']).all())

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = posteriortablegen.add_penalized_probability(route_table)

        self.assertTrue(result['Penalized_Probability'].isna().all())
        runtime_warnings = [
            warning for warning in caught
            if issubclass(warning.category, RuntimeWarning)
        ]
        self.assertEqual(runtime_warnings, [])


class SyntheticScientificOracleTest(unittest.TestCase):
    @staticmethod
    def make_sample():
        # A 4+1 comb history with cumulative major-gain times .20, .50, .75
        # has unphased mutation exposure [2.8, .25, .30, .20] for
        # multiplicities 1..4. Two hundred exact high-depth VAF observations
        # therefore give the rounded counts below. Multiplicity 3 makes the
        # balanced competing history scientifically incompatible with the data.
        multiplicities = np.repeat(np.arange(1, 5), [158, 14, 17, 11])
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Segment_Start': [0, 10_000],
            'Segment_End': [10_000, 11_000],
            'Major_CN': [1, 4],
            'Minor_CN': [1, 1],
        })
        base_mutation_count = 12
        alt_counts = np.concatenate((
            np.full(base_mutation_count, 150),
            multiplicities * 60,
        ))
        coverage = np.full(alt_counts.size, 300)
        mutation_table = pd.DataFrame({
            'Chromosome': ['1'] * alt_counts.size,
            'Position': (
                list(range(10, 10 + base_mutation_count))
                + list(range(10_010, 10_010 + multiplicities.size))
            ),
            'Tumor_Ref_Count': coverage - alt_counts,
            'Tumor_Alt_Count': alt_counts,
        })
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                message='There are no mutations with less than 10 alt reads',
                category=UserWarning,
            )
            return sampletools.Sample(
                mutation_table,
                copy_number_table,
                None,
                sample_id='SYNTHETIC_ORACLE',
                purity=1.0,
                merge_cn=False,
            )

    def test_process_sample_recovers_known_route_and_gain_times(self):
        classifier = gritictimer.RouteClassifier(4, 1, False, 'No_WGD')
        truth_route = next(
            route for route in classifier.routes.values()
            if [
                route.route_tree.node_attributes[node]['Multiplicity']
                for node in route.route_tree.timeable_nodes
            ] == [4, 3, 2]
        )
        truth_by_node = dict(zip(
            truth_route.route_tree.timeable_nodes,
            (0.20, 0.50, 0.75),
        ))

        with tempfile.TemporaryDirectory() as temporary_directory:
            gritictimer.process_sample(
                self.make_sample(),
                temporary_directory,
                wgd_count=0,
                random_seed=739,
            )
            output_directory = Path(temporary_directory) / 'SYNTHETIC_ORACLE'
            route_table = pd.read_csv(
                output_directory / 'SYNTHETIC_ORACLE_route_table.tsv',
                sep='\t',
            )
            timing_table = pd.read_csv(
                output_directory / 'SYNTHETIC_ORACLE_gain_timing_table.tsv',
                sep='\t',
            )

        gained_routes = route_table.loc[
            route_table['Segment_ID'].eq('1-10000-11000')
        ].sort_values('Probability', ascending=False)
        self.assertEqual(gained_routes.iloc[0]['Route'], truth_route.short_id)
        self.assertGreater(gained_routes.iloc[0]['Probability'], 0.99)
        self.assertLess(gained_routes.iloc[1:]['Probability'].sum(), 0.01)

        recovered = timing_table.loc[
            timing_table['Route'].eq(truth_route.short_id)
        ].set_index('Node')
        self.assertEqual(set(recovered.index), set(truth_by_node))
        for node, truth in truth_by_node.items():
            with self.subTest(node=node, truth=truth):
                self.assertAlmostEqual(
                    recovered.loc[node, 'Timing'],
                    truth,
                    delta=0.08,
                )
                self.assertLessEqual(recovered.loc[node, 'Timing_CI_Low'], truth)
                self.assertGreaterEqual(recovered.loc[node, 'Timing_CI_High'], truth)


if __name__ == '__main__':
    unittest.main()
