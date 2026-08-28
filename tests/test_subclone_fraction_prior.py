import contextlib
import inspect
import io
import unittest
from unittest import mock

import numpy as np

from gritic import cli, gritictimer


class SubcloneFractionPriorTest(unittest.TestCase):
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

    def test_modes_and_default_use_new_names(self):
        self.assertEqual(
            gritictimer.SUBCLONE_FRACTION_PRIOR_MODES,
            ('adjusted', 'supplied'),
        )
        self.assertEqual(
            gritictimer.DEFAULT_SUBCLONE_FRACTION_PRIOR,
            'adjusted',
        )

        args = cli.build_parser().parse_args(self.required_cli_arguments())
        self.assertEqual(args.subclone_fraction_prior, 'adjusted')

        args = cli.build_parser().parse_args([
            *self.required_cli_arguments(),
            '--subclone-fraction-prior',
            'supplied',
        ])
        self.assertEqual(args.subclone_fraction_prior, 'supplied')

    def test_old_cli_interface_is_rejected(self):
        parser = cli.build_parser()
        for arguments in (
            ['--subclone-prior', 'corrected'],
            ['--subclone-fraction-prior', 'corrected'],
            ['--subclone-fraction-prior', 'uncorrected'],
        ):
            with self.subTest(arguments=arguments):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        parser.parse_args([
                            *self.required_cli_arguments(),
                            *arguments,
                        ])

    def test_validator_rejects_old_modes(self):
        for mode in gritictimer.SUBCLONE_FRACTION_PRIOR_MODES:
            self.assertEqual(
                gritictimer.validate_subclone_fraction_prior(mode),
                mode,
            )
        for old_mode in ('corrected', 'uncorrected'):
            with self.subTest(old_mode=old_mode):
                with self.assertRaisesRegex(
                    ValueError,
                    'subclone_fraction_prior must be one of: '
                    'adjusted, supplied',
                ):
                    gritictimer.validate_subclone_fraction_prior(old_mode)

    def test_supplied_prior_uses_input_fractions_directly(self):
        clone_fractions = np.array([0.7, 0.2, 0.1])
        alpha = gritictimer.get_clone_share_prior_alpha(
            clone_fractions,
            subclone_fraction_prior='supplied',
        )
        np.testing.assert_allclose(alpha, 1.0 + clone_fractions)

    def test_adjusted_prior_inverse_weights_by_detection_probability(self):
        clone_fractions = np.array([0.7, 0.2, 0.1])
        detection_probabilities = np.array([1.0, 0.5, 0.25])
        expected_fractions = clone_fractions / detection_probabilities
        expected_fractions /= expected_fractions.sum()

        alpha = gritictimer.get_clone_share_prior_alpha(
            clone_fractions,
            detection_probabilities,
            subclone_fraction_prior='adjusted',
        )
        np.testing.assert_allclose(alpha, 1.0 + expected_fractions)

    def test_geometry_sampling_does_not_draw_clone_shares(self):
        classifier = gritictimer.RouteClassifier(
            2,
            1,
            False,
            'No_WGD',
        )
        route = next(iter(classifier.routes.values()))
        n_nodes = len(route.route_tree.non_phased_node_order)
        base_mult = np.arange(15, dtype=float).reshape(3, 5)
        timing = np.arange(n_nodes * 3, dtype=float).reshape(n_nodes, 3)

        with mock.patch.object(
            route,
            'sample_mults',
            return_value=(base_mult, timing),
        ) as sample_mults, mock.patch.object(
            route,
            'get_density_estimate',
            return_value=(1.0, 1.0),
        ), mock.patch.object(
            route,
            'simulate_clone_share',
            side_effect=AssertionError(
                'geometry sampling must be clone-prior independent'
            ),
        ):
            geometry = route.run_geometry_sampling(
                None,
                samples_per_run=3,
                max_samples=3,
            )

        sample_mults.assert_called_once_with(np.nan, 3)
        np.testing.assert_array_equal(geometry.mult_store, base_mult)
        np.testing.assert_array_equal(geometry.timing_store, timing)

    def test_clone_shares_are_materialized_per_segment_without_mutation(self):
        classifier = gritictimer.RouteClassifier(
            2,
            1,
            False,
            'No_WGD',
        )
        route = next(iter(classifier.routes.values()))
        base_mult = np.array([
            [0.1, 0.2, 0.3, 0.4, 0.5],
            [0.6, 0.7, 0.8, 0.9, 1.0],
        ])
        geometry = gritictimer.ProposalGeometry(
            mult_store=base_mult.copy(),
            timing_store=np.zeros((
                len(route.route_tree.non_phased_node_order),
                2,
            )),
            wgd_timing_store=np.full(2, np.nan),
            density=np.array([0.9, 0.8]),
        )
        first_shares = np.array([
            [0.7, 0.2, 0.1],
            [0.6, 0.3, 0.1],
        ])
        second_shares = np.array([
            [0.4, 0.4, 0.2],
            [0.5, 0.2, 0.3],
        ])
        first_alpha = np.array([1.7, 1.2, 1.1])
        second_alpha = np.array([1.4, 1.4, 1.2])

        with mock.patch.object(
            route,
            'simulate_clone_share',
            side_effect=[first_shares, second_shares],
        ) as simulate_clone_share:
            first = route.materialize_mult_store(
                geometry,
                first_alpha,
                n_subclones=2,
            )
            second = route.materialize_mult_store(
                geometry,
                second_alpha,
                n_subclones=2,
            )

        np.testing.assert_array_equal(geometry.mult_store, base_mult)
        np.testing.assert_allclose(
            first[:, :base_mult.shape[1]],
            base_mult * first_shares[:, :1],
        )
        np.testing.assert_allclose(first[:, -2:], first_shares[:, 1:])
        np.testing.assert_allclose(
            second[:, :base_mult.shape[1]],
            base_mult * second_shares[:, :1],
        )
        np.testing.assert_allclose(second[:, -2:], second_shares[:, 1:])
        self.assertFalse(np.array_equal(first, second))
        self.assertEqual(
            simulate_clone_share.call_args_list,
            [mock.call(first_alpha, 2), mock.call(second_alpha, 2)],
        )

    def test_no_subclones_reuses_geometry_multiplicities(self):
        classifier = gritictimer.RouteClassifier(
            2,
            1,
            False,
            'No_WGD',
        )
        route = next(iter(classifier.routes.values()))
        geometry = gritictimer.ProposalGeometry(
            mult_store=np.ones((2, 5)),
            timing_store=np.zeros((
                len(route.route_tree.non_phased_node_order),
                2,
            )),
            wgd_timing_store=np.full(2, np.nan),
            density=np.array([0.9, 0.8]),
        )

        with mock.patch.object(
            route,
            'simulate_clone_share',
            side_effect=AssertionError('no Dirichlet draw is needed'),
        ):
            observed = route.materialize_mult_store(geometry, None, 0)

        self.assertIs(observed, geometry.mult_store)

    def test_programmatic_interface_uses_new_keyword(self):
        parameters = inspect.signature(
            gritictimer.process_sample
        ).parameters
        self.assertIn('subclone_fraction_prior', parameters)
        self.assertNotIn('subclone_prior', parameters)
        with self.assertRaisesRegex(TypeError, 'unexpected keyword argument'):
            gritictimer.process_sample(
                None,
                'output',
                subclone_prior='corrected',
            )

if __name__ == '__main__':
    unittest.main()
