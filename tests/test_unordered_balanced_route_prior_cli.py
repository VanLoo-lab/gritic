import contextlib
import inspect
import io
import pathlib
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd

from gritic import cli, gritictimer


class UnorderedBalancedRoutePriorCliTest(unittest.TestCase):
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

    @staticmethod
    def classifier_mock():
        route = SimpleNamespace(
            get_node_timing=lambda _node: np.array([0.5]),
            route_tree=SimpleNamespace(
                major_cn=4,
                minor_cn=4,
                timeable_nodes=[],
                non_phased_node_order=[],
            ),
        )
        return SimpleNamespace(
            routes={'route-id': route},
            fit_routes=mock.Mock(),
            get_route_table=mock.Mock(
                return_value=pd.DataFrame({'Route': ['route-id']})
            ),
            get_timing_table=mock.Mock(
                return_value=pd.DataFrame({'Route': ['route-id']})
            ),
            get_timing_dict=mock.Mock(return_value={}),
            get_best_timing=mock.Mock(return_value=np.array([[0.5]])),
            plot_trees=mock.Mock(),
        )

    def test_parser_defaults_off_and_enables_semantic_flag(self):
        parser = cli.build_parser()
        default_args = parser.parse_args(self.required_cli_arguments())
        self.assertFalse(default_args.unordered_balanced_route_prior)

        enabled_args = parser.parse_args([
            *self.required_cli_arguments(),
            '--unordered-balanced-route-prior',
        ])
        self.assertTrue(enabled_args.unordered_balanced_route_prior)

        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args([
                    *self.required_cli_arguments(),
                    '--unordered_balanced_route_prior',
                ])

    def test_programmatic_process_sample_option_is_strict_boolean(self):
        parameter = inspect.signature(
            gritictimer.process_sample
        ).parameters['unordered_balanced_route_prior']
        self.assertIs(
            parameter.default,
            gritictimer.DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
        )

        with self.assertRaisesRegex(
            ValueError,
            'unordered_balanced_route_prior must be a boolean',
        ):
            gritictimer.process_sample(
                None,
                'output',
                unordered_balanced_route_prior='true',
            )

    def test_prior_reaches_batched_wgd_candidate_classifier(self):
        classifier = self.classifier_mock()
        segment = SimpleNamespace(
            segment_id='segment',
            major_cn=2,
            width=100,
            subclone_table=None,
            get_info_dict=lambda: {},
        )
        sample = SimpleNamespace(
            sample_id='sample',
            subclone_table=None,
            purity=0.8,
        )
        with mock.patch.object(
            gritictimer,
            'get_potential_wgd_segments',
            return_value=[segment],
        ), mock.patch.object(
            gritictimer,
            '_get_wgd_timing_model',
            return_value=(1, object()),
        ), mock.patch.object(
            gritictimer,
            'RouteClassifier',
            return_value=classifier,
        ) as route_classifier, mock.patch.object(
            gritictimer,
            'fit_route_classifiers',
        ), mock.patch.object(
            gritictimer,
            '_get_wgd_segment_result',
            side_effect=RuntimeError('stop after classifier construction'),
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                'stop after classifier construction',
            ):
                gritictimer.time_wgd_major_cn_2(
                    sample,
                    pathlib.Path('output'),
                    pathlib.Path('timing'),
                    unordered_balanced_route_prior=True,
                )

        self.assertIs(
            route_classifier.call_args.kwargs[
                'unordered_balanced_route_prior'
            ],
            True,
        )

    def test_prior_reaches_pooled_wgd_classifier(self):
        classifier = self.classifier_mock()
        pooled_segment = SimpleNamespace(
            major_cn=2,
            subclone_table=None,
        )
        with mock.patch.object(
            gritictimer,
            'Segment',
            return_value=pooled_segment,
        ), mock.patch.object(
            gritictimer,
            '_get_wgd_timing_model',
            return_value=(0, object()),
        ), mock.patch.object(
            gritictimer,
            'RouteClassifier',
            return_value=classifier,
        ) as route_classifier, mock.patch.object(
            gritictimer,
            'write_timing_archive',
        ):
            gritictimer._time_combined_wgd_segment(
                minor_cn=2,
                mutation_table=pd.DataFrame({'Mutation_ID': ['mutation']}),
                combined_width=100,
                subclone_table=None,
                sample_purity=0.8,
                apply_reads_correction=True,
                min_mutation_alt_count=3,
                coverage_vaf_quantile=0.9,
                timing_dict_dir='timing',
                subclone_fraction_prior=(
                    gritictimer.DEFAULT_SUBCLONE_FRACTION_PRIOR
                ),
                unordered_balanced_route_prior=True,
            )

        self.assertIs(
            route_classifier.call_args.kwargs[
                'unordered_balanced_route_prior'
            ],
            True,
        )

    def test_prior_reaches_batched_ordinary_segment_classifier(self):
        classifier = self.classifier_mock()
        segment = SimpleNamespace(
            segment_id='segment',
            major_cn=4,
            minor_cn=4,
            subclone_table=None,
            multiplicity_probabilities=object(),
            get_info_dict=lambda: {},
        )
        with mock.patch.object(
            gritictimer,
            'RouteClassifier',
            return_value=classifier,
        ) as route_classifier, mock.patch.object(
            gritictimer,
            '_initialize_table',
        ), mock.patch.object(
            gritictimer,
            'fit_route_classifiers',
        ), mock.patch.object(
            gritictimer,
            '_write_segment_results',
        ):
            gritictimer.process_segments(
                segments={(4, 4): [segment]},
                wgd_timing_distribution=None,
                output_dir=pathlib.Path('output'),
                timing_dict_dir=pathlib.Path('timing'),
                sample_id='sample',
                wgd_status=False,
                plot_trees=False,
                wgd_info={},
                interval_config=gritictimer.DEFAULT_TIMING_INTERVALS,
                subclone_fraction_prior=(
                    gritictimer.DEFAULT_SUBCLONE_FRACTION_PRIOR
                ),
                unordered_balanced_route_prior=True,
            )

        self.assertIs(
            route_classifier.call_args.kwargs[
                'unordered_balanced_route_prior'
            ],
            True,
        )


if __name__ == '__main__':
    unittest.main()
