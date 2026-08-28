import json
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd

from gritic import gritictimer, sampletools, timingio


def _candidate_timing_table(segment):
    """Return the minimum realistic table consumed by WGD orchestration."""
    return pd.DataFrame({
        'Route': [f'route-{segment.segment_id}'],
        'Node': [0],
        'Node_Phasing': ['Major'],
        'Timing': [0.5],
        'Timing_CI_Low': [0.4],
        'Timing_CI_High': [0.6],
        'Probability': [1.0],
        'Average_Pre_WGD_Losses': [0.0],
        'Average_Post_WGD_Losses': [0.0],
        'Segment_ID': [segment.segment_id],
        'Minor_CN': [segment.minor_cn],
    })


def _fake_candidate(segment_id, minor_cn, width):
    return SimpleNamespace(
        segment_id=segment_id,
        major_cn=2,
        minor_cn=minor_cn,
        width=width,
        subclone_table=None,
    )


class AutomaticWGDIntegrationTest(unittest.TestCase):
    @staticmethod
    def make_sample():
        copy_number_table = pd.DataFrame({
            'Chromosome': ['1', '1'],
            'Segment_Start': [0, 1_000],
            'Segment_End': [1_000, 1_200],
            'Major_CN': [2, 3],
            'Minor_CN': [1, 1],
        })
        mutation_table = pd.DataFrame({
            'Chromosome': ['1'] * 24,
            'Position': (
                list(range(10, 22)) + list(range(1_010, 1_022))
            ),
            'Tumor_Ref_Count': [20] * 24,
            'Tumor_Alt_Count': [5] * 24,
        })
        return sampletools.Sample(
            mutation_table,
            copy_number_table,
            None,
            sample_id='AUTO_WGD_COMPLEX',
            purity=0.8,
            merge_cn=False,
        )

    def test_automatic_call_times_a_real_wgd_constrained_3_plus_1_gain(self):
        sample = self.make_sample()
        original_geometry_sampling = (
            gritictimer.Route.run_geometry_sampling
        )
        original_sample_mults = gritictimer.Route.sample_mults
        constrained_geometry_calls = []

        def bounded_geometry_sampling(
            route,
            wgd_timing_distribution,
            *args,
            **kwargs,
        ):
            return original_geometry_sampling(
                route,
                wgd_timing_distribution,
                samples_per_run=32,
                max_samples=32,
                density_cut_off=0.0,
            )

        def track_constraint_use(route, wgd_timing, n_samples):
            if route.wgd_status and route.major_cn == 3:
                combined_matrix, combined_sum = (
                    route.route_tree.get_combined_constraints(wgd_timing)
                )
                wgd_matrix = route.route_tree.wgd_constraint_matrix
                constrained_geometry_calls.append((
                    wgd_timing,
                    combined_matrix.shape[0],
                    combined_sum.copy(),
                    0 if wgd_matrix is None else wgd_matrix.shape[0],
                    route.route_tree.sum_constraint_matrix.shape[0],
                ))
            return original_sample_mults(route, wgd_timing, n_samples)

        with tempfile.TemporaryDirectory() as temporary_directory, (
            mock.patch.object(
                gritictimer.Route,
                'run_geometry_sampling',
                autospec=True,
                side_effect=bounded_geometry_sampling,
            )
        ), mock.patch.object(
            gritictimer.Route,
            'sample_mults',
            autospec=True,
            side_effect=track_constraint_use,
        ):
            gritictimer.process_sample(
                sample,
                temporary_directory,
                wgd_count=None,
                random_seed=8128,
            )
            output_directory = Path(temporary_directory) / sample.sample_id

            with (
                output_directory
                / 'AUTO_WGD_COMPLEX_wgd_calling_info.json'
            ).open(encoding='utf-8') as input_file:
                calling_info = json.load(input_file)

            self.assertEqual(calling_info['Major_CN_Mode'], 2)
            self.assertIs(calling_info['WGD_Status'], True)
            self.assertEqual(calling_info['Overlap_Proportion'], 1.0)
            self.assertTrue(
                all(
                    0 <= calling_info[key] <= 1
                    for key in (
                        'WGD_Timing',
                        'WGD_Timing_CI_Low',
                        'WGD_Timing_CI_High',
                        'Best_Overlap_Timing',
                    )
                )
            )

            candidate_table = pd.read_csv(
                output_directory
                / 'AUTO_WGD_COMPLEX_gain_timing_table_wgd_segments.tsv',
                sep='\t',
            )
            self.assertEqual(candidate_table['Segment_ID'].tolist(), [
                '1-0-1000',
            ])
            self.assertTrue(candidate_table['Intersecting'].all())

            route_table = pd.read_csv(
                output_directory / 'AUTO_WGD_COMPLEX_route_table.tsv',
                sep='\t',
            )
            self.assertFalse(route_table.empty)
            self.assertEqual(
                set(route_table['Segment_ID']),
                {'1-1000-1200'},
            )
            self.assertTrue(route_table['WGD_Status'].all())
            self.assertTrue(
                route_table[
                    ['WGD_Timing', 'WGD_Timing_CI_Low', 'WGD_Timing_CI_High']
                ].apply(np.isfinite).all().all()
            )

            gain_table = pd.read_csv(
                output_directory / 'AUTO_WGD_COMPLEX_gain_timing_table.tsv',
                sep='\t',
            )
            self.assertFalse(gain_table.empty)
            self.assertEqual(set(gain_table['Segment_ID']), {'1-1000-1200'})

            archive_path, manifest_path = timingio.get_timing_archive_paths(
                output_directory / 'AUTO_WGD_COMPLEX_timing_dicts',
                '1-1000-1200',
            )
            timing_hierarchy = timingio.load_timing_archive(
                archive_path,
                manifest_path,
            )
            self.assertEqual(set(timing_hierarchy), set(route_table['Route']))
            for route_timing in timing_hierarchy.values():
                wgd_timing = route_timing['Timing']['WGD']
                self.assertEqual(
                    wgd_timing.shape,
                    (gritictimer.ROUTE_CONDITIONAL_SAMPLE_COUNT,),
                )
                self.assertTrue(np.isfinite(wgd_timing).all())
                self.assertTrue(((0 <= wgd_timing) & (wgd_timing <= 1)).all())

        self.assertTrue(constrained_geometry_calls)
        self.assertTrue(
            all(np.isfinite(observation[0]) for observation in constrained_geometry_calls)
        )
        actual_wgd_constraints = [
            observation
            for observation in constrained_geometry_calls
            if observation[3] > 0
        ]
        self.assertTrue(actual_wgd_constraints)
        for wgd_timing, row_count, constraint_sum, wgd_rows, sum_rows in (
            actual_wgd_constraints
        ):
            self.assertEqual(row_count, sum_rows + wgd_rows)
            np.testing.assert_allclose(
                constraint_sum[-wgd_rows:],
                wgd_timing,
            )


class WGDRunDecisionTest(unittest.TestCase):
    @staticmethod
    def run_sample_with_patches(*, major_cn_mode, wgd_count, wgd_result):
        sample = SimpleNamespace(sample_id='decision-sample')
        patchers = (
            mock.patch.object(
                gritictimer,
                'time_wgd_major_cn_2',
                return_value=wgd_result,
            ),
            mock.patch.object(gritictimer, 'write_wgd_calling_info'),
            mock.patch.object(
                gritictimer,
                'get_timeable_segments',
                return_value={},
            ),
            mock.patch.object(gritictimer, 'process_segments'),
            mock.patch.object(gritictimer.os.path, 'exists', return_value=False),
        )
        active = [patcher.start() for patcher in patchers]
        try:
            gritictimer._run_sample(
                sample,
                Path('output'),
                Path('timing'),
                False,
                0.6,
                wgd_count,
                gritictimer.DEFAULT_TIMING_INTERVALS,
                gritictimer.DEFAULT_SUBCLONE_FRACTION_PRIOR,
                gritictimer.DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
                major_cn_mode,
            )
        finally:
            for patcher in reversed(patchers):
                patcher.stop()
        return active

    def test_automatic_low_overlap_call_falls_back_to_non_wgd(self):
        wgd_distribution = np.array([0.2, 0.4, 0.6])
        with self.assertWarnsRegex(
            UserWarning,
            'overlap proportion is less than 0.6',
        ):
            (
                time_wgd,
                write_call,
                get_segments,
                process_segments,
                _,
            ) = self.run_sample_with_patches(
                major_cn_mode=2,
                wgd_count=None,
                wgd_result=(wgd_distribution, ['excluded'], 0.4, 0.35),
            )

        time_wgd.assert_called_once()
        get_segments.assert_called_once_with(
            mock.ANY,
            wgd_status=False,
            excluded_segment_ids=[],
        )
        process_args = process_segments.call_args.args
        self.assertIsNone(process_args[1])
        self.assertIs(process_args[5], False)
        write_args = write_call.call_args.args
        self.assertEqual(write_args[1:4], (0.4, 0.35, 2))
        self.assertIs(write_args[4], False)

    def test_automatic_mode_one_skips_wgd_timing_and_uses_non_wgd_model(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            (
                time_wgd,
                write_call,
                get_segments,
                process_segments,
                _,
            ) = self.run_sample_with_patches(
                major_cn_mode=1,
                wgd_count=None,
                wgd_result=(np.array([0.5]), ['unused'], 1.0, 0.5),
            )

        self.assertEqual(caught, [])
        time_wgd.assert_not_called()
        get_segments.assert_called_once_with(
            mock.ANY,
            wgd_status=False,
            excluded_segment_ids=[],
        )
        process_args = process_segments.call_args.args
        self.assertIsNone(process_args[1])
        self.assertIs(process_args[5], False)
        write_args = write_call.call_args.args
        self.assertTrue(np.isnan(write_args[1]))
        self.assertTrue(np.isnan(write_args[2]))
        self.assertEqual(write_args[3], 1)
        self.assertIs(write_args[4], False)

    def test_accepted_automatic_call_forwards_distribution_and_exclusions(self):
        wgd_distribution = np.array([0.2, 0.4, 0.6])
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter('always')
            (
                time_wgd,
                write_call,
                get_segments,
                process_segments,
                _,
            ) = self.run_sample_with_patches(
                major_cn_mode=2,
                wgd_count=None,
                wgd_result=(
                    wgd_distribution,
                    ['non-overlapping-segment'],
                    0.8,
                    0.35,
                ),
            )

        self.assertEqual(caught_warnings, [])
        time_wgd.assert_called_once()
        get_segments.assert_called_once_with(
            mock.ANY,
            wgd_status=True,
            excluded_segment_ids=['non-overlapping-segment'],
        )
        process_args = process_segments.call_args.args
        self.assertIs(process_args[1], wgd_distribution)
        self.assertIs(process_args[5], True)
        write_args = write_call.call_args.args
        self.assertEqual(write_args[1:4], (0.8, 0.35, 2))
        self.assertIs(write_args[4], True)

    def test_forced_wgd_warns_when_major_copy_number_mode_is_one(self):
        with self.assertWarnsRegex(
            UserWarning,
            'WGD count 1 but major CN mode is 1',
        ):
            (
                time_wgd,
                _,
                _,
                process_segments,
                _,
            ) = self.run_sample_with_patches(
                major_cn_mode=1,
                wgd_count=1,
                wgd_result=(np.array([0.5]), [], 1.0, 0.5),
            )

        time_wgd.assert_called_once()
        self.assertIs(process_segments.call_args.args[5], True)

    def test_forced_non_wgd_warns_when_major_copy_number_mode_is_two(self):
        with self.assertWarnsRegex(
            UserWarning,
            'WGD count 0 but major CN mode is 2',
        ):
            (
                time_wgd,
                _,
                _,
                process_segments,
                _,
            ) = self.run_sample_with_patches(
                major_cn_mode=2,
                wgd_count=0,
                wgd_result=(np.array([0.5]), [], 1.0, 0.5),
            )

        time_wgd.assert_not_called()
        self.assertIs(process_segments.call_args.args[5], False)


class WGDCandidateOrchestrationTest(unittest.TestCase):
    def test_multi_candidate_pooling_reuses_geometry_and_skips_nonfinite_fit(self):
        segments = [
            _fake_candidate('minor-1-a', minor_cn=1, width=100),
            _fake_candidate('minor-1-nonfinite', minor_cn=1, width=200),
            _fake_candidate('minor-2', minor_cn=2, width=30),
            _fake_candidate('minor-0', minor_cn=0, width=40),
        ]
        sample = SimpleNamespace(
            sample_id='multi',
            subclone_table=None,
            purity=0.75,
        )
        classifiers = []
        constructor_route_trees = []

        def classifier_factory(
            major_cn,
            minor_cn,
            wgd_status,
            wgd_trees_status,
            **kwargs,
        ):
            route_trees = kwargs['route_trees']
            constructor_route_trees.append(route_trees)
            if route_trees is None:
                route_trees = {f'pseudo-{minor_cn}': object()}
            classifier = SimpleNamespace(
                major_cn=major_cn,
                minor_cn=minor_cn,
                routes={
                    route_id: SimpleNamespace(route_tree=route_tree)
                    for route_id, route_tree in route_trees.items()
                },
            )
            classifiers.append(classifier)
            return classifier

        timing_by_id = {
            'minor-1-a': np.full(8, 0.3),
            'minor-1-nonfinite': np.full(8, np.nan),
            'minor-2': np.full(8, 0.3),
            'minor-0': np.full(8, 0.8),
        }

        def candidate_result(segment, classifier, interval_config):
            return (
                _candidate_timing_table(segment),
                timing_by_id[segment.segment_id],
            )

        observed_overlap_inputs = {}

        def maximum_overlap(segment_ci_store, segment_width_store):
            observed_overlap_inputs['ci'] = dict(segment_ci_store)
            observed_overlap_inputs['width'] = dict(segment_width_store)
            return ['minor-1-a', 'minor-2'], 130, 0.3

        pooled_distributions = [
            np.array([0.25, 0.3]),
            np.array([0.28, 0.32]),
        ]
        combined_distribution = np.array([0.275, 0.31])

        with tempfile.TemporaryDirectory() as directory, mock.patch.object(
            gritictimer,
            'get_potential_wgd_segments',
            return_value=segments,
        ), mock.patch.object(
            gritictimer,
            '_get_wgd_timing_model',
            side_effect=lambda segment: (
                0 if segment.minor_cn == 2 else segment.minor_cn,
                f'multiplicity-{segment.segment_id}',
            ),
        ), mock.patch.object(
            gritictimer,
            'RouteClassifier',
            side_effect=classifier_factory,
        ) as route_classifier, mock.patch.object(
            gritictimer,
            'fit_route_classifiers',
        ) as fit_classifiers, mock.patch.object(
            gritictimer,
            '_get_wgd_segment_result',
            side_effect=candidate_result,
        ), mock.patch.object(
            gritictimer.distributiontools,
            'get_ids_with_maximum_overlap',
            side_effect=maximum_overlap,
        ), mock.patch.object(
            gritictimer,
            'get_combined_segment_timing_cn_2',
            return_value=pooled_distributions,
        ) as pool_segments, mock.patch.object(
            gritictimer,
            'get_combined_distribution',
            return_value=combined_distribution,
        ) as combine_distributions:
            result = gritictimer.time_wgd_major_cn_2(
                sample,
                Path(directory),
                Path(directory) / 'timing',
            )
            candidate_output = pd.read_csv(
                Path(directory) / 'multi_gain_timing_table_wgd_segments.tsv',
                sep='\t',
            )

        np.testing.assert_array_equal(result[0], combined_distribution)
        self.assertEqual(result[1], [
            'minor-1-nonfinite',
            'minor-0',
        ])
        self.assertAlmostEqual(result[2], 130 / 170)
        self.assertEqual(result[3], 0.3)

        self.assertEqual(route_classifier.call_count, 4)
        self.assertIsNone(constructor_route_trees[0])
        self.assertIsNone(constructor_route_trees[2])
        self.assertEqual(
            constructor_route_trees[1],
            {
                route_id: route.route_tree
                for route_id, route in classifiers[0].routes.items()
            },
        )
        self.assertEqual(
            constructor_route_trees[3],
            {
                route_id: route.route_tree
                for route_id, route in classifiers[2].routes.items()
            },
        )
        self.assertEqual(fit_classifiers.call_count, 2)
        self.assertEqual(
            [len(call.args[0]) for call in fit_classifiers.call_args_list],
            [2, 2],
        )
        for call in fit_classifiers.call_args_list:
            self.assertIsNone(call.args[1])

        self.assertEqual(
            set(observed_overlap_inputs['ci']),
            {'minor-1-a', 'minor-2', 'minor-0'},
        )
        self.assertEqual(observed_overlap_inputs['width'], {
            'minor-1-a': 100,
            'minor-2': 30,
            'minor-0': 40,
        })
        self.assertEqual(
            pool_segments.call_args.args[0],
            [segments[0], segments[2]],
        )
        self.assertIs(pool_segments.call_args.args[1], sample.subclone_table)
        self.assertEqual(pool_segments.call_args.args[2], sample.purity)
        combine_distributions.assert_called_once_with(pooled_distributions)

        self.assertEqual(len(candidate_output), 4)
        self.assertEqual(
            candidate_output['Intersecting'].tolist(),
            [True, False, True, False],
        )
        np.testing.assert_allclose(
            candidate_output['Overlap_Proportion'],
            130 / 170,
        )
        np.testing.assert_allclose(candidate_output['Best_Overlap_Timing'], 0.3)

    def test_no_eligible_candidate_has_sample_specific_error(self):
        sample = SimpleNamespace(
            sample_id='no-candidate',
            autosomes=['1'],
            segments=[SimpleNamespace(
                major_cn=2,
                n_mutations=9,
                chromosome='1',
            )],
        )
        with self.assertRaisesRegex(
            ValueError,
            'No eligible segments.*sample no-candidate',
        ):
            gritictimer.time_wgd_major_cn_2(
                sample,
                Path('output'),
                Path('timing'),
            )

    def test_all_nonfinite_candidate_timings_have_explicit_error(self):
        segment = _fake_candidate('nonfinite', minor_cn=1, width=100)
        sample = SimpleNamespace(sample_id='no-finite')
        classifier = SimpleNamespace(
            routes={'route': SimpleNamespace(route_tree=object())},
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
        ), mock.patch.object(
            gritictimer,
            'fit_route_classifiers',
        ) as fit_classifiers, mock.patch.object(
            gritictimer,
            '_get_wgd_segment_result',
            return_value=(pd.DataFrame(), np.array([np.nan, np.inf])),
        ):
            with self.assertRaisesRegex(
                ValueError,
                'No eligible.*finite WGD timing interval.*sample no-finite',
            ):
                gritictimer.time_wgd_major_cn_2(
                    sample,
                    Path('output'),
                    Path('timing'),
                )
        fit_classifiers.assert_called_once()


class CombinedWGDSegmentTest(unittest.TestCase):
    @staticmethod
    def segment(
        segment_id,
        minor_cn,
        width,
        *,
        min_mutation_alt_count=3,
        coverage_vaf_quantile=0.95,
        apply_reads_correction=True,
    ):
        return SimpleNamespace(
            segment_id=segment_id,
            major_cn=2,
            minor_cn=minor_cn,
            width=width,
            min_mutation_alt_count=min_mutation_alt_count,
            coverage_vaf_quantile=coverage_vaf_quantile,
            apply_reads_correction=apply_reads_correction,
            mutation_table=pd.DataFrame({
                'Mutation_ID': [f'{segment_id}-mutation'],
                'Major_CN': [2],
                'Minor_CN': [minor_cn],
            }),
        )

    def test_pooling_groups_real_minor_cn_and_sums_segment_widths(self):
        segments = [
            self.segment('minor-1-a', 1, 100),
            self.segment('minor-0', 0, 30),
            self.segment('minor-1-b', 1, 40),
        ]
        first_result = np.array([0.1, 0.2])
        second_result = np.array([0.3, 0.4])
        with mock.patch.object(
            gritictimer,
            '_time_combined_wgd_segment',
            side_effect=[first_result, second_result],
        ) as time_combined:
            result = gritictimer.get_combined_segment_timing_cn_2(
                segments,
                subclone_table=None,
                sample_purity=0.8,
                timing_dict_dir=Path('timing'),
            )

        self.assertEqual(result, [first_result, second_result])
        self.assertEqual(time_combined.call_count, 2)
        minor_zero_call, minor_one_call = time_combined.call_args_list
        self.assertEqual(minor_zero_call.args[0], 0)
        self.assertEqual(minor_zero_call.args[2], 30)
        self.assertEqual(
            minor_zero_call.args[1]['Mutation_ID'].tolist(),
            ['minor-0-mutation'],
        )
        self.assertEqual(minor_one_call.args[0], 1)
        self.assertEqual(minor_one_call.args[2], 140)
        self.assertEqual(
            minor_one_call.args[1]['Mutation_ID'].tolist(),
            ['minor-1-a-mutation', 'minor-1-b-mutation'],
        )
        for call in time_combined.call_args_list:
            self.assertIs(call.args[5], True)
            self.assertEqual(call.args[6], 3)
            self.assertEqual(call.args[7], 0.95)

    def test_pooling_rejects_inconsistent_segment_preprocessing_metadata(self):
        cases = (
            (
                {'min_mutation_alt_count': 4},
                'one minimum mutation alternate-read count',
            ),
            (
                {'coverage_vaf_quantile': 0.8},
                'one coverage VAF quantile',
            ),
            (
                {'apply_reads_correction': False},
                'one reads-correction setting',
            ),
        )
        for changed_metadata, message in cases:
            with self.subTest(changed_metadata=changed_metadata):
                segments = [
                    self.segment('first', 1, 100),
                    self.segment('second', 1, 100, **changed_metadata),
                ]
                with self.assertRaisesRegex(ValueError, message):
                    gritictimer.get_combined_segment_timing_cn_2(
                        segments,
                        subclone_table=None,
                        sample_purity=0.8,
                        timing_dict_dir=Path('timing'),
                    )


class ProcessSegmentsBatchingTest(unittest.TestCase):
    def test_empty_copy_number_state_only_initializes_output_tables(self):
        with mock.patch.object(
            gritictimer,
            '_initialize_table',
        ) as initialize_table, mock.patch.object(
            gritictimer,
            'RouteClassifier',
        ) as route_classifier, mock.patch.object(
            gritictimer,
            'fit_route_classifiers',
        ) as fit_classifiers, mock.patch.object(
            gritictimer,
            '_write_segment_results',
        ) as write_results:
            gritictimer.process_segments(
                {(3, 1): []},
                np.array([0.4, 0.6]),
                Path('output'),
                Path('timing'),
                'sample',
                True,
                False,
                {'WGD_Timing': 0.5},
            )

        self.assertEqual(initialize_table.call_count, 2)
        self.assertEqual(
            initialize_table.call_args_list[0].args,
            (
                Path('output/sample_route_table.tsv'),
                gritictimer.ROUTE_TABLE_COLUMNS,
            ),
        )
        self.assertEqual(
            initialize_table.call_args_list[1].args,
            (
                Path('output/sample_gain_timing_table.tsv'),
                gritictimer.GAIN_TIMING_TABLE_COLUMNS,
            ),
        )
        route_classifier.assert_not_called()
        fit_classifiers.assert_not_called()
        write_results.assert_not_called()

    def test_second_same_cn_classifier_reuses_routes_across_real_batches(self):
        segments = [
            SimpleNamespace(
                segment_id=f'segment-{index}',
                multiplicity_probabilities=f'multiplicity-{index}',
                subclone_table=None,
            )
            for index in range(3)
        ]
        classifiers = []

        def classifier_factory(
            major_cn,
            minor_cn,
            wgd_status,
            wgd_trees_status,
            **kwargs,
        ):
            route_trees = kwargs.get('route_trees')
            if route_trees is None:
                route_trees = {'route': object()}
            classifier = SimpleNamespace(
                routes={
                    route_id: SimpleNamespace(route_tree=route_tree)
                    for route_id, route_tree in route_trees.items()
                },
            )
            classifiers.append(classifier)
            return classifier

        with mock.patch.object(
            gritictimer,
            'RouteClassifier',
            side_effect=classifier_factory,
        ) as route_classifier, mock.patch.object(
            gritictimer,
            'estimate_classifier_output_bytes',
            return_value=100,
        ), mock.patch.object(
            gritictimer,
            'SHARED_FIT_OUTPUT_MEMORY_BUDGET',
            100,
        ), mock.patch.object(
            gritictimer,
            '_initialize_table',
        ), mock.patch.object(
            gritictimer,
            'fit_route_classifiers',
        ) as fit_classifiers, mock.patch.object(
            gritictimer,
            '_write_segment_results',
        ) as write_results:
            gritictimer.process_segments(
                {(3, 1): segments},
                np.array([0.4, 0.6]),
                Path('output'),
                Path('timing'),
                'sample',
                True,
                False,
                {'WGD_Timing': 0.5},
            )

        self.assertEqual(route_classifier.call_count, 3)
        shared_route_trees = {
            route_id: route.route_tree
            for route_id, route in classifiers[0].routes.items()
        }
        self.assertNotIn('route_trees', route_classifier.call_args_list[0].kwargs)
        for call in route_classifier.call_args_list[1:]:
            self.assertEqual(call.kwargs['route_trees'], shared_route_trees)
        self.assertEqual(fit_classifiers.call_count, 3)
        for index, call in enumerate(fit_classifiers.call_args_list):
            jobs, wgd_distribution = call.args
            self.assertEqual(len(jobs), 1)
            self.assertIs(jobs[0][0], classifiers[index])
            self.assertEqual(jobs[0][1], f'multiplicity-{index}')
            np.testing.assert_array_equal(wgd_distribution, [0.4, 0.6])
        self.assertEqual(
            [call.args[0].segment_id for call in write_results.call_args_list],
            ['segment-0', 'segment-1', 'segment-2'],
        )


if __name__ == '__main__':
    unittest.main()
