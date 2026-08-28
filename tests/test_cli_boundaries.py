import argparse
import contextlib
import io
import logging
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from gritic import cli, gritictimer, intervaltools, sampletools


class CliArgumentFixture:
    @staticmethod
    def required_arguments():
        return [
            '--mutation-table',
            'mutations.tsv',
            '--copy-number-table',
            'copy-number.tsv',
            '--purity',
            '0.8',
            '--sample-id',
            'TEST',
            '--output',
            'output',
        ]

    def assert_parse_error(self, arguments, message):
        parser = cli.build_parser()
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            with self.assertRaisesRegex(SystemExit, '^2$'):
                parser.parse_args(arguments)
        self.assertIn(message, stderr.getvalue())


class CliTypeValidationTest(unittest.TestCase):
    def test_integer_types_accept_boundaries_and_reject_nonintegers(self):
        self.assertEqual(cli.nonnegative_integer('0'), 0)
        self.assertEqual(cli.nonnegative_integer('17'), 17)
        self.assertEqual(cli.positive_integer('1'), 1)
        self.assertEqual(cli.positive_integer('17'), 17)

        for value in ('-1', '1.5', 'nan', '', None):
            with self.subTest(function='nonnegative', value=value):
                with self.assertRaisesRegex(
                    argparse.ArgumentTypeError,
                    'must be a non-negative integer',
                ):
                    cli.nonnegative_integer(value)
        for value in ('0', '-1', '1.5', 'nan', '', None):
            with self.subTest(function='positive', value=value):
                with self.assertRaisesRegex(
                    argparse.ArgumentTypeError,
                    'must be a positive integer',
                ):
                    cli.positive_integer(value)

        self.assertEqual(cli.random_seed('0'), 0)
        self.assertEqual(cli.random_seed(str(2**32 - 1)), 2**32 - 1)
        for value in ('-1', str(2**32), '1.5', 'nan', '', None):
            with self.subTest(function='random_seed', value=value):
                with self.assertRaisesRegex(
                    argparse.ArgumentTypeError,
                    r'integer between 0 and 2\*\*32 - 1',
                ):
                    cli.random_seed(value)

    def test_unit_interval_type_accepts_endpoints_and_rejects_nonfinite_range(self):
        self.assertEqual(cli.unit_interval_number('0'), 0.0)
        self.assertEqual(cli.unit_interval_number('.25'), 0.25)
        self.assertEqual(cli.unit_interval_number('1'), 1.0)

        for value in ('-0.01', '1.01', 'nan', 'inf', '-inf', 'bad', None):
            with self.subTest(value=value):
                with self.assertRaisesRegex(
                    argparse.ArgumentTypeError,
                    'finite number between 0 and 1',
                ):
                    cli.unit_interval_number(value)

    def test_domain_specific_number_types_translate_validation_errors(self):
        valid_cases = (
            (cli.minimum_subclone_ccf, '0.01', 0.01),
            (cli.minimum_subclone_ccf, '1', 1.0),
            (cli.quantile, '0', 0.0),
            (cli.quantile, '1', 1.0),
            (cli.interval_width, '0.01', 0.01),
            (cli.interval_width, '1', 1.0),
        )
        for converter, value, expected in valid_cases:
            with self.subTest(converter=converter.__name__, value=value):
                self.assertEqual(converter(value), expected)

        invalid_cases = (
            (cli.minimum_subclone_ccf, '0', 'greater than 0'),
            (cli.minimum_subclone_ccf, 'nan', 'greater than 0'),
            (cli.minimum_subclone_ccf, None, 'greater than 0'),
            (cli.quantile, '-0.1', 'between 0 and 1'),
            (cli.quantile, 'nan', 'between 0 and 1'),
            (cli.quantile, None, 'between 0 and 1'),
            (cli.interval_width, '0', 'greater than 0'),
            (cli.interval_width, 'inf', 'greater than 0'),
            (cli.interval_width, None, 'greater than 0'),
        )
        for converter, value, message in invalid_cases:
            with self.subTest(converter=converter.__name__, value=value):
                with self.assertRaisesRegex(argparse.ArgumentTypeError, message):
                    converter(value)


class CliParserBoundaryTest(CliArgumentFixture, unittest.TestCase):
    def test_required_arguments_are_enforced(self):
        self.assert_parse_error([], 'the following arguments are required:')
        self.assert_parse_error(
            self.required_arguments()[:-2],
            '--output',
        )

    def test_long_option_abbreviations_and_underscore_aliases_are_rejected(self):
        self.assert_parse_error(
            [*self.required_arguments(), '--plot'],
            'unrecognized arguments: --plot',
        )
        self.assert_parse_error(
            [*self.required_arguments(), '--wgd_count', '1'],
            'unrecognized arguments: --wgd_count 1',
        )

    def test_choice_constrained_options_reject_unknown_values(self):
        cases = (
            (['--sample-sex', 'XO'], "invalid choice: 'XO'"),
            (['--wgd-count', '2'], "invalid choice: '2'"),
            (['--route-gain-interval-method', 'central'], "invalid choice: 'central'"),
            (['--subclone-fraction-prior', 'uniform'], "invalid choice: 'uniform'"),
        )
        for extra_arguments, message in cases:
            with self.subTest(extra_arguments=extra_arguments):
                self.assert_parse_error(
                    [*self.required_arguments(), *extra_arguments],
                    message,
                )

    def test_parser_defaults_match_library_defaults(self):
        args = cli.build_parser().parse_args(self.required_arguments())

        self.assertIsNone(args.sample_sex)
        self.assertEqual(args.autosome_count, sampletools.DEFAULT_AUTOSOME_COUNT)
        self.assertFalse(args.drop_unmatched_chromosomes)
        self.assertFalse(args.drop_unmatched_snvs)
        self.assertFalse(args.drop_unrecognized_phasing)
        self.assertTrue(args.merge_adjacent_segments)
        self.assertEqual(
            args.min_mutation_alt_count,
            sampletools.DEFAULT_MIN_MUTATION_ALT_COUNT,
        )
        self.assertEqual(
            args.min_mutation_coverage,
            sampletools.DEFAULT_MIN_MUTATION_COVERAGE,
        )
        self.assertEqual(
            args.coverage_vaf_quantile,
            sampletools.DEFAULT_COVERAGE_VAF_QUANTILE,
        )
        self.assertIsNone(args.subclone_table)
        self.assertFalse(args.clip_subclone_ccf)
        self.assertEqual(args.min_subclone_ccf, sampletools.DEFAULT_MIN_SUBCLONE_CCF)
        self.assertEqual(args.max_subclone_ccf, sampletools.DEFAULT_MAX_SUBCLONE_CCF)
        self.assertEqual(
            args.min_subclone_fraction,
            sampletools.DEFAULT_MIN_SUBCLONE_FRACTION,
        )
        self.assertEqual(
            args.subclone_fraction_prior,
            gritictimer.DEFAULT_SUBCLONE_FRACTION_PRIOR,
        )
        self.assertIsNone(args.wgd_count)
        self.assertIsNone(args.random_seed)
        self.assertFalse(args.plot_trees)
        self.assertIs(
            args.unordered_balanced_route_prior,
            gritictimer.DEFAULT_UNORDERED_BALANCED_ROUTE_PRIOR,
        )

    def test_interval_configuration_uses_every_custom_width_and_method(self):
        args = cli.build_parser().parse_args([
            *self.required_arguments(),
            '--route-gain-interval-width', '0.81',
            '--route-gain-interval-method', 'equal-tailed',
            '--tree-gain-interval-width', '0.82',
            '--tree-gain-interval-method', 'hpd',
            '--wgd-overlap-interval-width', '0.83',
            '--wgd-overlap-interval-method', 'equal-tailed',
            '--wgd-timing-interval-width', '0.84',
            '--wgd-timing-interval-method', 'hpd',
            '--posterior-summary-interval-width', '0.85',
            '--posterior-summary-interval-method', 'equal-tailed',
        ])

        config = cli.build_interval_config(args)

        self.assertEqual(config.route_gain, intervaltools.IntervalSpec(0.81, 'equal-tailed'))
        self.assertEqual(config.tree_gain, intervaltools.IntervalSpec(0.82, 'hpd'))
        self.assertEqual(config.wgd_overlap, intervaltools.IntervalSpec(0.83, 'equal-tailed'))
        self.assertEqual(config.sample_wgd, intervaltools.IntervalSpec(0.84, 'hpd'))
        self.assertEqual(
            config.posterior_summary,
            intervaltools.IntervalSpec(0.85, 'equal-tailed'),
        )

    def test_subclone_loader_preserves_cluster_text(self):
        table = pd.DataFrame({
            'Cluster': ['001', 'A'],
            'Subclone_CCF': [0.5, 0.25],
            'Subclone_Fraction': [0.2, 0.1],
        })
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'subclones.tsv'
            table.to_csv(path, sep='\t', index=False)
            loaded = cli.load_subclone_table(path)

        self.assertEqual(loaded['Cluster'].tolist(), ['001', 'A'])
        self.assertEqual(loaded['Subclone_CCF'].tolist(), [0.5, 0.25])


class CliProgressLoggingTest(unittest.TestCase):
    def setUp(self):
        self.logger = gritictimer.logger
        self.original_handlers = self.logger.handlers[:]
        self.original_level = self.logger.level
        self.original_propagate = self.logger.propagate
        self.logger.handlers.clear()

    def tearDown(self):
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)
            if handler not in self.original_handlers:
                handler.close()
        self.logger.handlers[:] = self.original_handlers
        self.logger.setLevel(self.original_level)
        self.logger.propagate = self.original_propagate

    def test_context_installs_plain_progress_handler_then_restores_on_exception(self):
        self.logger.setLevel(logging.ERROR)
        self.logger.propagate = True
        stderr = io.StringIO()

        with contextlib.redirect_stderr(stderr):
            with self.assertRaisesRegex(RuntimeError, 'stop'):
                with cli._cli_progress_logging():
                    self.assertEqual(len(self.logger.handlers), 1)
                    self.assertEqual(self.logger.level, logging.INFO)
                    self.assertFalse(self.logger.propagate)
                    self.logger.info('working')
                    raise RuntimeError('stop')

        self.assertEqual(stderr.getvalue(), 'working\n')
        self.assertEqual(self.logger.handlers, [])
        self.assertEqual(self.logger.level, logging.ERROR)
        self.assertTrue(self.logger.propagate)

    def test_existing_handler_is_respected_without_state_changes(self):
        output = io.StringIO()
        existing_handler = logging.StreamHandler(output)
        self.logger.addHandler(existing_handler)
        self.logger.setLevel(logging.WARNING)
        self.logger.propagate = True

        with cli._cli_progress_logging():
            self.assertEqual(self.logger.handlers, [existing_handler])
            self.assertEqual(self.logger.level, logging.WARNING)
            self.assertTrue(self.logger.propagate)
            self.logger.warning('existing')

        self.assertEqual(output.getvalue(), 'existing\n')
        self.assertEqual(self.logger.handlers, [existing_handler])
        self.assertEqual(self.logger.level, logging.WARNING)
        self.assertTrue(self.logger.propagate)


class CliMainTest(CliArgumentFixture, unittest.TestCase):
    def test_cross_argument_subclone_bounds_fail_before_input_loading(self):
        stderr = io.StringIO()
        with mock.patch.object(cli.dataloader, 'load_input_tables') as loader:
            with contextlib.redirect_stderr(stderr):
                with self.assertRaisesRegex(SystemExit, '^2$'):
                    cli.main([
                        *self.required_arguments(),
                        '--min-subclone-ccf', '0.8',
                        '--max-subclone-ccf', '0.7',
                    ])

        self.assertIn(
            '--min-subclone-ccf must be less than or equal to --max-subclone-ccf',
            stderr.getvalue(),
        )
        loader.assert_not_called()

    @mock.patch.object(cli.gritictimer, 'process_sample')
    @mock.patch.object(cli.sampletools.Sample, '_from_validated_input_tables')
    @mock.patch.object(cli, 'load_subclone_table')
    @mock.patch.object(cli.dataloader, 'load_input_tables')
    def test_success_path_forwards_every_input_model_and_interval_option(
        self,
        load_input_tables,
        load_subclone_table,
        make_sample,
        process_sample,
    ):
        copy_number_table = object()
        mutation_table = object()
        subclone_table = object()
        sample = object()
        load_input_tables.return_value = (copy_number_table, mutation_table)
        load_subclone_table.return_value = subclone_table
        make_sample.return_value = sample

        cli.main([
            *self.required_arguments(),
            '--autosome-count', '12',
            '--sample-sex', 'ZW',
            '--drop-unmatched-chromosomes',
            '--drop-unmatched-snvs',
            '--drop-unrecognized-phasing',
            '--no-merge-adjacent-segments',
            '--min-mutation-alt-count', '4',
            '--min-mutation-coverage', '11',
            '--coverage-vaf-quantile', '0.8',
            '--subclone-table', 'subclones.tsv',
            '--clip-subclone-ccf',
            '--min-subclone-ccf', '0.2',
            '--max-subclone-ccf', '0.7',
            '--min-subclone-fraction', '0.15',
            '--subclone-fraction-prior', 'supplied',
            '--wgd-count', '1',
            '--random-seed', '4294967295',
            '--unordered-balanced-route-prior',
            '--plot-trees',
            '--route-gain-interval-width', '0.81',
            '--route-gain-interval-method', 'equal-tailed',
            '--tree-gain-interval-width', '0.82',
            '--tree-gain-interval-method', 'hpd',
            '--wgd-overlap-interval-width', '0.83',
            '--wgd-overlap-interval-method', 'equal-tailed',
            '--wgd-timing-interval-width', '0.84',
            '--wgd-timing-interval-method', 'hpd',
            '--posterior-summary-interval-width', '0.85',
            '--posterior-summary-interval-method', 'equal-tailed',
        ])

        load_input_tables.assert_called_once_with(
            'copy-number.tsv',
            'mutations.tsv',
            drop_unmatched_snvs=True,
            sex='ZW',
            autosome_count=12,
            drop_unmatched_chromosomes=True,
            drop_unrecognized_phasing=True,
        )
        load_subclone_table.assert_called_once_with('subclones.tsv')
        make_sample.assert_called_once_with(
            mutation_table,
            copy_number_table,
            subclone_table,
            'TEST',
            0.8,
            sex='ZW',
            merge_cn=False,
            min_mutation_alt_count=4,
            min_mutation_coverage=11,
            coverage_vaf_quantile=0.8,
            min_subclone_ccf=0.2,
            max_subclone_ccf=0.7,
            min_subclone_fraction=0.15,
            clip_subclone_ccf=True,
            autosome_count=12,
            drop_unmatched_snvs=True,
            drop_unmatched_chromosomes=True,
            drop_unrecognized_phasing=True,
        )
        process_sample.assert_called_once_with(
            sample,
            'output',
            plot_trees=True,
            wgd_count=1,
            interval_config=intervaltools.TimingIntervalConfig(
                route_gain=intervaltools.IntervalSpec(0.81, 'equal-tailed'),
                tree_gain=intervaltools.IntervalSpec(0.82, 'hpd'),
                wgd_overlap=intervaltools.IntervalSpec(0.83, 'equal-tailed'),
                sample_wgd=intervaltools.IntervalSpec(0.84, 'hpd'),
                posterior_summary=intervaltools.IntervalSpec(0.85, 'equal-tailed'),
            ),
            subclone_fraction_prior='supplied',
            unordered_balanced_route_prior=True,
            random_seed=2**32 - 1,
        )

    @mock.patch.object(cli.gritictimer, 'process_sample')
    @mock.patch.object(cli.sampletools.Sample, '_from_validated_input_tables')
    @mock.patch.object(
        cli.dataloader,
        'load_input_tables',
        side_effect=ValueError('invalid input'),
    )
    def test_input_failure_propagates_without_constructing_or_processing_sample(
        self,
        _load_input_tables,
        make_sample,
        process_sample,
    ):
        with self.assertRaisesRegex(ValueError, 'invalid input'):
            cli.main(self.required_arguments())

        make_sample.assert_not_called()
        process_sample.assert_not_called()


if __name__ == '__main__':
    unittest.main()
