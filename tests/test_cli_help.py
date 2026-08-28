import unittest

from gritic import cli


class CliHelpTest(unittest.TestCase):
    EXPECTED_GROUPS = (
        (
            'Required run arguments',
            (
                '--mutation-table',
                '--copy-number-table',
                '--purity',
                '--sample-id',
                '--output',
            ),
        ),
        (
            'Genome and input handling',
            (
                '--autosome-count',
                '--sample-sex',
                '--drop-unmatched-chromosomes',
                '--drop-unmatched-snvs',
                '--drop-unrecognized-phasing',
                '--no-merge-adjacent-segments',
            ),
        ),
        (
            'Mutation filtering and detection correction',
            (
                '--min-mutation-alt-count',
                '--min-mutation-coverage',
                '--coverage-vaf-quantile',
            ),
        ),
        (
            'Subclone handling',
            (
                '--subclone-table',
                '--clip-subclone-ccf',
                '--min-subclone-ccf',
                '--max-subclone-ccf',
                '--min-subclone-fraction',
                '--subclone-fraction-prior',
            ),
        ),
        (
            'Inference model and WGD calling',
            (
                '--wgd-count',
                '--random-seed',
                '--wgd-overlap-interval-width',
                '--wgd-overlap-interval-method',
                '--unordered-balanced-route-prior',
            ),
        ),
        (
            'Reported timing intervals',
            (
                '--route-gain-interval-width',
                '--route-gain-interval-method',
                '--wgd-timing-interval-width',
                '--wgd-timing-interval-method',
                '--posterior-summary-interval-width',
                '--posterior-summary-interval-method',
            ),
        ),
        (
            'Tree plots',
            (
                '--plot-trees',
                '--tree-gain-interval-width',
                '--tree-gain-interval-method',
            ),
        ),
    )

    def test_custom_argument_groups_and_options_have_canonical_order(self):
        parser = cli.build_parser()
        custom_groups = parser._action_groups[2:]
        actual_groups = tuple(
            (
                group.title,
                tuple(
                    action.option_strings[0]
                    for action in group._group_actions
                ),
            )
            for group in custom_groups
        )

        self.assertEqual(actual_groups, self.EXPECTED_GROUPS)

    def test_help_remains_the_only_general_option(self):
        parser = cli.build_parser()
        general_options = tuple(
            tuple(action.option_strings)
            for action in parser._optionals._group_actions
        )

        self.assertEqual(general_options, (('-h', '--help'),))

if __name__ == '__main__':
    unittest.main()
