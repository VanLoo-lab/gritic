import argparse
import logging
import math
from contextlib import contextmanager

import pandas as pd

from gritic import dataloader, gritictimer, intervaltools, sampletools


@contextmanager
def _cli_progress_logging():
    logger = gritictimer.logger
    if logger.handlers:
        yield
        return

    previous_level = logger.level
    previous_propagate = logger.propagate
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    try:
        yield
    finally:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(previous_level)
        logger.propagate = previous_propagate


def nonnegative_integer(value):
    try:
        parsed_value = int(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError('must be a non-negative integer') from error
    if parsed_value < 0:
        raise argparse.ArgumentTypeError('must be a non-negative integer')
    return parsed_value


def positive_integer(value):
    try:
        parsed_value = int(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError('must be a positive integer') from error
    if parsed_value <= 0:
        raise argparse.ArgumentTypeError('must be a positive integer')
    return parsed_value


def unit_interval_number(value):
    try:
        parsed_value = float(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(
            'must be a finite number between 0 and 1'
        ) from error
    if not math.isfinite(parsed_value) or not 0 <= parsed_value <= 1:
        raise argparse.ArgumentTypeError(
            'must be a finite number between 0 and 1'
        )
    return parsed_value


def minimum_subclone_ccf(value):
    try:
        parsed_value = float(value)
        return sampletools.validate_min_subclone_ccf(parsed_value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(
            'must be a finite number greater than 0 and at most 1'
        ) from error


def quantile(value):
    try:
        parsed_value = float(value)
        return sampletools.validate_coverage_vaf_quantile(parsed_value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(
            'must be a finite number between 0 and 1'
        ) from error


def interval_width(value):
    try:
        parsed_value = float(value)
        return intervaltools.validate_interval_width(parsed_value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(
            'must be a finite proportion greater than 0 and at most 1'
        ) from error


def _add_interval_arguments(
    parser,
    option_name,
    destination,
    default_interval,
    description,
):
    parser.add_argument(
        f'--{option_name}-interval-width',
        dest=f'{destination}_interval_width',
        type=interval_width,
        default=default_interval.width,
        metavar='PROPORTION',
        help=(
            f'{description} interval width as a proportion greater than 0 '
            'and at most 1 '
            f'(default: {default_interval.width:g}).'
        ),
    )
    parser.add_argument(
        f'--{option_name}-interval-method',
        dest=f'{destination}_interval_method',
        choices=intervaltools.INTERVAL_METHODS,
        default=default_interval.method,
        help=(
            f'{description} interval method; HPD is the shortest contiguous '
            'empirical interval (default: hpd).'
        ),
    )


def build_interval_config(args):
    return intervaltools.TimingIntervalConfig(
        route_gain=intervaltools.IntervalSpec(
            args.route_gain_interval_width,
            args.route_gain_interval_method,
        ),
        tree_gain=intervaltools.IntervalSpec(
            args.tree_gain_interval_width,
            args.tree_gain_interval_method,
        ),
        wgd_overlap=intervaltools.IntervalSpec(
            args.wgd_overlap_interval_width,
            args.wgd_overlap_interval_method,
        ),
        sample_wgd=intervaltools.IntervalSpec(
            args.wgd_timing_interval_width,
            args.wgd_timing_interval_method,
        ),
        posterior_summary=intervaltools.IntervalSpec(
            args.posterior_summary_interval_width,
            args.posterior_summary_interval_method,
        ),
    )


def load_subclone_table(subclone_table_path):
    return pd.read_csv(
        subclone_table_path,
        sep='\t',
        dtype={
            'Cluster': str,
        },
    )


def build_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            'Time copy-number gains from mutation, copy-number, and optional '
            'subclone tables.'
        ),
    )
    parser.add_argument(
        '--mutation-table',
        required=True,
        help=(
            'A tab separated file containing Chromosome, Tumor_Ref_Count, '
            'Tumor_Alt_Count, and either Mutation_ID or Position. Mutation_ID '
            'takes precedence; otherwise the non-negative integer Position '
            'is used to generate GRITIC_Mutation_ID. Position is additionally '
            'required '
            'to assign copy number unless both input tables contain '
            'Segment_ID. '
            'The selected values must be unique within each source segment; '
            'GRITIC emits a sample-unique GRITIC_Mutation_ID.'
        ),
    )
    parser.add_argument(
        '--copy-number-table',
        required=True,
        help=(
            'A tab separated file containing Chromosome, Segment_Start, '
            'Segment_End, Major_CN, Minor_CN and optionally Segment_ID. Copy '
            'numbers must be non-negative integers with Major_CN greater than '
            'or equal to Minor_CN.'
        ),
    )
    parser.add_argument(
        '--purity',
        type=sampletools.validate_purity,
        required=True,
        help=(
            'The purity of the sample as determined by the copy number '
            'caller; must be greater than 0 and at most 1.'
        ),
    )
    parser.add_argument(
        '--sample-id',
        required=True,
        type=sampletools.validate_sample_id,
        help=(
            'A cross-platform-safe filename component used for the sample '
            'output directory and filename prefixes.'
        ),
    )
    parser.add_argument(
        '--output',
        required=True,
        help=(
            'The output directory. GRITIC stores sample output in '
            'OUTPUT/SAMPLE_ID.'
        ),
    )
    parser.add_argument(
        '--drop-unmatched-snvs',
        action='store_true',
        help=(
            'Drop mutations that cannot be matched to a copy-number segment, '
            'whether matching by Segment_ID or Position, and emit one warning '
            'with the number dropped. By default, unmatched mutations raise '
            'an error.'
        ),
    )
    parser.add_argument(
        '--drop-unmatched-chromosomes',
        action='store_true',
        help=(
            'Drop copy-number and mutation rows whose Chromosome is outside '
            'the configured autosomes and sex karyotype, and warn with the '
            'number dropped. By default, unmatched chromosomes raise an '
            'error.'
        ),
    )
    parser.add_argument(
        '--merge-adjacent-segments',
        action='store_true',
        help=(
            'Merge adjacent copy-number segments with identical Major_CN and '
            'Minor_CN. Disabled by default.'
        ),
    )
    parser.add_argument(
        '--min-mutation-alt-count',
        type=nonnegative_integer,
        default=sampletools.DEFAULT_MIN_MUTATION_ALT_COUNT,
        help=(
            'Minimum Tumor_Alt_Count required to retain a mutation '
            '(default: 3).'
        ),
    )
    parser.add_argument(
        '--min-mutation-coverage',
        type=nonnegative_integer,
        default=sampletools.DEFAULT_MIN_MUTATION_COVERAGE,
        help=(
            'Minimum Tumor_Ref_Count + Tumor_Alt_Count required to retain a '
            'mutation (default: 10).'
        ),
    )
    parser.add_argument(
        '--coverage-vaf-quantile',
        type=quantile,
        default=sampletools.DEFAULT_COVERAGE_VAF_QUANTILE,
        metavar='QUANTILE',
        help=(
            'Observed-SNV VAF quantile between 0 and 1 used to select '
            'mutations for the '
            'mean-coverage estimate in the exact Poisson-thinning detection '
            'correction (default: 0.9). '
            'This high-VAF coverage heuristic differs from the publication/'
            'thesis method, which averages detection power over every SNV\'s '
            'observed depth in the segment; 0.9 is the unit-interval form of '
            'the historical code value (the 90th percentile). '
            'It affects the likelihood correction in both subclone-fraction '
            'prior modes and the prior adjustment only in adjusted mode.'
        ),
    )
    parser.add_argument(
        '--min-subclone-ccf',
        type=minimum_subclone_ccf,
        default=sampletools.DEFAULT_MIN_SUBCLONE_CCF,
        help=(
            'Minimum Subclone_CCF retained as a subclone, inclusive; it must '
            'not exceed --max-subclone-ccf (default: 0.01).'
        ),
    )
    parser.add_argument(
        '--max-subclone-ccf',
        type=unit_interval_number,
        default=sampletools.DEFAULT_MAX_SUBCLONE_CCF,
        help=(
            'Maximum Subclone_CCF retained as a subclone, inclusive; it must '
            'not be below --min-subclone-ccf (default: 0.9).'
        ),
    )
    parser.add_argument(
        '--min-subclone-fraction',
        type=unit_interval_number,
        default=sampletools.DEFAULT_MIN_SUBCLONE_FRACTION,
        help=(
            'Minimum normalized share of subclonal mutations required to '
            'retain a subclone; retention is strictly above this threshold '
            '(default: 0.1).'
        ),
    )
    parser.add_argument(
        '--clip-subclone-ccf',
        action='store_true',
        help=(
            'Clip finite Subclone_CCF values outside [0, 1] to the nearest '
            'boundary before validation and filtering. Non-numeric and '
            'non-finite values remain errors. Disabled by default.'
        ),
    )
    parser.add_argument(
        '--wgd-count',
        type=int,
        choices=(0, 1),
        default=None,
        help=(
            'Override the inferred whole-genome-duplication count with 0 or '
            '1. By default, GRITIC infers the count.'
        ),
    )
    parser.add_argument(
        '--plot-trees',
        action='store_true',
        help='Plot copy number trees. Default is False.',
    )
    parser.add_argument(
        '--subclone-table',
        default=None,
        help=(
            'An optional tab separated file containing Cluster, Subclone_CCF '
            'and Subclone_Fraction. Subclone_CCF is cellular prevalence; '
            'Subclone_Fraction is the fraction of called/input SNVs assigned '
            'to the subclone.'
        ),
    )
    parser.add_argument(
        '--subclone-fraction-prior',
        choices=gritictimer.SUBCLONE_FRACTION_PRIOR_MODES,
        default=gritictimer.DEFAULT_SUBCLONE_FRACTION_PRIOR,
        help=(
            'How supplied Subclone_Fraction values are used to construct the '
            'prior over clonal and subclonal mutation shares. adjusted '
            '(default) divides the supplied called-SNV shares by estimated '
            'detection probabilities separately for each copy-number segment '
            'and renormalizes them; supplied uses the values directly and '
            'reproduces the publication/thesis sample-wide prior. This option '
            'does not modify Subclone_CCF. The likelihood detection correction '
            'remains active per segment in both modes.'
        ),
    )
    parser.add_argument(
        '--sample-sex',
        default=None,
        choices=('XX', 'XY', 'ZZ', 'ZW'),
        help=(
            'Override the sample sex with XX, XY, ZZ, or ZW. If omitted, '
            'GRITIC infers the sex-chromosome system and karyotype from X, Y, '
            'Z, and W copy-number segments. XX accepts X, XY accepts X/Y, ZZ '
            'accepts Z, and ZW accepts Z/W.'
        ),
    )
    parser.add_argument(
        '--autosome-count',
        type=positive_integer,
        default=sampletools.DEFAULT_AUTOSOME_COUNT,
        help=(
            'Number of numbered autosomes in the organism. This defines the '
            'accepted numbered chromosome labels and the chromosomes '
            'eligible for WGD inference (default: 22).'
        ),
    )
    defaults = intervaltools.DEFAULT_TIMING_INTERVALS
    _add_interval_arguments(
        parser,
        'route-gain',
        'route_gain',
        defaults.route_gain,
        'Route-conditional gain-table and WGD-candidate display',
    )
    _add_interval_arguments(
        parser,
        'tree-gain',
        'tree_gain',
        defaults.tree_gain,
        'Blue gain-node tree label',
    )
    _add_interval_arguments(
        parser,
        'wgd-overlap',
        'wgd_overlap',
        defaults.wgd_overlap,
        'Hidden WGD-candidate overlap; changing this can change WGD inference',
    )
    _add_interval_arguments(
        parser,
        'wgd-timing',
        'wgd_timing',
        defaults.sample_wgd,
        'Final sample-level WGD timing',
    )
    _add_interval_arguments(
        parser,
        'posterior-summary',
        'posterior_summary',
        defaults.posterior_summary,
        'Route-marginalized posterior summary',
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.min_subclone_ccf > args.max_subclone_ccf:
        parser.error(
            '--min-subclone-ccf must be less than or equal to '
            '--max-subclone-ccf'
        )
    interval_config = build_interval_config(args)

    copy_number_table, mutation_table = dataloader.load_input_tables(
        args.copy_number_table,
        args.mutation_table,
        drop_unmatched_snvs=args.drop_unmatched_snvs,
        sex=args.sample_sex,
        autosome_count=args.autosome_count,
        drop_unmatched_chromosomes=args.drop_unmatched_chromosomes,
    )
    subclone_table = (
        None
        if args.subclone_table is None
        else load_subclone_table(args.subclone_table)
    )

    sample = sampletools.Sample._from_validated_input_tables(
        mutation_table,
        copy_number_table,
        subclone_table,
        args.sample_id,
        args.purity,
        sex=args.sample_sex,
        merge_cn=args.merge_adjacent_segments,
        min_mutation_alt_count=args.min_mutation_alt_count,
        min_mutation_coverage=args.min_mutation_coverage,
        coverage_vaf_quantile=args.coverage_vaf_quantile,
        min_subclone_ccf=args.min_subclone_ccf,
        max_subclone_ccf=args.max_subclone_ccf,
        min_subclone_fraction=args.min_subclone_fraction,
        clip_subclone_ccf=args.clip_subclone_ccf,
        autosome_count=args.autosome_count,
        drop_unmatched_snvs=args.drop_unmatched_snvs,
        drop_unmatched_chromosomes=args.drop_unmatched_chromosomes,
    )
    with _cli_progress_logging():
        gritictimer.process_sample(
            sample,
            args.output,
            plot_trees=args.plot_trees,
            wgd_count=args.wgd_count,
            interval_config=interval_config,
            subclone_fraction_prior=args.subclone_fraction_prior,
        )


if __name__ == '__main__':
    raise SystemExit(main())
