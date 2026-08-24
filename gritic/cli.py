import argparse
import math

import pandas as pd

from gritic import dataloader, gritictimer, sampletools


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


def load_subclone_table(subclone_table_path):
    return pd.read_csv(
        subclone_table_path,
        sep='\t',
        dtype={
            'Cluster': str,
            'Subclone_CCF': float,
            'Subclone_Fraction': float,
        },
    )


def build_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False)
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
        '--max-subclone-ccf',
        type=unit_interval_number,
        default=sampletools.DEFAULT_MAX_SUBCLONE_CCF,
        help=(
            'Maximum Subclone_CCF retained as a subclone (default: 0.9).'
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
            'and Subclone_Fraction.'
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
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

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

    sample = sampletools.Sample(
        mutation_table,
        copy_number_table,
        subclone_table,
        args.sample_id,
        args.purity,
        sex=args.sample_sex,
        merge_cn=args.merge_adjacent_segments,
        min_mutation_alt_count=args.min_mutation_alt_count,
        min_mutation_coverage=args.min_mutation_coverage,
        max_subclone_ccf=args.max_subclone_ccf,
        min_subclone_fraction=args.min_subclone_fraction,
        autosome_count=args.autosome_count,
        drop_unmatched_snvs=args.drop_unmatched_snvs,
        drop_unmatched_chromosomes=args.drop_unmatched_chromosomes,
    )
    gritictimer.process_sample(
        sample,
        args.output,
        plot_trees=args.plot_trees,
        wgd_count=args.wgd_count,
    )


if __name__ == '__main__':
    raise SystemExit(main())
