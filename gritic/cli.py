import argparse

import pandas as pd

from gritic import dataloader, gritictimer, sampletools


def str2bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError('Boolean value expected.')


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
        '--mutation_table',
        required=True,
        help=(
            'A tab separated file containing Chromosome, Tumor_Ref_Count, '
            'Tumor_Alt_Count and either Position or a Segment_ID column also '
            'present in the copy number table'
        ),
    )
    parser.add_argument(
        '--copy_number_table',
        required=True,
        help=(
            'A tab separated file containing Chromosome, Segment_Start, '
            'Segment_End, Major_CN, Minor_CN and optionally Segment_ID. Copy '
            'number should be integer values.'
        ),
    )
    parser.add_argument(
        '--purity',
        type=float,
        required=True,
        help='The purity of the sample as determined by the copy number caller.',
    )
    parser.add_argument('--sample_id', required=True)
    parser.add_argument(
        '--output',
        required=True,
        help=(
            'The output directory. GRITIC stores sample output in '
            'OUTPUT/SAMPLE_ID.'
        ),
    )
    parser.add_argument(
        '--wgd_status',
        type=str2bool,
        default=None,
        help=(
            'Override whole genome duplication status with True or False. '
            'By default, GRITIC determines the status.'
        ),
    )
    parser.add_argument(
        '--non_parsimony_penalty',
        type=str2bool,
        default=False,
        help=(
            'Apply an additional penalty to non-parsimonious copy number '
            'route histories. Default is False.'
        ),
    )
    parser.add_argument(
        '--plot_trees',
        type=str2bool,
        default=True,
        help='Plot copy number trees. Default is True.',
    )
    parser.add_argument(
        '--subclone_table',
        default=None,
        help=(
            'An optional tab separated file containing Cluster, Subclone_CCF '
            'and Subclone_Fraction.'
        ),
    )
    parser.add_argument(
        '--sample_sex',
        default=None,
        help=(
            'Override the sample sex with XX or XY. If omitted, GRITIC '
            'imputes XY when the copy number table contains a Y segment and '
            'XX otherwise.'
        ),
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    copy_number_table, mutation_table = dataloader.load_input_tables(
        args.copy_number_table,
        args.mutation_table,
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
    )
    gritictimer.process_sample(
        sample,
        args.output,
        plot_trees=args.plot_trees,
        wgd_override=args.wgd_status,
        non_parsimony_penalty=args.non_parsimony_penalty,
    )


if __name__ == '__main__':
    raise SystemExit(main())
