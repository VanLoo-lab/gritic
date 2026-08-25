import warnings
from decimal import Decimal, InvalidOperation
from numbers import Integral
from urllib.parse import quote

import numpy as np
import pandas as pd

COPY_NUMBER_REQUIRED_COLUMNS = (
    'Chromosome',
    'Segment_Start',
    'Segment_End',
    'Major_CN',
    'Minor_CN',
)
MUTATION_REQUIRED_COLUMNS = (
    'Chromosome',
    'Tumor_Ref_Count',
    'Tumor_Alt_Count',
)
MUTATION_COPY_NUMBER_ANNOTATION_COLUMNS = (
    'Segment_Start',
    'Segment_End',
    'Major_CN',
    'Minor_CN',
)
SUBCLONE_REQUIRED_COLUMNS = (
    'Cluster',
    'Subclone_CCF',
    'Subclone_Fraction',
)
SUBCLONE_OUTPUT_COLUMNS = SUBCLONE_REQUIRED_COLUMNS + ('N_SNVs',)

VALID_SEX_KARYOTYPES = frozenset({'XX', 'XY', 'ZZ', 'ZW'})
SEX_CHROMOSOME_PAIR = {
    'XX': ('X', 'Y'),
    'XY': ('X', 'Y'),
    'ZZ': ('Z', 'W'),
    'ZW': ('Z', 'W'),
}
PRESENT_SEX_CHROMOSOMES = {
    'XX': ('X',),
    'XY': ('X', 'Y'),
    'ZZ': ('Z',),
    'ZW': ('Z', 'W'),
}


def validate_input_table_headers(copy_number_columns, mutation_columns):
    """Validate both input schemas and return whether supplied IDs can be used."""
    copy_number_columns = set(copy_number_columns)
    mutation_columns = set(mutation_columns)

    missing_copy_number_columns = [
        column for column in COPY_NUMBER_REQUIRED_COLUMNS
        if column not in copy_number_columns
    ]
    if missing_copy_number_columns:
        raise ValueError(
            "Copy number table is missing required column(s): "
            + ", ".join(missing_copy_number_columns)
        )

    missing_mutation_columns = [
        column for column in MUTATION_REQUIRED_COLUMNS
        if column not in mutation_columns
    ]
    if missing_mutation_columns:
        raise ValueError(
            "Mutation table is missing required column(s): "
            + ", ".join(missing_mutation_columns)
        )

    use_supplied_segment_ids = (
        'Segment_ID' in copy_number_columns
        and 'Segment_ID' in mutation_columns
    )
    if not {'Mutation_ID', 'Position'} & mutation_columns:
        raise ValueError(
            "Mutation table must contain either Mutation_ID or Position to "
            "identify each mutation"
        )
    if not use_supplied_segment_ids and 'Position' not in mutation_columns:
        raise ValueError(
            "Mutation table must contain Position for copy-number segment "
            "assignment unless Segment_ID is present in both the copy number "
            "and mutation tables"
        )

    return use_supplied_segment_ids


def validate_sex_karyotype(sex):
    """Validate a supported diploid sex-chromosome karyotype."""
    if sex is not None and sex not in VALID_SEX_KARYOTYPES:
        raise ValueError(
            'sex must be None or one of XX, XY, ZZ, or ZW'
        )
    return sex


def validate_autosome_count(autosome_count):
    """Validate and return the number of numbered autosomes."""
    if (
        isinstance(autosome_count, bool)
        or not isinstance(autosome_count, Integral)
        or autosome_count <= 0
    ):
        raise ValueError('autosome_count must be a positive integer')
    return int(autosome_count)


def validate_purity(value):
    """Return a finite purity in the biologically meaningful ``(0, 1]``."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(
            'Purity must be finite, greater than 0, and less than or equal to 1'
        )
    try:
        parsed_value = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(
            'Purity must be finite, greater than 0, and less than or equal to 1'
        ) from error
    if not np.isfinite(parsed_value) or not 0 < parsed_value <= 1:
        raise ValueError(
            'Purity must be finite, greater than 0, and less than or equal to 1'
        )
    return parsed_value


def get_autosome_labels(autosome_count):
    autosome_count = validate_autosome_count(autosome_count)
    return tuple(map(str, range(1, autosome_count + 1)))


def get_normal_total_copy_number(chromosome, sex):
    """Return normal-cell copy number for a chromosome and karyotype."""
    if sex is None:
        return 2
    validate_sex_karyotype(sex)
    first_chromosome, second_chromosome = SEX_CHROMOSOME_PAIR[sex]
    if chromosome == first_chromosome:
        return 2 if sex in {'XX', 'ZZ'} else 1
    if chromosome == second_chromosome:
        return 0 if sex in {'XX', 'ZZ'} else 1
    return 2


def normalize_chromosome_labels(table):
    """Canonicalize chromosome labels, removing one optional ``chr`` prefix."""
    table = table.copy()
    table['Chromosome'] = (
        table['Chromosome']
        .astype(str)
        .str.replace(r'^chr', '', regex=True)
    )
    return table


def get_allowed_chromosome_labels(autosome_count, sex):
    """Return canonical chromosome labels permitted by a genome/karyotype."""
    autosome_count = validate_autosome_count(autosome_count)
    validate_sex_karyotype(sex)
    if sex is None:
        raise ValueError(
            'sex must be resolved before chromosome labels can be validated'
        )
    allowed_chromosomes = set(get_autosome_labels(autosome_count))
    allowed_chromosomes.update(PRESENT_SEX_CHROMOSOMES[sex])
    return allowed_chromosomes


def validate_or_drop_chromosome_rows(
    table,
    autosome_count,
    sex,
    *,
    table_name,
    drop_unmatched_chromosomes=False,
    warn=True,
):
    """Validate chromosome labels, optionally warning and dropping bad rows."""
    table = normalize_chromosome_labels(table)
    allowed_chromosomes = get_allowed_chromosome_labels(
        autosome_count,
        sex,
    )
    normalized_chromosomes = table['Chromosome']
    matched_rows = normalized_chromosomes.isin(allowed_chromosomes)
    invalid_chromosomes = sorted(
        set(normalized_chromosomes.loc[~matched_rows])
    )
    if not invalid_chromosomes:
        return table

    invalid_description = ', '.join(invalid_chromosomes)
    allowed_sex_chromosomes = '/'.join(PRESENT_SEX_CHROMOSOMES[sex])
    requirement = (
        f'numbered autosomes 1 through {autosome_count} or '
        f'{allowed_sex_chromosomes} for sex {sex}'
    )
    if not drop_unmatched_chromosomes:
        raise ValueError(
            f'{table_name} Chromosome values must be {requirement}; '
            f'unmatched value(s): {invalid_description}'
        )
    if warn:
        row_count = int((~matched_rows).sum())
        row_label = 'row' if row_count == 1 else 'rows'
        warnings.warn(
            f'Dropping {row_count} {table_name.lower()} {row_label} with '
            f'unmatched Chromosome value(s): {invalid_description}.',
            UserWarning,
            stacklevel=2,
        )
    return table.loc[matched_rows].copy()


def validate_or_drop_input_chromosomes(
    copy_number_table,
    mutation_table,
    autosome_count,
    sex,
    *,
    drop_unmatched_chromosomes=False,
):
    """Resolve sex, then validate/filter chromosome rows in both inputs."""
    resolved_sex = (
        infer_sex_from_copy_number_table(copy_number_table)
        if sex is None
        else validate_sex_karyotype(sex)
    )
    copy_number_table = validate_or_drop_chromosome_rows(
        copy_number_table,
        autosome_count,
        resolved_sex,
        table_name='Copy number table',
        drop_unmatched_chromosomes=drop_unmatched_chromosomes,
    )
    mutation_table = validate_or_drop_chromosome_rows(
        mutation_table,
        autosome_count,
        resolved_sex,
        table_name='Mutation table',
        drop_unmatched_chromosomes=drop_unmatched_chromosomes,
    )
    if copy_number_table.empty and mutation_table.empty:
        raise ValueError(
            'No copy-number segments or mutations remain after chromosome '
            'filtering'
        )
    if copy_number_table.empty:
        raise ValueError(
            'No copy-number segments remain after chromosome filtering'
        )
    if mutation_table.empty:
        raise ValueError(
            'No mutations remain after chromosome filtering'
        )
    return copy_number_table, mutation_table, resolved_sex


def get_gritic_mutation_id_components(mutation_table):
    """Select the Mutation_ID or canonical Position used for GRITIC IDs."""
    identity_column = (
        'Mutation_ID'
        if 'Mutation_ID' in mutation_table.columns
        else 'Position'
    )
    id_components = mutation_table[identity_column].astype('string')
    missing_ids = (
        id_components.isna()
        | id_components.str.strip().eq('').fillna(True)
    )
    if missing_ids.any():
        raise ValueError(
            f"Every mutation table row must have a {identity_column} value "
            "that is not null, empty, or whitespace-only"
        )
    return id_components.astype(str)


def validate_source_scoped_mutation_components(
    source_segment_ids,
    id_components,
):
    """Require selected identity components to be unique within their source."""
    identity_table = pd.DataFrame({
        'Source_Segment_ID': source_segment_ids.astype(str),
        'Identity_Component': id_components.astype(str),
    })
    duplicate_identities = identity_table[
        identity_table.duplicated(
            subset=['Source_Segment_ID', 'Identity_Component'],
            keep=False,
        )
    ].drop_duplicates()
    if not duplicate_identities.empty:
        duplicate_pairs = ', '.join(
            f"({row.Source_Segment_ID}, {row.Identity_Component})"
            for row in duplicate_identities.itertuples(index=False)
        )
        raise ValueError(
            "The selected Mutation_ID or Position values must be unique "
            "within each source Segment_ID; "
            "duplicate pair(s): " + duplicate_pairs
        )


def generate_gritic_mutation_ids(
    source_segment_ids,
    id_components,
):
    """Build opaque GRITIC IDs from source-scoped mutation IDs."""
    return pd.Series(
        [
            f"{quote(str(source_id), safe='')}:"
            f"{quote(str(mutation_id), safe='')}"
            for source_id, mutation_id in zip(
                source_segment_ids,
                id_components,
            )
        ],
        index=source_segment_ids.index,
        dtype='string',
    )


def parse_positions(positions):
    """Parse non-negative genomic positions without float precision loss."""
    maximum_position = np.iinfo(np.int64).max
    parsed_positions = []
    invalid_values = []
    for value in positions:
        try:
            decimal_value = Decimal(str(value))
            if (
                not decimal_value.is_finite()
                or decimal_value != decimal_value.to_integral_value()
                or not 0 <= decimal_value <= maximum_position
            ):
                raise ValueError
            integer_value = int(decimal_value)
        except (InvalidOperation, ValueError, OverflowError):
            invalid_values.append(str(value))
            parsed_positions.append(0)
        else:
            parsed_positions.append(integer_value)

    if invalid_values:
        unique_invalid_values = list(dict.fromkeys(invalid_values))
        raise ValueError(
            "Position must be a finite non-negative integer within the "
            "signed 64-bit range; invalid value(s): "
            + ", ".join(unique_invalid_values)
        )
    return pd.Series(
        parsed_positions,
        index=positions.index,
        dtype=np.int64,
    )


def validate_segment_coordinates(copy_number_table):
    """Validate nonempty half-open segments and canonicalize coordinates."""
    copy_number_table = copy_number_table.copy()
    minimum_coordinate = np.iinfo(np.int64).min
    maximum_coordinate = np.iinfo(np.int64).max
    violations = []
    parsed_coordinates = {}
    for column in ('Segment_Start', 'Segment_End'):
        parsed_values = []
        invalid_values = []
        for value in copy_number_table[column]:
            try:
                decimal_value = Decimal(str(value))
                if (
                    not decimal_value.is_finite()
                    or decimal_value != decimal_value.to_integral_value()
                    or not minimum_coordinate
                    <= decimal_value
                    <= maximum_coordinate
                ):
                    raise ValueError
                parsed_value = int(decimal_value)
            except (InvalidOperation, ValueError, OverflowError):
                invalid_values.append(str(value))
                parsed_values.append(0)
            else:
                parsed_values.append(parsed_value)

        if invalid_values:
            violations.append(
                f'{column} must contain finite integers within the signed '
                '64-bit range; invalid value(s): '
                + ', '.join(dict.fromkeys(invalid_values))
            )
        parsed_coordinates[column] = pd.Series(
            parsed_values,
            index=copy_number_table.index,
            dtype=np.int64,
        )

    if violations:
        raise ValueError(
            'Invalid copy number table: ' + '; '.join(violations)
        )

    segment_starts = parsed_coordinates['Segment_Start']
    segment_ends = parsed_coordinates['Segment_End']
    if (segment_starts < 0).any():
        violations.append('Segment_Start must be non-negative')
    if (segment_ends < 0).any():
        violations.append('Segment_End must be non-negative')
    if (segment_ends <= segment_starts).any():
        violations.append(
            'Segment_End must be greater than Segment_Start'
        )
    segment_widths = [
        int(segment_end) - int(segment_start)
        for segment_start, segment_end in zip(
            segment_starts,
            segment_ends,
        )
    ]
    if any(width > maximum_coordinate for width in segment_widths):
        violations.append(
            'segment width must not exceed the maximum signed '
            '64-bit integer'
        )
    if violations:
        raise ValueError(
            'Invalid copy number table: ' + '; '.join(violations)
        )
    for column, values in parsed_coordinates.items():
        copy_number_table[column] = values
    return copy_number_table


def get_segment_widths(copy_number_table):
    """Return exact Python-integer widths for validated half-open segments."""
    return pd.Series(
        [
            int(segment_end) - int(segment_start)
            for segment_start, segment_end in zip(
                copy_number_table['Segment_Start'],
                copy_number_table['Segment_End'],
            )
        ],
        index=copy_number_table.index,
        dtype=object,
    )


def validate_mutation_read_counts(mutation_table):
    """Validate and canonicalize mutation read counts as int64 values."""
    mutation_table = mutation_table.copy()
    parsed_counts = {}
    violations = []
    maximum_count = np.iinfo(np.int64).max

    for column in ('Tumor_Ref_Count', 'Tumor_Alt_Count'):
        parsed_values = []
        invalid_values = []
        for value in mutation_table[column]:
            try:
                decimal_value = Decimal(str(value))
                if (
                    not decimal_value.is_finite()
                    or decimal_value != decimal_value.to_integral_value()
                    or not 0 <= decimal_value <= maximum_count
                ):
                    raise ValueError
                parsed_value = int(decimal_value)
            except (InvalidOperation, ValueError, OverflowError):
                invalid_values.append(str(value))
                parsed_values.append(0)
            else:
                parsed_values.append(parsed_value)

        if invalid_values:
            unique_invalid_values = list(dict.fromkeys(invalid_values))
            violations.append(
                f'{column} must contain finite non-negative integers within '
                'the signed 64-bit range; invalid value(s): '
                + ', '.join(unique_invalid_values)
            )
        parsed_counts[column] = pd.Series(
            parsed_values,
            index=mutation_table.index,
            dtype=np.int64,
        )

    if violations:
        raise ValueError('Invalid mutation table: ' + '; '.join(violations))

    total_counts = [
        int(reference_count) + int(alternate_count)
        for reference_count, alternate_count in zip(
            parsed_counts['Tumor_Ref_Count'],
            parsed_counts['Tumor_Alt_Count'],
        )
    ]
    if any(total_count == 0 for total_count in total_counts):
        raise ValueError(
            'Invalid mutation table: Tumor_Ref_Count + Tumor_Alt_Count must '
            'be greater than 0'
        )
    coverage_too_large = any(
        total_count > maximum_count for total_count in total_counts
    )
    if coverage_too_large:
        raise ValueError(
            'Invalid mutation table: Tumor_Ref_Count + Tumor_Alt_Count must '
            'not exceed the maximum signed 64-bit integer'
        )

    for column, parsed_values in parsed_counts.items():
        mutation_table[column] = parsed_values
    return mutation_table


def validate_copy_number_values(copy_number_table):
    """Validate and canonicalize allele-specific copy numbers."""
    copy_number_table = copy_number_table.copy()
    parsed_copy_numbers = {}
    violations = []
    maximum_copy_number = np.iinfo(np.int64).max
    for column in ('Major_CN', 'Minor_CN'):
        parsed_values = []
        invalid_values = []
        for value in copy_number_table[column]:
            try:
                decimal_value = Decimal(str(value))
                if (
                    not decimal_value.is_finite()
                    or decimal_value != decimal_value.to_integral_value()
                    or not -maximum_copy_number - 1
                    <= decimal_value
                    <= maximum_copy_number
                ):
                    raise ValueError
                parsed_value = int(decimal_value)
            except (InvalidOperation, ValueError, OverflowError):
                invalid_values.append(str(value))
                parsed_values.append(0)
            else:
                parsed_values.append(parsed_value)

        if invalid_values:
            unique_invalid_values = list(dict.fromkeys(invalid_values))
            violations.append(
                f'{column} must contain finite integers within the signed '
                '64-bit range; invalid value(s): '
                + ', '.join(unique_invalid_values)
            )
        parsed_copy_numbers[column] = pd.Series(
            parsed_values,
            index=copy_number_table.index,
            dtype=np.int64,
        )

    if violations:
        raise ValueError(
            'Invalid copy number table: ' + '; '.join(violations)
        )

    major_copy_number = parsed_copy_numbers['Major_CN']
    minor_copy_number = parsed_copy_numbers['Minor_CN']
    invalid_major_copy_number = major_copy_number < 0
    invalid_minor_copy_number = minor_copy_number < 0
    inverted_copy_number = minor_copy_number > major_copy_number

    violations = []
    if invalid_major_copy_number.any():
        violations.append('Major_CN must be non-negative')
    if invalid_minor_copy_number.any():
        violations.append('Minor_CN must be non-negative')
    if inverted_copy_number.any():
        violations.append('Minor_CN must not exceed Major_CN')
    if violations:
        raise ValueError(
            'Invalid copy number table: ' + '; '.join(violations)
        )
    copy_number_table['Major_CN'] = major_copy_number
    copy_number_table['Minor_CN'] = minor_copy_number
    return copy_number_table


def validate_subclone_values(subclone_table):
    """Validate and canonicalize the subclone table used as a mixture prior."""
    missing_columns = [
        column for column in SUBCLONE_REQUIRED_COLUMNS
        if column not in subclone_table.columns
    ]
    if missing_columns:
        raise ValueError(
            'Subclone table is missing required column(s): '
            + ', '.join(missing_columns)
        )

    subclone_table = subclone_table.copy()
    for column in ('Subclone_CCF', 'Subclone_Fraction'):
        boolean_values = subclone_table[column].map(
            lambda value: isinstance(value, (bool, np.bool_))
        )
        numeric_values = pd.to_numeric(
            subclone_table[column],
            errors='coerce',
        )
        valid_values = (
            ~boolean_values
            & numeric_values.notna()
            & np.isfinite(numeric_values)
            & numeric_values.between(0, 1)
        )
        if not valid_values.all():
            raise ValueError(
                f'{column} must contain finite values between 0 and 1'
            )
        subclone_table[column] = numeric_values.astype(float)

    total_subclone_fraction = subclone_table['Subclone_Fraction'].sum()
    if total_subclone_fraction > 1:
        if total_subclone_fraction - 1 > 1e-12:
            raise ValueError(
                'Subclone_Fraction values must sum to no more than 1'
            )
        subclone_table['Subclone_Fraction'] /= total_subclone_fraction
    return subclone_table


def validate_supplied_segment_ids(
    copy_number_table,
    mutation_table,
    allow_unmatched=False,
):
    """Validate the value-level contract for supplied segment IDs."""
    copy_number_id_strings = copy_number_table['Segment_ID'].astype('string')
    mutation_id_strings = mutation_table['Segment_ID'].astype('string')

    missing_copy_number_ids = (
        copy_number_id_strings.isna()
        | copy_number_id_strings.str.strip().eq('').fillna(True)
    )
    if missing_copy_number_ids.any():
        raise ValueError(
            "Every copy number table row must have a Segment_ID that is not "
            "null, empty, or whitespace-only"
        )

    missing_mutation_ids = (
        mutation_id_strings.isna()
        | mutation_id_strings.str.strip().eq('').fillna(True)
    )
    if missing_mutation_ids.any():
        raise ValueError(
            "Every mutation table row must have a Segment_ID that is not "
            "null, empty, or whitespace-only"
        )

    copy_number_ids = copy_number_id_strings.astype(str)
    mutation_ids = mutation_id_strings.astype(str)
    duplicate_copy_number_ids = copy_number_ids[
        copy_number_ids.duplicated(keep=False)
    ].unique()
    if duplicate_copy_number_ids.size:
        raise ValueError(
            "Copy number table Segment_ID values must be globally unique; "
            "duplicate value(s): " + ", ".join(duplicate_copy_number_ids)
        )

    matched_mutation_ids = mutation_ids.isin(copy_number_ids)
    unknown_mutation_ids = mutation_ids[~matched_mutation_ids].unique()
    if unknown_mutation_ids.size and not allow_unmatched:
        raise ValueError(
            "Mutation table Segment_ID value(s) not present in the copy "
            "number table: " + ", ".join(unknown_mutation_ids)
        )

    chromosome_by_segment_id = pd.Series(
        copy_number_table['Chromosome'].astype(str).str.replace(
            r'^chr',
            '',
            regex=True,
        ).to_numpy(),
        index=copy_number_ids,
    )
    mutation_chromosomes = mutation_table['Chromosome'].astype(str).str.replace(
        r'^chr',
        '',
        regex=True,
    )
    expected_chromosomes = mutation_ids.map(chromosome_by_segment_id)
    chromosome_mismatches = matched_mutation_ids & (
        mutation_chromosomes != expected_chromosomes
    )
    if chromosome_mismatches.any():
        mismatched_segment_ids = mutation_ids[
            chromosome_mismatches
        ].unique()
        raise ValueError(
            "Chromosome does not match the copy number table for Segment_ID "
            "value(s): " + ", ".join(mismatched_segment_ids)
        )


def drop_unmatched_segment_id_mutations(copy_number_table, mutation_table):
    """Drop mutations whose supplied segment ID has no copy-number match."""
    copy_number_ids = copy_number_table['Segment_ID'].astype('string')
    mutation_ids = mutation_table['Segment_ID'].astype('string')
    matched_mutations = mutation_ids.isin(copy_number_ids)
    unmatched_snv_count = int((~matched_mutations).sum())

    if unmatched_snv_count:
        snv_label = 'SNV' if unmatched_snv_count == 1 else 'SNVs'
        warnings.warn(
            f'Dropping {unmatched_snv_count} unmatched {snv_label} because '
            'their Segment_ID is not present in the copy number table.',
            UserWarning,
            stacklevel=2,
        )

    return mutation_table.loc[matched_mutations].copy()


def load_input_tables(
    copy_number_path,
    mutation_table_path,
    drop_unmatched_snvs=False,
    sex=None,
    autosome_count=22,
    drop_unmatched_chromosomes=False,
):
    """Load paired copy number and mutation tables using their joint schema."""
    validate_sex_karyotype(sex)
    copy_number_header = pd.read_csv(copy_number_path, sep='\t', nrows=0)
    mutation_header = pd.read_csv(mutation_table_path, sep='\t', nrows=0)
    use_supplied_segment_ids = validate_input_table_headers(
        copy_number_header.columns,
        mutation_header.columns,
    )

    copy_number_dtypes = {'Chromosome': str}
    mutation_dtypes = {'Chromosome': str}
    copy_number_converters = {}
    mutation_converters = {}
    if use_supplied_segment_ids:
        copy_number_converters['Segment_ID'] = str
        mutation_converters['Segment_ID'] = str
    for column in (
        'Segment_Start',
        'Segment_End',
        'Major_CN',
        'Minor_CN',
    ):
        copy_number_converters[column] = str
    if 'Position' in mutation_header.columns:
        mutation_converters['Position'] = str
    if 'Mutation_ID' in mutation_header.columns:
        mutation_converters['Mutation_ID'] = str
    for column in ('Tumor_Ref_Count', 'Tumor_Alt_Count'):
        mutation_converters[column] = str

    copy_number_table = pd.read_csv(
        copy_number_path,
        sep='\t',
        dtype=copy_number_dtypes,
        converters=copy_number_converters,
    )
    mutation_table = pd.read_csv(
        mutation_table_path,
        sep='\t',
        dtype=mutation_dtypes,
        converters=mutation_converters,
        usecols=lambda column: (
            column not in MUTATION_COPY_NUMBER_ANNOTATION_COLUMNS
        ),
    )
    copy_number_table, mutation_table, _ = validate_or_drop_input_chromosomes(
        copy_number_table,
        mutation_table,
        autosome_count,
        sex,
        drop_unmatched_chromosomes=drop_unmatched_chromosomes,
    )
    copy_number_table = validate_segment_coordinates(copy_number_table)
    copy_number_table = validate_non_overlapping_segments(copy_number_table)
    copy_number_table = validate_copy_number_values(copy_number_table)
    mutation_table = validate_mutation_read_counts(mutation_table)
    if 'Position' in mutation_table.columns:
        mutation_table['Position'] = parse_positions(
            mutation_table['Position']
        )
    gritic_id_components = get_gritic_mutation_id_components(mutation_table)
    if use_supplied_segment_ids:
        validate_supplied_segment_ids(
            copy_number_table,
            mutation_table,
            allow_unmatched=drop_unmatched_snvs,
        )
        validate_source_scoped_mutation_components(
            mutation_table['Segment_ID'],
            gritic_id_components,
        )
        if drop_unmatched_snvs:
            mutation_table = drop_unmatched_segment_id_mutations(
                copy_number_table,
                mutation_table,
            )
    return copy_number_table, mutation_table

def calculate_ploidy(
    cn_table,
    autosome_count=22,
    sex=None,
    drop_unmatched_chromosomes=False,
):
    if sex is None:
        sex = infer_sex_from_copy_number_table(cn_table)
    else:
        validate_sex_karyotype(sex)
    cn_table = validate_or_drop_chromosome_rows(
        cn_table,
        autosome_count,
        sex,
        table_name='Copy number table',
        drop_unmatched_chromosomes=drop_unmatched_chromosomes,
    )
    if cn_table.empty:
        raise ValueError(
            'No copy-number segments remain after chromosome filtering'
        )
    cn_table = validate_segment_coordinates(cn_table)
    cn_table = validate_non_overlapping_segments(cn_table)
    cn_table = validate_copy_number_values(cn_table)
    cn_width = get_segment_widths(cn_table).astype(float)
    cn_total = cn_table['Major_CN'] + cn_table['Minor_CN']
    ploidy = np.average(cn_total,weights=cn_width)
    return ploidy


def infer_sex_from_copy_number_table(cn_table):
    """Infer an X/Y or Z/W karyotype from exact chromosome labels.

    A Y or W segment identifies the corresponding heterogametic karyotype.
    With only a homogametic chromosome present, X identifies XX and Z
    identifies ZZ. When no sex chromosome is represented, XX remains the
    default; callers can override this explicitly.
    """
    chromosomes = set(
        cn_table['Chromosome']
        .astype(str)
        .str.replace(r'^chr', '', regex=True)
    )
    xy_chromosomes = chromosomes & {'X', 'Y'}
    zw_chromosomes = chromosomes & {'Z', 'W'}
    if xy_chromosomes and zw_chromosomes:
        raise ValueError(
            'Cannot infer sample sex because the copy number table mixes '
            'X/Y and Z/W chromosome systems'
        )
    if 'Y' in chromosomes:
        return 'XY'
    if 'W' in chromosomes:
        return 'ZW'
    if 'Z' in chromosomes:
        return 'ZZ'
    return 'XX'


def calculate_normal_ploidy(
    cn_table,
    sex=None,
    autosome_count=22,
    drop_unmatched_chromosomes=False,
):
    if sex is None:
        sex = infer_sex_from_copy_number_table(cn_table)
    else:
        validate_sex_karyotype(sex)
    normal_table = validate_or_drop_chromosome_rows(
        cn_table,
        autosome_count,
        sex,
        table_name='Copy number table',
        drop_unmatched_chromosomes=drop_unmatched_chromosomes,
    )
    if normal_table.empty:
        raise ValueError(
            'No copy-number segments remain after chromosome filtering'
        )
    normal_table = validate_segment_coordinates(normal_table)
    normal_table = validate_non_overlapping_segments(normal_table)
    chromosome_width = get_segment_widths(normal_table).astype(
        float
    )
    normal_copy_number = normal_table['Chromosome'].map(
        lambda chromosome: get_normal_total_copy_number(chromosome, sex)
    )
    return np.average(normal_copy_number, weights=chromosome_width)


def calculate_nrpcc(
    cn_table,
    mutation_table,
    purity,
    sex=None,
    autosome_count=22,
    drop_unmatched_chromosomes=False,
):
    purity = validate_purity(purity)
    cn_table, mutation_table, sex = validate_or_drop_input_chromosomes(
        cn_table,
        mutation_table,
        autosome_count,
        sex,
        drop_unmatched_chromosomes=drop_unmatched_chromosomes,
    )
    tumor_ploidy = calculate_ploidy(
        cn_table,
        autosome_count=autosome_count,
        sex=sex,
    )
    normal_ploidy = calculate_normal_ploidy(
        cn_table,
        sex=sex,
        autosome_count=autosome_count,
    )
    sample_ploidy = purity*tumor_ploidy+normal_ploidy*(1-purity)
    if mutation_table.empty:
        raise ValueError(
            'Cannot calculate NRPCC because there are no mutations on the '
            'configured chromosomes'
        )
    mutation_table = validate_mutation_read_counts(mutation_table)
    coverage = (mutation_table['Tumor_Ref_Count']+mutation_table['Tumor_Alt_Count']).mean()
    nrpcc = purity*coverage/sample_ploidy
    return nrpcc,coverage,tumor_ploidy
def generate_segment_ids(cn_table):
    """Replace segment IDs with IDs generated from genomic coordinates."""
    cn_table = cn_table.copy()
    cn_table.loc[:, 'Segment_ID'] = (
        cn_table['Chromosome'].astype(str)
        + "-"
        + cn_table['Segment_Start'].astype(str)
        + "-"
        + cn_table['Segment_End'].astype(str)
    )
    return cn_table


def sort_copy_number_segments(copy_number_table):
    """Return copy-number rows in deterministic genomic coordinate order."""
    copy_number_table = copy_number_table.copy().reset_index(drop=True)
    chromosome_labels = (
        copy_number_table['Chromosome']
        .astype(str)
        .str.replace(r'^chr', '', regex=True)
    )
    numbered_chromosomes = chromosome_labels.str.fullmatch(r'[1-9]\d*')
    sex_chromosome_rank = chromosome_labels.map({
        'X': 0,
        'Y': 1,
        'Z': 2,
        'W': 3,
    })
    chromosome_group = pd.Series(
        np.select(
            [numbered_chromosomes, sex_chromosome_rank.notna()],
            [0, 1],
            default=2,
        ),
        index=copy_number_table.index,
    )
    numbered_chromosome_rank = pd.to_numeric(
        chromosome_labels.where(numbered_chromosomes),
        errors='coerce',
    )
    sort_keys = pd.DataFrame({
        'Chromosome_Group': chromosome_group,
        'Numbered_Chromosome_Rank': numbered_chromosome_rank,
        'Sex_Chromosome_Rank': sex_chromosome_rank,
        'Chromosome_Label': chromosome_labels,
        'Segment_Start': copy_number_table['Segment_Start'],
        'Segment_End': copy_number_table['Segment_End'],
    })
    ordered_indices = sort_keys.sort_values(
        [
            'Chromosome_Group',
            'Numbered_Chromosome_Rank',
            'Sex_Chromosome_Rank',
            'Chromosome_Label',
            'Segment_Start',
            'Segment_End',
        ],
        kind='mergesort',
        na_position='last',
    ).index
    return copy_number_table.iloc[ordered_indices].reset_index(drop=True)


def validate_non_overlapping_segments(copy_number_table):
    """Sort half-open copy-number intervals and reject any overlap."""
    copy_number_table = sort_copy_number_segments(copy_number_table)
    overlapping_intervals = []
    for chromosome, chromosome_table in copy_number_table.groupby(
        'Chromosome',
        sort=False,
    ):
        active_start = None
        active_end = None
        for segment in chromosome_table.itertuples(index=False):
            segment_start = int(segment.Segment_Start)
            segment_end = int(segment.Segment_End)
            if active_end is not None and segment_start < active_end:
                overlapping_intervals.append(
                    f'{chromosome}:{active_start}-{active_end} and '
                    f'{chromosome}:{segment_start}-{segment_end}'
                )
            if active_end is None or segment_end > active_end:
                active_start = segment_start
                active_end = segment_end

    if overlapping_intervals:
        raise ValueError(
            'Copy number segments must not overlap within a chromosome; '
            'overlapping interval(s): '
            + ', '.join(overlapping_intervals)
        )
    return copy_number_table


def merge_segments(
    cn_table,
    return_segment_id_map=False,
    *,
    _validated=False,
):
    """Merge equal-CN segments that share a half-open interval boundary."""

    if _validated:
        cn_table = cn_table.copy()
    else:
        cn_table = validate_segment_coordinates(cn_table)
        cn_table = validate_non_overlapping_segments(cn_table)
    cn_table.loc[:,'Gain_Type'] = cn_table['Major_CN'].astype(str) + "_"+cn_table['Minor_CN'].astype(str)
    source_ids_by_index = None
    if return_segment_id_map:
        source_ids_by_index = {
            index: [str(segment_id)]
            for index, segment_id in cn_table['Segment_ID'].items()
        }
    indexes_to_delete = []
    for _, chr_data in cn_table.groupby("Chromosome"):
        for i in range(len(chr_data.index)-1):
            index = chr_data.index[i]
            forward_index= chr_data.index[i+1]

            equal_copy_number = (
                chr_data.loc[index, 'Gain_Type']
                == chr_data.loc[forward_index, 'Gain_Type']
            )
            coordinate_adjacent = (
                int(chr_data.loc[forward_index, 'Segment_Start'])
                == int(chr_data.loc[index, 'Segment_End'])
            )
            if equal_copy_number and coordinate_adjacent:
                
                indexes_to_delete.append(index)
                cn_table.loc[forward_index,'Segment_Start'] = cn_table.loc[index,'Segment_Start']
                if return_segment_id_map:
                    source_ids_by_index[forward_index] = (
                        source_ids_by_index[index]
                        + source_ids_by_index[forward_index]
                    )

    merged_table = cn_table.loc[~cn_table.index.isin(indexes_to_delete)].copy()
    if not _validated:
        merged_table = validate_segment_coordinates(merged_table)
    merged_table = generate_segment_ids(merged_table)
    if not return_segment_id_map:
        return merged_table.reset_index(drop=True)

    segment_id_map = {}
    for index, merged_segment_id in merged_table['Segment_ID'].items():
        for source_segment_id in source_ids_by_index[index]:
            segment_id_map[source_segment_id] = merged_segment_id
    return merged_table.reset_index(drop=True), segment_id_map
def assign_cn_to_snv(
    snv_table,
    cn_table,
    use_supplied_segment_ids,
    drop_unmatched_snvs=False,
    *,
    _validated=False,
):
    """Attach copy number annotations using supplied IDs or genomic position."""
    if _validated:
        cn_table = cn_table.copy()
    else:
        cn_table = validate_segment_coordinates(cn_table)
        cn_table = validate_non_overlapping_segments(cn_table)
    snv_table = snv_table.copy()

    if use_supplied_segment_ids:
        cn_table.loc[:, 'Segment_ID'] = cn_table['Segment_ID'].astype(str)
        snv_table.loc[:, 'Segment_ID'] = snv_table['Segment_ID'].astype(str)
    else:
        cn_table = generate_segment_ids(cn_table)
        snv_table.loc[:, 'Segment_ID'] = "None"
        numeric_positions = (
            snv_table['Position']
            if _validated
            else parse_positions(snv_table['Position'])
        )
        for _, segment in cn_table.iterrows():
            matching_segments = (
                (snv_table['Chromosome'] == segment['Chromosome'])
                & (numeric_positions >= segment['Segment_Start'])
                & (numeric_positions < segment['Segment_End'])
            )
            snv_table.loc[matching_segments, 'Segment_ID'] = segment['Segment_ID']

        matched_mutations = snv_table['Segment_ID'] != 'None'
        unmatched_snv_count = int((~matched_mutations).sum())
        if unmatched_snv_count:
            snv_label = 'SNV' if unmatched_snv_count == 1 else 'SNVs'
            message = (
                f'{unmatched_snv_count} unmatched {snv_label} whose '
                'Chromosome and Position do not match any copy number segment'
            )
            if not drop_unmatched_snvs:
                raise ValueError(
                    'Mutation table contains ' + message + '.'
                )
            warnings.warn(
                'Dropping ' + message + '.',
                UserWarning,
                stacklevel=2,
            )
            snv_table = snv_table.loc[matched_mutations].copy()

    copy_number_annotation_columns = [
        'Chromosome',
        'Segment_ID',
        'Segment_Start',
        'Segment_End',
        'Major_CN',
        'Minor_CN',
    ]
    existing_annotation_columns = [
        column for column in copy_number_annotation_columns[2:]
        if column in snv_table.columns
    ]
    snv_table = snv_table.drop(columns=existing_annotation_columns)
    snv_table = snv_table.merge(
        cn_table[copy_number_annotation_columns],
        on=['Chromosome', 'Segment_ID'],
        how='inner',
        validate='many_to_one',
    )
    snv_table.loc[:,'Total_CN'] = snv_table['Major_CN'] + snv_table['Minor_CN']
    return snv_table

def get_valid_subclones(
    subclone_table,
    min_ccf=0.01,
    max_ccf=0.9,
    min_fraction=0.1,
):
    retained_ccf = (
        (subclone_table['Subclone_CCF'] >= min_ccf)
        & (subclone_table['Subclone_CCF'] <= max_ccf)
    )
    subclone_table = subclone_table[retained_ccf].copy()
    subclone_fraction_norm = subclone_table['Subclone_Fraction']/subclone_table['Subclone_Fraction'].sum()
    subclone_table = subclone_table[subclone_fraction_norm > min_fraction]
    return subclone_table.reset_index(drop=True).copy()
def filter_excess_subclones(subclone_table):
    if subclone_table is None or len(subclone_table.index)==1:
        return subclone_table
    subclone_table = subclone_table.sort_values(by=['Subclone_CCF'],ascending=False).copy()

    top_clone = subclone_table.iloc[0:1]
    
    other_clones = subclone_table.iloc[1:]
    
    combined_clone_ccf = np.average(other_clones['Subclone_CCF'],weights=other_clones['Subclone_Fraction'])
    combined_clone_fraction = np.sum(other_clones['Subclone_Fraction'])
    
    combined_clone_table =pd.DataFrame({'Subclone_CCF':combined_clone_ccf,'Subclone_Fraction':combined_clone_fraction,'Cluster':other_clones['Cluster'].iloc[0]},index=[other_clones.index[0]])
    new_subclone_table = pd.concat([top_clone,combined_clone_table])
    return new_subclone_table
