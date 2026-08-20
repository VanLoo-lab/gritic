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
    if not use_supplied_segment_ids and 'Position' not in mutation_columns:
        raise ValueError(
            "Mutation table must contain Position unless Segment_ID is present "
            "in both the copy number and mutation tables"
        )

    return use_supplied_segment_ids


def validate_supplied_segment_ids(copy_number_table, mutation_table):
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

    unknown_mutation_ids = mutation_ids[
        ~mutation_ids.isin(copy_number_ids)
    ].unique()
    if unknown_mutation_ids.size:
        raise ValueError(
            "Mutation table Segment_ID value(s) not present in the copy "
            "number table: " + ", ".join(unknown_mutation_ids)
        )

    chromosome_by_segment_id = pd.Series(
        copy_number_table['Chromosome'].astype(str).str.replace(
            'chr',
            '',
            regex=False,
        ).to_numpy(),
        index=copy_number_ids,
    )
    mutation_chromosomes = mutation_table['Chromosome'].astype(str).str.replace(
        'chr',
        '',
        regex=False,
    )
    expected_chromosomes = mutation_ids.map(chromosome_by_segment_id)
    chromosome_mismatches = mutation_chromosomes != expected_chromosomes
    if chromosome_mismatches.any():
        mismatched_segment_ids = mutation_ids[
            chromosome_mismatches
        ].unique()
        raise ValueError(
            "Chromosome does not match the copy number table for Segment_ID "
            "value(s): " + ", ".join(mismatched_segment_ids)
        )


def load_input_tables(copy_number_path, mutation_table_path):
    """Load paired copy number and mutation tables using their joint schema."""
    copy_number_header = pd.read_csv(copy_number_path, sep='\t', nrows=0)
    mutation_header = pd.read_csv(mutation_table_path, sep='\t', nrows=0)
    use_supplied_segment_ids = validate_input_table_headers(
        copy_number_header.columns,
        mutation_header.columns,
    )

    copy_number_dtypes = {
        'Chromosome': str,
        'Segment_Start': int,
        'Segment_End': int,
        'Major_CN': int,
        'Minor_CN': int,
    }
    mutation_dtypes = {
        'Chromosome': str,
        'Tumor_Ref_Count': int,
        'Tumor_Alt_Count': int,
    }
    copy_number_converters = {}
    mutation_converters = {}
    if use_supplied_segment_ids:
        copy_number_converters['Segment_ID'] = str
        mutation_converters['Segment_ID'] = str
    else:
        mutation_dtypes['Position'] = int

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
    )
    if use_supplied_segment_ids:
        validate_supplied_segment_ids(copy_number_table, mutation_table)
    return copy_number_table, mutation_table

def calculate_ploidy(cn_table):
    cn_table = cn_table[cn_table['Chromosome'].isin(list(map(str,range(1,23)))+['X','Y'])].copy()
    cn_width = cn_table['Segment_End']-cn_table['Segment_Start']
    cn_total = cn_table['Major_CN'] + cn_table['Minor_CN']
    ploidy = np.average(cn_total,weights=cn_width)
    return ploidy
def calculate_normal_ploidy(sex):
    sex_table = pd.read_csv('/camp/home/bakert/secure-working/GRITIC/Analysis_Pipeline/resources/chrom_arm_positions.tsv',sep='\t').groupby('Chromosome').agg({'Arm_End':'max'}).reset_index()
    sex_table = sex_table.rename(columns={'Arm_End':'Chromosome_Size'})
    
    if sex is None or sex == 'XX':
        x_ploidy = 2
        sex_table = sex_table[sex_table['Chromosome']!='Y']
    elif sex == 'XY':
        x_ploidy = 1
    else:
        raise ValueError(f'invalid sex {sex}')
    
    sex_table['Total_Copy_Number'] = np.where(sex_table['Chromosome'] == 'X',x_ploidy,2)
    sex_table['Total_Copy_Number'] = np.where(sex_table['Chromosome'] == 'Y',1,sex_table['Total_Copy_Number'])

    return np.average(sex_table['Total_Copy_Number'],weights=sex_table['Chromosome_Size'])
def calculate_nrpcc(cn_table,mutation_table,purity,sex=None):
    tumor_ploidy = calculate_ploidy(cn_table)
    normal_ploidy = calculate_normal_ploidy(sex)
    sample_ploidy = purity*tumor_ploidy+normal_ploidy*(1-purity)
    mutation_table = mutation_table[mutation_table['Chromosome'].isin(list(map(str,range(1,23)))+['X','Y'])].copy()
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


def merge_segments(cn_table, return_segment_id_map=False):
    """Merge equal-CN segments, optionally mapping supplied IDs to merged IDs."""
    
    cn_table = cn_table.copy().reset_index(drop=True)
    cn_table.loc[:,'Gain_Type'] = cn_table['Major_CN'].astype(str) + "_"+cn_table['Minor_CN'].astype(str)
    source_ids_by_index = None
    if return_segment_id_map:
        source_ids_by_index = {
            index: [str(segment_id)]
            for index, segment_id in cn_table['Segment_ID'].items()
        }
    indexes_to_delete = []
    for chromosome,chr_data in cn_table.groupby("Chromosome"):
        #print(chr_data)
        for i in range(len(chr_data.index)-1):
            index = chr_data.index[i]
            forward_index= chr_data.index[i+1]

            if chr_data.loc[index,'Gain_Type'] == chr_data.loc[forward_index,'Gain_Type']:
                
                indexes_to_delete.append(index)
                cn_table.loc[forward_index,'Segment_Start'] = cn_table.loc[index,'Segment_Start']
                if return_segment_id_map:
                    source_ids_by_index[forward_index] = (
                        source_ids_by_index[index]
                        + source_ids_by_index[forward_index]
                    )

    merged_table = cn_table.loc[~cn_table.index.isin(indexes_to_delete)].copy()
    merged_table = generate_segment_ids(merged_table)
    if not return_segment_id_map:
        return merged_table.reset_index(drop=True)

    segment_id_map = {}
    for index, merged_segment_id in merged_table['Segment_ID'].items():
        for source_segment_id in source_ids_by_index[index]:
            segment_id_map[source_segment_id] = merged_segment_id
    return merged_table.reset_index(drop=True), segment_id_map
def filter_sex_chromosomes(mutation_table):
    autosomes = list(map(str,range(1,23)))
    mutation_table = mutation_table[mutation_table['Chromosome'].isin(autosomes)]
    return mutation_table
def assign_cn_to_snv(snv_table, cn_table, use_supplied_segment_ids):
    """Attach copy number annotations using supplied IDs or genomic position."""
    cn_table = cn_table.copy()
    snv_table = snv_table.copy()

    if use_supplied_segment_ids:
        cn_table.loc[:, 'Segment_ID'] = cn_table['Segment_ID'].astype(str)
        snv_table.loc[:, 'Segment_ID'] = snv_table['Segment_ID'].astype(str)
    else:
        cn_table = generate_segment_ids(cn_table)
        snv_table.loc[:, 'Segment_ID'] = "None"
        for _, segment in cn_table.iterrows():
            matching_segments = (
                (snv_table['Chromosome'] == segment['Chromosome'])
                & (snv_table['Position'] >= segment['Segment_Start'])
                & (snv_table['Position'] <= segment['Segment_End'])
            )
            snv_table.loc[matching_segments, 'Segment_ID'] = segment['Segment_ID']

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

def get_major_cn_mode_from_cn_table(cn_table):
    cn_table = cn_table.copy()
    cn_table.loc[:,'Segment_Width'] = cn_table['Segment_End']-cn_table['Segment_Start']
    major_cn_widths = cn_table.groupby('Major_CN').agg({'Segment_Width':'sum'})
    max_width = major_cn_widths['Segment_Width'].max()

    major_cn_mode = major_cn_widths[major_cn_widths['Segment_Width']==max_width].index[0]
    return major_cn_mode
    
def get_major_cn_mode(sample):
    return get_major_cn_mode_from_cn_table(sample.cn_table)
def get_valid_subclones(subclone_table,max_ccf=0.9,min_fraction=0.1):
    subclone_table = subclone_table[subclone_table['Subclone_CCF'] <= max_ccf].copy()
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
