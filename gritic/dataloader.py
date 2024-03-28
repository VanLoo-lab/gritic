import numpy as np
import pandas as pd

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
def merge_segments(cn_table):
    
    cn_table = cn_table.copy().reset_index(drop=True)
    cn_table.loc[:,'Gain_Type'] = cn_table['Major_CN'].astype(str) + "_"+cn_table['Minor_CN'].astype(str)
    indexes_to_delete = []
    for chromosome,chr_data in cn_table.groupby("Chromosome"):
        #print(chr_data)
        for i in range(len(chr_data.index)-1):
            index = chr_data.index[i]
            forward_index= chr_data.index[i+1]

            if chr_data.loc[index,'Gain_Type'] == chr_data.loc[forward_index,'Gain_Type']:
                
                indexes_to_delete.append(index)
                cn_table.loc[forward_index,'Segment_Start'] = cn_table.loc[index,'Segment_Start']

    merged_table = cn_table.loc[~cn_table.index.isin(indexes_to_delete)].copy()
    return merged_table.reset_index(drop=True)
def filter_sex_chromosomes(mutation_table):
    autosomes = list(map(str,range(1,23)))
    mutation_table = mutation_table[mutation_table['Chromosome'].isin(autosomes)]
    return mutation_table
def assign_cn_to_snv(snv_table,cn_table):
    cn_table = cn_table.copy()
    snv_table = snv_table.copy()
    cn_table.loc[:,'Segment_ID'] = cn_table['Chromosome'].astype(str) + "-" + cn_table['Segment_Start'].astype(str) + "-" + cn_table['Segment_End'].astype(str)
    snv_table.loc[:,'Segment_ID'] = "None"
    for index,segment in cn_table.iterrows():
        matching_segments = (snv_table['Chromosome']==segment['Chromosome']) & (snv_table['Position'] >= segment['Segment_Start']) & (snv_table['Position'] <= segment['Segment_End'])
        snv_table.loc[matching_segments,'Segment_ID']=segment['Segment_ID']
    snv_table = snv_table.merge(cn_table,how='inner')
    snv_table.loc[:,'Total_CN'] = snv_table['Major_CN'] + snv_table['Minor_CN']
    #snv_table = snv_table[['Chromosome','Arm','Position','Tumor_Ref_Count','Tumor_Alt_Count','Segment_Start','Segment_End','Total_CN']]
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


