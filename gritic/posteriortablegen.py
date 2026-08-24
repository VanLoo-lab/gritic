import os
import bz2
import pickle

import numpy as np
import pandas as pd
import warnings


NON_PARSIMONY_PENALTY_COEFFICIENT = 2.7
ROUTE_KEY_COLUMNS = ['Sample_ID', 'Segment_ID', 'Route']


def apply_non_parsimony_penalty(sample_table):
    """Return a copy with route probabilities penalized and renormalized."""
    required_columns = ROUTE_KEY_COLUMNS + [
        'Average_N_Events',
        'Probability',
    ]
    missing_columns = [
        column for column in required_columns
        if column not in sample_table.columns
    ]
    if missing_columns:
        raise ValueError(
            'Cannot apply the non-parsimony penalty because the gain timing '
            f"table is missing columns: {', '.join(missing_columns)}"
        )

    route_table = sample_table[required_columns].drop_duplicates()
    inconsistent_routes = route_table.duplicated(
        subset=ROUTE_KEY_COLUMNS,
        keep=False,
    )
    if inconsistent_routes.any():
        inconsistent_keys = (
            route_table.loc[inconsistent_routes, ROUTE_KEY_COLUMNS]
            .drop_duplicates()
            .astype(str)
            .agg(':'.join, axis=1)
            .tolist()
        )
        raise ValueError(
            'Each route must have exactly one Probability and '
            'Average_N_Events value; inconsistent routes: '
            f"{', '.join(inconsistent_keys)}"
        )

    route_table = route_table.copy()
    if route_table[ROUTE_KEY_COLUMNS].isnull().any().any():
        raise ValueError(
            'Sample_ID, Segment_ID, and Route must be present for every '
            'gain timing table row'
        )
    route_table['Probability'] = pd.to_numeric(
        route_table['Probability'],
        errors='coerce',
    )
    route_table['Average_N_Events'] = pd.to_numeric(
        route_table['Average_N_Events'],
        errors='coerce',
    )
    if (route_table['Probability'] < 0).any():
        raise ValueError('Route probabilities must be non-negative')
    if (route_table['Average_N_Events'] < 0).any():
        raise ValueError('Average_N_Events values must be non-negative')

    route_table['_Route_Data_Valid'] = (
        np.isfinite(route_table['Probability'])
        & np.isfinite(route_table['Average_N_Events'])
    )
    group_columns = ['Sample_ID', 'Segment_ID']
    route_table['_Segment_Data_Valid'] = route_table.groupby(
        group_columns,
        sort=False,
    )['_Route_Data_Valid'].transform('all')
    valid_routes = route_table['_Segment_Data_Valid']
    route_table['_Penalized_Weight'] = np.nan
    route_table.loc[valid_routes, '_Penalized_Weight'] = (
        route_table.loc[valid_routes, 'Probability']
        * np.exp(
            -NON_PARSIMONY_PENALTY_COEFFICIENT
            * route_table.loc[valid_routes, 'Average_N_Events']
        )
    )
    normalizers = route_table.groupby(
        group_columns,
        sort=False,
    )['_Penalized_Weight'].transform('sum')
    valid_normalizers = (
        route_table['_Segment_Data_Valid']
        & np.isfinite(normalizers)
        & (normalizers > 0)
    )
    route_table.loc[valid_normalizers, 'Probability'] = (
        route_table.loc[valid_normalizers, '_Penalized_Weight']
        / normalizers[valid_normalizers]
    )
    route_table.loc[~valid_normalizers, 'Probability'] = np.nan

    probability_lookup = route_table.set_index(
        ROUTE_KEY_COLUMNS
    )['Probability']
    row_route_keys = pd.MultiIndex.from_frame(
        sample_table[ROUTE_KEY_COLUMNS]
    )
    penalized_probabilities = probability_lookup.reindex(
        row_route_keys
    ).to_numpy()
    penalized_table = sample_table.copy()
    penalized_table['Probability'] = penalized_probabilities
    return penalized_table


def load_route_table(path):
    read_cols = ['Sample_ID','Segment_ID','Route','Average_N_Events','Average_Pre_WGD_Losses','Average_Post_WGD_Losses','Probability','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','WGD_Status','N_Mutations']
    route_table = pd.read_csv(path,sep="\t",usecols=read_cols,dtype={'Chromosome':str})
    baseline_cn = np.where(route_table['WGD_Status'],2.1,1.1)
    route_table = route_table[route_table['Major_CN']>baseline_cn]
    route_table = route_table[route_table['N_Mutations']>=20]
    route_table = route_table.drop_duplicates()
    return route_table
def shorten_dict_routes(timing_dict):
    new_dict = {}
    for route,route_dict in timing_dict.items():
        new_dict[route[0:9]] = route_dict
    return new_dict
def get_random_timing(route_dict):
    n_samples = route_dict['WGD'].shape[0]
    sample_index = np.random.randint(0,n_samples)
    wgd_timing = route_dict['WGD'][sample_index]
    #gain_timing = np.sort(np.array([route_dict[node][sample_index] for node in route_dict.keys() if node!='WGD']))
    nodes = [node for node in route_dict.keys() if node!='WGD']

    gain_timing = np.array([route_dict[node][sample_index] for node in nodes])
    gain_timing_indices = np.argsort(gain_timing)
    gain_timing = gain_timing[gain_timing_indices]
    nodes = np.array(nodes)[gain_timing_indices]

    return gain_timing,wgd_timing,nodes
def produce_timing_segment_table(segment_table,timing_dict,segment_id,n_samples=100):
    segment_data= {'Segment_ID':segment_id,'Node':[],'Gain_Timing':[],'WGD_Timing':[],'Posterior_Sample_Index':[],'Route':[]}
    if segment_table['Probability'].isnull().any():
        return None
    for i in range(n_samples):
        #occasionally we get some routes with probability very close to zero, so we just skip them
        while True:
            probs = segment_table.Probability/np.sum(segment_table.Probability)
            route = np.random.choice(segment_table.Route,p=probs)
            if route in timing_dict.keys():
                break
    
        gain_timing,wgd_timing,nodes = get_random_timing(timing_dict[route]['Timing'])
        segment_data['Gain_Timing'].extend(gain_timing)
        segment_data['WGD_Timing'].extend([wgd_timing]*gain_timing.size)
        segment_data['Posterior_Sample_Index'].extend([i]*gain_timing.size)
        segment_data['Route'].extend([route]*gain_timing.size)
        segment_data['Node'].extend(nodes)
    return pd.DataFrame(segment_data)

def load_timing_from_dict(segment_path):
    input_file = bz2.BZ2File(segment_path,'rb')
    timing_dict = pickle.load(input_file)
    input_file.close()
    return timing_dict
def get_sample_posterior_table(
    sample_table_path,
    input_dir,
    sample_id: str,
    n_posterior_samples: int = 100,
    apply_penalty: bool = False,
):

    if not isinstance(apply_penalty, (bool, np.bool_)):
        raise ValueError('apply_penalty must be a boolean')

    sample_table = pd.read_csv(
        sample_table_path,
        sep='\t',
        dtype={'Chromosome': str},
        converters={
            'Sample_ID': str,
            'Segment_ID': str,
            'Route': str,
        },
    )
    if apply_penalty:
        sample_table = apply_non_parsimony_penalty(sample_table)
    full_segment_table = sample_table[['Segment_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','N_Mutations']].drop_duplicates()
    node_table = sample_table[['Route','Node','Node_Phasing','Major_CN','Minor_CN','WGD_Status']].drop_duplicates()

    dict_dir = f'{input_dir}/{sample_id}_timing_dicts'
    segment_frames = []
    for segment_file in os.listdir(dict_dir):
        
        if not segment_file.endswith("_timing_dict.bz2"):
            continue
        
        segment_path = f'{dict_dir}/{segment_file}'
        timing_dict = load_timing_from_dict(segment_path)
        timing_dict = shorten_dict_routes(timing_dict)
        
        segment_id = segment_file.removesuffix('_timing_dict.bz2')
        
        segment_table = sample_table[sample_table['Segment_ID']==segment_id]

        if len(segment_table)==0:
            continue
        segment_frame = produce_timing_segment_table(segment_table,timing_dict,segment_id,n_samples=n_posterior_samples)
        if segment_frame is None:
            warnings.warn(f'WARNING {segment_id} has NAs in probability, skipping in posterior.')
            continue
        segment_frames.append(segment_frame)
    if not segment_frames:
        return pd.DataFrame()
    segment_frame = pd.concat(segment_frames)
    segment_frame = pd.merge(full_segment_table,segment_frame,on=['Segment_ID'],how='inner')
    segment_frame = pd.merge(segment_frame,node_table,how='inner')
    segment_frame = segment_frame.sort_values(by=['Segment_ID','Posterior_Sample_Index','Gain_Timing'])
    segment_frame['Gain_Index'] = segment_frame.groupby(['Segment_ID','Posterior_Sample_Index']).cumcount()
    segment_frame = segment_frame.sort_values(by=['Segment_ID','Posterior_Sample_Index','Gain_Index'])
    return segment_frame

def get_segment_posterior_table_summary(segment_posterior_table):
    n_samples = segment_posterior_table['Posterior_Sample_Index'].max()+1

    segment_posterior_summary = {'Gain_Index':[],'Proportion':[],'Timing_Median':[],'Timing_Low_CI':[],'Timing_High_CI':[],'Pre_WGD_Probability':[],'Post_WGD_Probability':[],'WGD_Timing_Median':[],'WGD_Timing_Low_CI':[],'WGD_Timing_High_CI':[]}
    for gain_index,gain_index_table  in segment_posterior_table.groupby('Gain_Index'):
        gain_index_proportion = len(gain_index_table)/n_samples
        gain_index_low_ci = np.percentile(gain_index_table['Gain_Timing'],2.5)
        gain_index_high_ci = np.percentile(gain_index_table['Gain_Timing'],97.5)
        gain_index_median = np.median(gain_index_table['Gain_Timing'])
        if np.isnan(segment_posterior_table['WGD_Timing']).any():
            wgd_timing_median = np.nan
            wgd_timing_low_ci = np.nan
            wgd_timing_high_ci = np.nan
            pre_wgd_probability = np.nan
        else:
            wgd_timing_median = np.median(gain_index_table['WGD_Timing'])
            wgd_timing_low_ci = np.percentile(gain_index_table['WGD_Timing'],2.5)
            wgd_timing_high_ci = np.percentile(gain_index_table['WGD_Timing'],97.5)
            pre_wgd_probability = np.sum(gain_index_table['WGD_Timing']>gain_index_table['Gain_Timing'])/len(gain_index_table)
        
        segment_posterior_summary['Gain_Index'].append(gain_index+1)
        segment_posterior_summary['Proportion'].append(gain_index_proportion)
        segment_posterior_summary['Timing_Median'].append(gain_index_median)
        segment_posterior_summary['Timing_Low_CI'].append(gain_index_low_ci)
        segment_posterior_summary['Timing_High_CI'].append(gain_index_high_ci)

        segment_posterior_summary['Pre_WGD_Probability'].append(pre_wgd_probability)
        segment_posterior_summary['Post_WGD_Probability'].append(1-pre_wgd_probability)

        segment_posterior_summary['WGD_Timing_Median'].append(wgd_timing_median)
        segment_posterior_summary['WGD_Timing_Low_CI'].append(wgd_timing_low_ci)
        segment_posterior_summary['WGD_Timing_High_CI'].append(wgd_timing_high_ci)
    segment_posterior_summary = pd.DataFrame(segment_posterior_summary)
    return segment_posterior_summary
def get_sample_posterior_table_summary(sample_posterior_table,min_proportion_threshold=0.8):
    segment_summary_data = sample_posterior_table[['Segment_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','N_Mutations']].drop_duplicates()
    sample_posterior_summary_store = []
    for segment_id,sample_segment_table in sample_posterior_table.groupby('Segment_ID'):
        segment_posterior_summary = get_segment_posterior_table_summary(sample_segment_table)
        segment_posterior_summary['Segment_ID'] = segment_id
        sample_posterior_summary_store.append(segment_posterior_summary)
    sample_posterior_summary = pd.concat(sample_posterior_summary_store)
    sample_posterior_summary = pd.merge(segment_summary_data,sample_posterior_summary,on=['Segment_ID'])
    sample_posterior_summary = sample_posterior_summary[sample_posterior_summary['Proportion']>=min_proportion_threshold]
    return sample_posterior_summary
