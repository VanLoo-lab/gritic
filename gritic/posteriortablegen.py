import os
import bz2
import pickle

import numpy as np
import pandas as pd



def load_route_table(path):
    read_cols = ['Sample_ID','Segment_ID','Route','Average_N_Events','Average_Pre_WGD_Losses','Average_Post_WGD_Losses','Probability','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','WGD_Status','N_Mutations']
    route_table = pd.read_csv(path,sep="\t",usecols=read_cols,dtype={'Chromosome':str})
    print(route_table)
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
def get_sample_posterior_table(sample_table_path,input_dir,sample_id):

    sample_table = pd.read_csv(sample_table_path,sep='\t',dtype={'Chromosome':str})
    full_segment_table = sample_table[['Segment_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN']].drop_duplicates()


    dict_dir = f'{input_dir}/{sample_id}_timing_dicts'
    segment_frames = []
    for segment_file in os.listdir(dict_dir):
        
        if not segment_file.endswith(".bz2"):
            continue
        
        segment_path = f'{dict_dir}/{segment_file}'
        timing_dict = load_timing_from_dict(segment_path)
        timing_dict = shorten_dict_routes(timing_dict)
        
        segment_id = segment_file.replace('_timing_dict.bz2','')
        
        segment_table = sample_table[sample_table['Segment_ID']==segment_id]

        if len(segment_table)==0:
            continue
        segment_frame = produce_timing_segment_table(segment_table,timing_dict,segment_id,n_samples=250)
        segment_frames.append(segment_frame)
    segment_frame = pd.concat(segment_frames)
    segment_frame = pd.merge(full_segment_table,segment_frame,on=['Segment_ID'])
    return segment_frame

def get_segment_posterior_table_summary(segment_posterior_table):
    segment_posterior_table['Gain_Index'] = segment_posterior_table.groupby(['Segment_ID','Posterior_Sample_Index']).cumcount()
    n_samples = segment_posterior_table['Posterior_Sample_Index'].max()+1

    segment_posterior_summary = {'Gain_Index':[],'Proportion':[],'Timing_Median':[],'Timing_Low_CI':[],'Timing_High_CI':[]}
    for gain_index,gain_index_table  in segment_posterior_table.groupby('Gain_Index'):
        gain_index_proportion = len(gain_index_table)/n_samples
        gain_index_low_ci = np.percentile(gain_index_table['Gain_Timing'],2.5)
        gain_index_high_ci = np.percentile(gain_index_table['Gain_Timing'],97.5)
        gain_index_median = np.median(gain_index_table['Gain_Timing'])

        segment_posterior_summary['Gain_Index'].append(gain_index+1)
        segment_posterior_summary['Proportion'].append(gain_index_proportion)
        segment_posterior_summary['Timing_Median'].append(gain_index_median)
        segment_posterior_summary['Timing_Low_CI'].append(gain_index_low_ci)
        segment_posterior_summary['Timing_High_CI'].append(gain_index_high_ci)
    segment_posterior_summary = pd.DataFrame(segment_posterior_summary)
    return segment_posterior_summary
def get_sample_posterior_table_summary(sample_posterior_table,min_proportion_threshold=0.5):
    segment_summary_data = sample_posterior_table[['Segment_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN']].drop_duplicates()
    sample_posterior_summary_store = []
    for segment_id,sample_segment_table in sample_posterior_table.groupby('Segment_ID'):
        segment_posterior_summary = get_segment_posterior_table_summary(sample_segment_table)
        segment_posterior_summary['Segment_ID'] = segment_id
        sample_posterior_summary_store.append(segment_posterior_summary)
    sample_posterior_summary = pd.concat(sample_posterior_summary_store)
    sample_posterior_summary = pd.merge(segment_summary_data,sample_posterior_summary,on=['Segment_ID'])
    sample_posterior_summary = sample_posterior_summary[sample_posterior_summary['Proportion']>=min_proportion_threshold]
    return sample_posterior_summary


'''if apply_prior:
    prior_corrected_probabilities = apply_prior_quick(route_table,2.7)
    route_table['Probability'] = prior_corrected_probabilities'''

