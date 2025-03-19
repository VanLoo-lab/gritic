import warnings
import numpy as np
import pandas as pd
from scipy.stats import binom

from numba import njit,prange

import time
import gritic.multiplicityoptimiser as multiplicityoptimiser
import gritic.dataloader as dataloader
def get_major_cn_mode_from_cn_table(cn_table):
    cn_table['Segment_Width'] = cn_table['Segment_End']-cn_table['Segment_Start']
    major_cn_widths = cn_table.groupby('Major_CN').agg({'Segment_Width':'sum'})
    max_width = major_cn_widths['Segment_Width'].max()

    major_cn_mode = major_cn_widths[major_cn_widths['Segment_Width']==max_width].index[0]
    return major_cn_mode

def get_major_cn_mode(sample):
    cn_table = sample.cn_table.copy()
    return get_major_cn_mode_from_cn_table(cn_table)

@njit(parallel=True)                           
def log_likelihood_numba_parallel(mult_array,mult_states):
    log_likelihood_store = np.zeros(mult_states.shape[0])                                
    for i in prange(mult_states.shape[0]):
        mult_state = mult_states[i]
        mult_state = np.clip(mult_state,0.0,1.0)
        log_likelihood = np.multiply(mult_state, mult_array)
    
        log_likelihood = np.sum(log_likelihood, axis=1)
        log_likelihood = np.sum(np.log(log_likelihood + 2.2e-300))
        log_likelihood_store[i] = log_likelihood
    return log_likelihood_store

@njit                           
def log_likelihood_numba(mult_array,mult_states):
    log_likelihood_store = np.zeros(mult_states.shape[0])                                  
    for i in range(mult_states.shape[0]):
        mult_state = mult_states[i]
        mult_state = np.clip(mult_state,0.0,1.0)
        log_likelihood = np.multiply(mult_state, mult_array)
    
        log_likelihood = np.sum(log_likelihood, axis=1)
        log_likelihood = np.sum(np.log(log_likelihood + 2.2e-300))
        log_likelihood_store[i] = log_likelihood
    return log_likelihood_store


@njit  
def evaluate_likelihood_array_numba(full_states,non_phased_array,reads_correction_array,tolerance=1e-8):
    
    full_states = np.multiply(full_states,reads_correction_array)
    full_states = full_states/np.sum(full_states,axis=1).reshape(-1,1)
    
    prob_sums = np.sum(non_phased_array,axis=1)
    probs_sums_good = (prob_sums>(1.0-tolerance)) & (prob_sums<(1.0+tolerance))

    non_phased_array = non_phased_array[probs_sums_good,:]
    if non_phased_array.shape[0] ==0:
        return np.nan*np.ones_like(full_states[:,0])
    likelihood_array = log_likelihood_numba_parallel(non_phased_array,full_states)

    return likelihood_array

'''@njit 
def evaluate_likelihood_array_one_subclone_numba(clonal_states,non_phased_array,reads_correction_array,n_subclone_points=25):
      ll_store = np.zeros(clonal_states.shape[0])
      for index in range(clonal_states.shape[0]):
        clonal_state = clonal_states[index]
        subclonal_share = np.linspace(0,1,n_subclone_points)
        
        mult_states = np.zeros((n_subclone_points,clonal_state.size+1))
        for i in range(mult_states.shape[0]):
          mult_states[i,0:clonal_state.size] = clonal_state*(1-subclonal_share[i])
          mult_states[i,-1] = subclonal_share[i]
    
       
        likelihood_array = evaluate_likelihood_array_numba(mult_states,non_phased_array,reads_correction_array)
        exponential_constant = np.max(likelihood_array)
        ll_sum = np.log(np.sum(np.exp(likelihood_array-exponential_constant)))+exponential_constant
        ll_store[index] = ll_sum
      return ll_store'''

class MultProbabilityStore:
    def __init__(self,mult_array_store,reads_correction_store,major_cn,minor_cn,n_subclones):
        self.non_phased_array = mult_array_store['Non_Phased']
        self.major_array = mult_array_store['Major']
        self.minor_array = mult_array_store['Minor']
        self.combined_array = mult_array_store['All']

        
        self.reads_correction_non_phased_array = reads_correction_store['Non_Phased']
        self.reads_correction_major_array = reads_correction_store['Major']
        self.reads_correction_minor_array = reads_correction_store['Minor']
        self.reads_correction_combined_array = reads_correction_store['All']

        self.use_non_phased = self.non_phased_array is not None and self.non_phased_array.shape[0]>0
        self.use_major = self.major_array is not None and self.major_array.shape[0]>0
        self.use_minor = self.minor_array is not None and self.minor_array.shape[0]>0

       

        if self.use_non_phased:

            assert self.non_phased_array.shape[1] == self.reads_correction_non_phased_array.size
        if self.use_major:

            assert self.major_array.shape[1] == self.reads_correction_major_array.size
        if self.use_minor:

            assert self.minor_array.shape[1] == self.reads_correction_minor_array.size

        if sum([self.use_non_phased,self.use_major,self.use_minor]) == 0:
            print(self.non_phased_array,self.major_array,self.minor_array)
            raise ValueError('No arrays provided')

        self.n_mutations = self.non_phased_array.shape[0] if self.use_non_phased else self.major_array.shape[0] if self.use_major else self.minor_array.shape[0]
        
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.n_subclones = n_subclones
    
    def sample_array(self,array):
        return array[np.random.choice(array.shape[0],size=array.shape[0],replace=True),:]
    
    def __str__(self):
        return str(np.average(self.non_phased_array,axis=0))
    
    def get_subclonal_correction_array(self,subclone_table):
        weights = []
        arrays = []
        if self.use_non_phased:
            non_phased_array = np.concatenate([[self.reads_correction_non_phased_array[0]],self.reads_correction_non_phased_array[-(len(subclone_table.index)):]])
            arrays.append(non_phased_array)
            weights.append(self.non_phased_array.shape[0])
        if self.use_major:
            major_array = np.concatenate([[self.reads_correction_major_array[0]],self.reads_correction_major_array[-(len(subclone_table.index)):]])
            arrays.append(major_array)
            weights.append(self.major_array.shape[0])
        if self.use_minor:
            minor_array = np.concatenate([[self.reads_correction_minor_array[0]],self.reads_correction_minor_array[-(len(subclone_table.index)):]])
            arrays.append(minor_array)
            weights.append(self.minor_array.shape[0])
        arrays = np.stack(arrays)
        weights = np.array(weights)
        weights = weights/np.sum(weights)
        return np.average(arrays,axis=0,weights=weights)
    def sample(self):
        non_phased_array = self.sample_array(self.non_phased_array)
        major_array = self.sample_array(self.major_array)
        minor_array = self.sample_array(self.minor_array)
        return MultProbabilityStore(non_phased_array,major_array,minor_array,self.major_cn,self.minor_cn,self.n_subclones)

    def array_log_likelihood(self,clonal_mult_state,subclonal_mult_state,mult_array):

        if mult_array.shape[0] ==1 and np.isclose(np.sum(mult_array),0):
            return 0
        if not np.isclose(np.sum(clonal_mult_state),0):
            clonal_mult_state = clonal_mult_state/np.sum(clonal_mult_state)
        mult_state = np.concatenate((clonal_mult_state,subclonal_mult_state))
        mult_state = np.clip(mult_state,0.0,1.0)
        mult_state = mult_state/np.sum(mult_state)
        
        log_likelihood = np.multiply(mult_state, mult_array)
        
        log_likelihood = np.sum(log_likelihood, axis=1)
        log_likelihood = np.sum(np.log(log_likelihood + 2.2e-300))

        return log_likelihood
    

    def array_likelihood(self,mult_state,array):
        likelihood = np.multiply(mult_state, array)
        
        likelihood = np.sum(likelihood, axis=1)
        likelihood = np.sum(np.log(likelihood + 2.2e-300))
        return likelihood
    def evaluate_likelihood(self,mult_state):
        mult_state = np.clip(mult_state,0.0,1.0)
        mult_state = mult_state/np.sum(mult_state)
        
        log_likelihood= self.array_likelihood(mult_state,self.combined_array)

        return log_likelihood
    
    def evaluate_likelihood_array(self,full_states):
        start_time = time.time()
        #to do - implement phasing
        #added three reads correction
        subclonal_indicies = np.arange(full_states.shape[1]-self.n_subclones,full_states.shape[1])
        non_phased_indicies = np.arange(0,self.major_cn)
        major_indicies = np.arange(self.major_cn,2*self.major_cn)
        minor_indicies = np.arange(2*self.major_cn,2*self.major_cn+self.minor_cn)
        assert np.all(np.concatenate((non_phased_indicies,major_indicies,minor_indicies,subclonal_indicies)) == np.arange(full_states.shape[1]))
        

        
        non_phased_indicies = np.concatenate((non_phased_indicies,subclonal_indicies))
        major_indicies = np.concatenate((major_indicies,subclonal_indicies))
        minor_indicies = np.concatenate((minor_indicies,subclonal_indicies))

        ll = np.zeros(full_states.shape[0])
        
        if self.use_non_phased:
            full_states_non_phased = full_states[:,non_phased_indicies]
            ll_non_phased = evaluate_likelihood_array_numba(full_states_non_phased,self.non_phased_array,self.reads_correction_non_phased_array)
            
            ll += ll_non_phased
        if self.use_major:
            full_states_major = full_states[:,major_indicies]
            ll_major = evaluate_likelihood_array_numba(full_states_major,self.major_array,self.reads_correction_major_array)
           
            ll += ll_major
        if self.use_minor:
            full_states_minor = full_states[:,minor_indicies]
            
            ll_minor = evaluate_likelihood_array_numba(full_states_minor,self.minor_array,self.reads_correction_minor_array)
            ll += ll_minor
        return ll
    

class Sample:
   
    def __init__(self,mutation_table,cn_table,subclone_table,sample_id,purity,sex=None,merge_cn=True,apply_reads_correction=True):
        
        self.chromosome_order  = list(map(str,range(1,23)))
        self.chromosome_order.extend(['X','Y'])
        self.sample_id = sample_id
        self.purity = purity
        self.sex = sex
        assert sex in [None,'XX','XY']
        self.merge_cn = merge_cn
        self.apply_reads_correction = apply_reads_correction
        self.copy_number_table = self.process_raw_copy_number_table(cn_table)
        self.mutation_table = self.process_raw_mutation_table(mutation_table,self.copy_number_table)

        self.cn_table = cn_table
        self.n_snvs = len(self.mutation_table.index)
        self.subclone_table = self.format_subclone_table(subclone_table)
        self.segments = self.get_segments()
    
    def process_raw_mutation_table(self,mutation_table,cn_table):
        mutation_table = mutation_table.copy()
        
        mutation_table['Chromosome'] = mutation_table['Chromosome'].str.replace('chr','')
        

        mutation_table = dataloader.assign_cn_to_snv(mutation_table,cn_table)
       
        if self.sex is None:
            mutation_table = dataloader.filter_sex_chromosomes(mutation_table)
            
        mutation_table = self.filter_mutation_table(mutation_table)
     
        
        mutation_table = self.format_mutation_table(mutation_table)

        
        return mutation_table
    def process_raw_copy_number_table(self,cn_table):
        cn_table = cn_table.copy()
        cn_table['Chromosome'] = cn_table['Chromosome'].str.replace('chr','')
        if self.merge_cn:
            cn_table = dataloader.merge_segments(cn_table)
        cn_table['Total_CN'] = cn_table['Major_CN']+cn_table['Minor_CN']
        cn_table = cn_table[cn_table['Major_CN']>0]
        return cn_table
    def format_subclone_table(self,subclone_table):
        if subclone_table is None:
            return None
        
        subclone_table = subclone_table.copy()
        subclone_table = dataloader.get_valid_subclones(subclone_table)

        
        subclone_table = dataloader.filter_excess_subclones(subclone_table)
        if not 'N_SNVs' in list(subclone_table.columns):
            n_snvs = np.round(len(self.mutation_table.index)*subclone_table['Subclone_Fraction']).astype(int)
            subclone_table['N_SNVs'] = n_snvs
        return subclone_table
    
    def format_mutation_table(self,mutation_table):
        if 'Phasing' not in mutation_table.columns:
            mutation_table['Phasing'] = np.nan
        if 'Context' not in mutation_table.columns:
            mutation_table['Context'] = np.nan
        mutation_table = mutation_table[['Segment_ID', 'Chromosome','Segment_Start', 'Segment_End', 'Major_CN','Minor_CN', 'Total_CN','Tumor_Ref_Count', 'Tumor_Alt_Count', 'Position','Phasing','Context']].copy()
        mutation_table['Chromosome'] = mutation_table['Chromosome'].astype(str)
        mutation_table["Gain_Type"] = mutation_table["Major_CN"].astype(str)+ "_"+ mutation_table["Minor_CN"].astype(str)
        return mutation_table
    
    def filter_mutation_table(self,mutation_table,min_alt_count=3,min_coverage=10):
        if mutation_table['Tumor_Alt_Count'].min() > 10:
            warnings.warn("There are no mutations with less than 10 alt reads. This may indicate a higher threshold for mutation calling than the default of 3 alt reads assumed by GRITIC")
        if len(mutation_table.index) ==0:
            raise ValueError("There are no mutations in the mutation table. Please check the mutation table is formatted correctly")
        mutation_table['Total_Count'] = mutation_table['Tumor_Ref_Count']+mutation_table['Tumor_Alt_Count']
        

        mutation_table = mutation_table[mutation_table['Tumor_Alt_Count']>=min_alt_count]
        mutation_table = mutation_table[(mutation_table['Tumor_Ref_Count']+mutation_table['Tumor_Alt_Count'])>=min_coverage]
        
        return mutation_table
    
    def get_segments(self,min_mutations=0):
        segments = []
        for segment_id,segment_table in self.mutation_table.groupby('Segment_ID'):
            if len(segment_table.index) >= min_mutations:
                segment = Segment(segment_table,self.subclone_table,self.purity,self.sex,apply_reads_correction=self.apply_reads_correction)
                segments.append(segment)
        return segments
    
    def get_segments_with_copy_number(self,major_cn=None,minor_cn=None):
        segments_with_copy_number = []
        for segment in self.segments:
            if major_cn is not None and segment.major_cn != major_cn:
                continue
            if minor_cn is not None and segment.minor_cn != minor_cn:
                continue
            segments_with_copy_number.append(segment)
        return segments_with_copy_number
    
    def get_mutation_table(self):
        mutation_table = pd.concat([segment.mutation_table for segment in self.segments])
        return self.sort_table_by_chromosome(mutation_table)
    
    def get_subclone_table(self):
        return self.subclone_table
    
    
    def sort_table_by_chromosome(self,table):
        table =table.copy()
        table['Chromosome'] = pd.Categorical(table['Chromosome'],categories=self.chromosome_order)
        return table.sort_values(by=['Chromosome','Segment_Start']).reset_index(drop=True)

class Segment:
    
    def __init__(self,mutation_table,subclone_table,purity,sex,apply_reads_correction=True):
        self.mutation_table = mutation_table
        self.n_snvs = len(mutation_table.index)
        self.n_mutations = len(self.mutation_table.index)
        self.subclone_table = subclone_table
        self.sex = sex
        self.sample_clone_fractions = self.get_sample_clone_fractions()
        self.n_subclones = self.get_n_subclones()
        self.apply_reads_correction = apply_reads_correction
        self.purity = purity

        self.segment_id = self.get_unique_attribute_from_table('Segment_ID')

        
        
        self.major_cn = self.get_unique_attribute_from_table('Major_CN')
        self.minor_cn = self.get_unique_attribute_from_table('Minor_CN')
        self.total_cn = self.major_cn + self.minor_cn

        self.possible_clonal_multiplicities = np.arange(1,self.major_cn+1)

        self.all_possible_clonal_multiplicities = np.arange(self.major_cn)+1
        self.all_possible_subclonal_multiplicities = self.get_all_possible_subclonal_multiplicities()
        
        self.chromosome = self.get_unique_attribute_from_table('Chromosome')
        self.start = self.get_unique_attribute_from_table('Segment_Start')
        self.end = self.get_unique_attribute_from_table('Segment_End')
        self.width = self.end - self.start
        
        self.timeable = False
        self.in_wgd_range = False
        self.gain_timing_distribution = None
        self.complex_timing_distributions = None
        
        self.phased = not self.mutation_table['Phasing'].isnull().all()
        self.assign_multiplicity_probabilities()
        
        self.multiplicity_probabilities = self.get_multiplicity_probabilities()

        
    
    def __str__(self):
        return "{}- {}+{} - {} Mutations".format(self.segment_id,self.major_cn,self.minor_cn,self.n_mutations)
    
    def get_all_possible_subclonal_multiplicities(self):
        if self.n_subclones ==0:
            return np.zeros(0)
        return self.subclone_table['Subclone_CCF'].to_numpy()
    
    def get_n_subclones(self):
        if self.subclone_table is None:
            return 0
        return len(self.subclone_table.index)
    
    def get_sample_clone_fractions(self):
        if self.subclone_table is None:
            return np.array([1])
        subclone_fractions = self.subclone_table['Subclone_Fraction'].to_numpy()
        total_subclone_fraction = np.sum(subclone_fractions)
        clone_fraction =  1 - total_subclone_fraction
        sample_clone_fractions = np.insert(subclone_fractions,0,clone_fraction)
        return sample_clone_fractions
    
    def get_unique_attribute_from_table(self,attribute):
        if len(self.mutation_table[attribute].unique()) >1:
            raise ValueError("Attribute {} is not unique for segment".format(attribute))
        return self.mutation_table[attribute].iloc[0]

    
    def get_multiplicity_names(self):
        sample_peak_names = ['Mult_{}'.format(mult) for mult in self.all_possible_clonal_multiplicities]
        for i in range(self.n_subclones):
            sample_peak_names.append("Subclone_{}".format(i))
        return sample_peak_names
    
    def get_mutation_rate(self):
        mult_proportions = multiplicityoptimiser.unconstrained_mult_optimisation(self.multiplicity_probabilities,self.n_subclones)
        if mult_proportions is None:
            return np.nan
        mult_states = np.concatenate([np.arange(self.major_cn)+1,np.ones(self.n_subclones)])
        
        mult_states = np.divide(mult_states,self.multiplicity_probabilities.reads_correction_combined_array)
        mult_correction_factor = np.sum(np.multiply(mult_proportions,mult_states))
        
        mutation_rate = self.n_mutations*mult_correction_factor
        return mutation_rate/self.width
    def assign_multiplicity_probabilities(self):
        #remove mult cols if they already exist
        mult_cols = [col for col in self.mutation_table.columns if 'Prob_' in col or 'Three_Reads_Correction_' in col]
        self.mutation_table = self.mutation_table.drop(columns=mult_cols)
        if 'X' in set(self.mutation_table['Chromosome'].unique()) or 'Y' in set(self.mutation_table['Chromosome'].unique()):
            assert self.sex is not None
        normal_total_cn = self.mutation_table['Chromosome'].map(lambda x: 1 if (x =='Y') or (x=='X' and self.sex=='XY') else 2)
        
        three_reads_correction_factors = []
        multiplicity_probabilities = []

        mult_one_vaf = self.purity / (self.total_cn* self.purity + normal_total_cn * (1 - self.purity))
        
        combined_all_possible_multiplicities = np.concatenate((self.all_possible_clonal_multiplicities,self.all_possible_subclonal_multiplicities))
        peak_names = self.get_multiplicity_names()
        
        vaf = self.mutation_table["Tumor_Alt_Count"] / (self.mutation_table["Tumor_Alt_Count"] + self.mutation_table["Tumor_Ref_Count"])
        highest_vaf_m_table = self.mutation_table[vaf>np.percentile(vaf,90)-0.01]
        highest_vaf_average_coverage = np.average(highest_vaf_m_table['Tumor_Alt_Count']+highest_vaf_m_table['Tumor_Ref_Count'])
        
        if not self.apply_reads_correction:
            highest_vaf_average_coverage = (self.mutation_table['Tumor_Alt_Count']+self.mutation_table['Tumor_Ref_Count']).mean()

        for multiplicity in combined_all_possible_multiplicities:
            # if multiplicity is bigger than total cn then the vaf can be > 1 which causes a crash
            #any mult > major cn is set to zero subsequently
            mult_vaf = np.minimum(1, mult_one_vaf * multiplicity)

            total_counts = self.mutation_table["Tumor_Alt_Count"] + self.mutation_table["Tumor_Ref_Count"]
            mult_probability = binom.pmf(self.mutation_table["Tumor_Alt_Count"],total_counts,mult_vaf)
            multiplicity_probabilities.append(mult_probability)

            three_reads_correction_factor = 1- binom.cdf(2,np.random.poisson(highest_vaf_average_coverage,size=mult_vaf.size),mult_vaf)
            #three_reads_correction_factor = np.random.binomial(np.random.poisson(highest_vaf_average_coverage,size=mult_vaf.size),mult_vaf,size=mult_vaf.size)
            three_reads_correction_factors.append(three_reads_correction_factor)

        multiplicity_probabilities = np.array(multiplicity_probabilities)
        three_reads_correction_factors = np.array(three_reads_correction_factors)
        
      
        normalising_sums = np.sum(multiplicity_probabilities, axis=0)
        normalising_sums = np.where(np.isclose(normalising_sums,0),1,normalising_sums)
        multiplicity_probabilities = multiplicity_probabilities / normalising_sums
        new_cols = {}
        for index, peak_name in enumerate(peak_names):
            new_cols[f"Prob_{peak_name}"] = multiplicity_probabilities[index, :]
            new_cols[f'Three_Reads_Correction_{peak_name}'] = three_reads_correction_factors[index,:]
        new_cols_table = pd.DataFrame(new_cols,index=self.mutation_table.index)
        self.mutation_table = pd.concat((self.mutation_table,new_cols_table),axis=1)
        
    def get_reads_correction_array(self,allele=None):
        if allele == 'Minor':
            clonal_multiplicities = np.arange(1,self.minor_cn+1)
        else:
            clonal_multiplicities = np.arange(1,self.major_cn+1)
        mult_names = ["Three_Reads_Correction_Mult_{}".format(mult) for mult in clonal_multiplicities]
        mult_names.extend(["Three_Reads_Correction_Subclone_{}".format(subclone) for subclone in range(self.n_subclones)])
 
        reads_correction = self.mutation_table[mult_names].to_numpy()

        if allele is None:
            reads_correction= reads_correction[self.mutation_table['Phasing'].isna(),:]
        elif allele =='All':
            reads_correction = reads_correction
        else:
            reads_correction = reads_correction[self.mutation_table['Phasing']==allele,:]
        
        if reads_correction.size ==0:
            return None

        
        
        
        reads_correction = np.average(reads_correction,axis=0)
        if not self.apply_reads_correction:
            reads_correction = np.ones_like(reads_correction)
        
        
        return reads_correction
    
    def get_multiplicity_probabilities_array(self,allele=None):

        if allele == 'Minor':
            clonal_multiplicities = np.arange(1,self.minor_cn+1)
        else:
            clonal_multiplicities = np.arange(1,self.major_cn+1)
        
        mult_names = [f"Prob_Mult_{mult}" for mult in clonal_multiplicities]
        
        mult_names.extend(["Prob_Subclone_{}".format(subclone) for subclone in range(self.n_subclones)])

        multiplicity_probabilities = self.mutation_table[mult_names].to_numpy()

        normalising_sums = np.sum(multiplicity_probabilities, axis=1)[:,None]
        normalising_sums = np.where(np.isclose(normalising_sums,0),1,normalising_sums)
        multiplicity_probabilities = multiplicity_probabilities / normalising_sums

   
  

        if allele is None:     
            multiplicity_probabilities= multiplicity_probabilities[self.mutation_table['Phasing'].isna(),:]
        elif allele =='All':
            multiplicity_probabilities = multiplicity_probabilities
        else:
            multiplicity_probabilities = multiplicity_probabilities[self.mutation_table['Phasing']==allele,:]

        if multiplicity_probabilities.size ==0:
            return None
        
        return multiplicity_probabilities
    

    def get_multiplicity_probabilities(self):
        reads_correction_store = {}
        reads_correction_store['Non_Phased'] = self.get_reads_correction_array()
        reads_correction_store['Major'] = self.get_reads_correction_array('Major')
        reads_correction_store['Minor'] = self.get_reads_correction_array('Minor')
        reads_correction_store['All'] = self.get_reads_correction_array('All')

        mult_array_store = {}
        mult_array_store['Non_Phased'] = self.get_multiplicity_probabilities_array()

        mult_array_store['Major'] = self.get_multiplicity_probabilities_array('Major')
        
        mult_array_store['Minor'] = self.get_multiplicity_probabilities_array('Minor')
        mult_array_store['All'] = self.get_multiplicity_probabilities_array('All')
        
        return MultProbabilityStore(mult_array_store,reads_correction_store,self.major_cn,self.minor_cn,self.n_subclones)
    
    def get_info_dict(self):
        
        info_dict = {'Segment_ID':self.segment_id,'Chromosome':self.chromosome,'Segment_Start':self.start,'Segment_End':self.end}
        info_dict['Major_CN'] = self.major_cn
        info_dict['Minor_CN'] = self.minor_cn
        info_dict['Total_CN'] = self.total_cn
        info_dict['N_Mutations']=self.n_mutations
        info_dict['Mutation_Rate'] = self.get_mutation_rate()
        return info_dict
    
  
   
