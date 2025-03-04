import os
import bz2
import pickle

import warnings

import pandas as pd
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

from scipy.special import logsumexp


from scipy.optimize import nnls
from scipy.linalg import null_space

from gritic.sampletools import Segment, get_major_cn_mode


from sklearn.neighbors import NearestNeighbors

import time
import shutil

import gritic.distributiontools as distributiontools
import gritic.treetools as treetools
import gritic.hitandrun as hitandrun

import gritic.posteriortablegen as posteriortablegen

import pathlib
global_test_vars = {'density_cut_off':0.9,'prior_weight':1}
class RouteTree:
    def __init__(self,tree,major_cn,minor_cn,wgd_status):
        self.main_tree = tree
        
        self.major_cn = major_cn
        self.minor_cn = minor_cn

        major_tree,minor_tree = treetools.split_tree(tree)
        self.major_tree = major_tree
        self.minor_tree = minor_tree

        self.wgd_status = wgd_status
        self.node_attributes = treetools.get_node_attributes(self.main_tree,self.wgd_status)

        self.non_phased_node_order = self.get_node_order()
        
        self.timeable_nodes = self.get_timeable_nodes()
        self.wgd_nodes = self.get_wgd_nodes()
        self.sum_constraint_matrix = self.get_sum_constraint_matrix()
        self.wgd_constraint_matrix = self.get_wgd_constraint_matrix()

        self.timing_matrix = self.get_timing_matrix()
        self.unphased_timing_matrix= self.get_unphased_timing_matrix()
    

    def get_timeable_nodes(self):
        timeable_nodes = []
        for node in self.non_phased_node_order:
            if len(self.main_tree.out_edges(node))==2 and not self.main_tree.nodes[node]['WGD_Symbol']:
                timeable_nodes.append(node)
        return timeable_nodes
    def get_non_terminal_nodes(self):
        non_terminal_nodes = []
        for node in self.non_phased_node_order:
            if len(self.main_tree.out_edges(node))>0:
                non_terminal_nodes.append(node)
        return non_terminal_nodes
    def get_terminal_nodes(self):
        terminal_nodes = []
        for node in self.non_phased_node_order:
            if len(self.main_tree.out_edges(node))==0:
                terminal_nodes.append(node)
        return terminal_nodes
    def get_wgd_nodes(self):
        wgd_nodes = []
        for node in self.non_phased_node_order:
            if self.main_tree.nodes[node]['WGD_Symbol']:
                wgd_nodes.append(node)
        return wgd_nodes
    #all nodes
    def get_node_order(self,allele=None):
        if allele is None:
            return list(nx.dfs_preorder_nodes(self.main_tree))
        if allele == 'Major':
            return list(nx.dfs_preorder_nodes(self.major_tree))
        if allele == 'Minor':
            return list(nx.dfs_preorder_nodes(self.minor_tree))
        raise ValueError(f'allele must be None, Major or Minor, not {allele}')
        
       
    #an mxn matrix
    #m = number of timing nodes
    #n is the number of multiplicity states
    def get_timing_matrix(self):

        n_mults = self.major_cn*2+self.minor_cn

        mult_offset = 1
        
        timing_matrix = np.zeros((len(self.non_phased_node_order),n_mults))
        
        for i,node in enumerate(self.non_phased_node_order):
            final_mult = self.node_attributes[node]['Multiplicity']
            timing_matrix[i,final_mult-mult_offset]=1

            if node in self.major_tree.nodes():
                timing_matrix[i,final_mult-mult_offset+self.major_cn] = 1
            if node in self.minor_tree.nodes():
                timing_matrix[i,final_mult-mult_offset+2*self.major_cn] = 1
        return timing_matrix
    
    def get_unphased_timing_matrix(self):
        n_mults = self.major_cn
        mult_offset = 1
        
        timing_matrix = np.zeros((len(self.non_phased_node_order),n_mults))
        
        for i,node in enumerate(self.non_phased_node_order):
            final_mult = self.node_attributes[node]['Multiplicity']
            timing_matrix[i,final_mult-mult_offset]=1
        
        return timing_matrix
    #an mxn matrix
    #m = number of paths through the tree
    #n number of timing nodes
    def get_sum_constraint_matrix(self):
        possible_paths = treetools.get_possible_paths(self.main_tree)
        sum_constraint_matrix = np.zeros((len(possible_paths),len(self.non_phased_node_order)))
        
        for row,path in enumerate(possible_paths):
            for node in path:
                sum_constraint_matrix[row,self.non_phased_node_order.index(node)] = 1
        
        return sum_constraint_matrix
    
    def get_wgd_constraint_matrix(self):
        
        wgd_paths = treetools.get_wgd_paths(self.main_tree)
        if len(wgd_paths)==0:
            return None
        
        constraint_matrix = np.zeros((len(wgd_paths),len(self.non_phased_node_order)))
        
        for row,path in enumerate(wgd_paths):
            for node in path:
                constraint_matrix[row,self.non_phased_node_order.index(node)] = 1
        return constraint_matrix
    
    def get_combined_constraints(self,wgd_timing):
        combined_constraint_matrix = self.sum_constraint_matrix
        combined_constraints_sum = np.ones(combined_constraint_matrix.shape[0])

        if self.wgd_constraint_matrix is not None:
            combined_constraint_matrix = np.vstack([combined_constraint_matrix,self.wgd_constraint_matrix])
            wgd_constraints_sum = np.ones(self.wgd_constraint_matrix.shape[0])*wgd_timing
            combined_constraints_sum = np.concatenate([combined_constraints_sum,wgd_constraints_sum])

        return combined_constraint_matrix,combined_constraints_sum

    def get_n_events(self,node_timing,wgd_timing):
        n_events_store = []
        if wgd_timing is None:
            #if no wgd the number of events is the number of gains
            n_events = len([node for node in self.main_tree.nodes if len(list(self.main_tree.successors(node)))==2])
            return n_events,np.nan,np.nan
        
        pre_wgd_losses =0
        plotting_tree = treetools.convert_to_plotting_tree(self.main_tree,wgd_timing,node_timing,self.non_phased_node_order)
        post_wgd_losses = sum([int(plotting_tree.nodes[node]['Loss_Symbol']) for node in plotting_tree.nodes])

        n_wgds = sum([int(self.main_tree.nodes[node]['WGD_Symbol']) for node in self.main_tree.nodes])
        n_gains = len([node for node in self.main_tree.nodes if len(list(self.main_tree.successors(node)))==2])
        #each wgd is only one event
        n_events = n_gains-(n_wgds-1)+post_wgd_losses
        #add extra loss not accounted for
        if self.minor_cn ==0:
            n_events +=1
            pre_wgd_losses = 1
        return n_events,pre_wgd_losses,post_wgd_losses

class Route:
    def __init__(self,route_id,tree,major_cn,minor_cn,wgd_status,mult_store_dir):
        self.route_id = route_id
        self.short_id = route_id[:9]
        self.route_tree = RouteTree(tree,major_cn,minor_cn,wgd_status)
        
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.total_cn = self.major_cn + self.minor_cn

        self.wgd_status = wgd_status
        self.mult_store_dir = mult_store_dir
        self.ll_store = None
        self.node_timing = None
        self.loss_timing = None
        self.wgd_timing_store = None
        self.n_events_store = None
        self.mult_store = None
        self.density = None
        self.density_high = None

        self.run_time = np.nan
    
    
    def get_average_events(self,event_type):
        if self.n_events_store is None:
            return np.nan
        return np.mean(self.n_events_store[event_type])

    
    def get_node_timing(self,node):
        if self.node_timing is None:
            return np.nan
        return self.node_timing[self.route_tree.non_phased_node_order.index(node),:]
    
    def get_cumulative_timing(self,timing_periods):
        cumulative_timing = []
        for i,node in enumerate(self.route_tree.non_phased_node_order):

            predecessor = self.route_tree.node_attributes[node]['Predecessor']
            if predecessor is None:
                cumulative_timing.append(timing_periods[:,i])
            else:
                cumulative_timing.append(cumulative_timing[self.route_tree.non_phased_node_order.index(predecessor)]+timing_periods[:,i])
        return np.array(cumulative_timing)
        
    
    def get_weighted_arrays(self,cumulative_timing,wgd_timing_store,mult_store,weights,n_samples=1000):
        #allow some tolerance
        if np.isnan(weights).any():
            cumulative_timing = np.ones_like(cumulative_timing)[:,:n_samples]*np.nan
            wgd_timing_store = np.ones_like(wgd_timing_store)[:n_samples]*np.nan
            mult_store = np.ones_like(mult_store)[:n_samples,:]*np.nan
            return cumulative_timing,wgd_timing_store,mult_store
        assert (weights >-1e-80).all()
        weights = np.clip(weights,0,1)
        weights = weights/np.sum(weights)
        #weighted_sample = np.random.choice(np.arange(cumulative_timing.shape[1]),size=cumulative_timing.shape[1],replace=True,p=weights)
        weighted_sample = np.random.choice(np.arange(cumulative_timing.shape[1]),size=n_samples,replace=True,p=weights)
        return np.array(cumulative_timing)[:,weighted_sample],wgd_timing_store[weighted_sample],mult_store[weighted_sample]
    
    def sample_mults(self,wgd_timing,n_samples):
        constraints_matrix,constraints_sum = self.route_tree.get_combined_constraints(wgd_timing)
        nnls_time = time.perf_counter()
        start_sol = nnls(constraints_matrix, constraints_sum)[0]
        constraints_null = null_space(constraints_matrix)
        null_space_time = time.perf_counter()
        solutions = hitandrun.hit_and_run(constraints_null,start_sol,n_samples=n_samples)

        timing = self.get_cumulative_timing(solutions)
        
        mult = np.matmul(solutions,self.route_tree.timing_matrix)
        
        unphased_mult_sum = np.tile(np.sum(mult[:,:self.major_cn],axis=1),(self.major_cn,1)).T
        major_cn_mult_sum = np.tile(np.sum(mult[:,self.major_cn:2*self.major_cn],axis=1),(self.major_cn,1)).T
        minor_cn_mult_sum = np.tile(np.sum(mult[:,2*self.major_cn:],axis=1),(self.minor_cn,1)).T
        combined_mult_sum = np.concatenate([unphased_mult_sum,major_cn_mult_sum,minor_cn_mult_sum],axis=1)
        mult = mult/combined_mult_sum
    
        
       
    

        return mult,timing

    
    def get_n_events_estimate(self,node_timing,wgd_timing,n_samples = 300):
        n_events_estimate = {'N_Events':[],'Pre_WGD_Losses':[],'Post_WGD_Losses':[]}
        for i in range(n_samples):
            random_index = np.random.choice(node_timing.shape[1])
            n_events,pre_wgd_losses,post_wgd_losses = self.route_tree.get_n_events(node_timing[:,random_index],wgd_timing[random_index])
            n_events_estimate['N_Events'].append(n_events)
            n_events_estimate['Pre_WGD_Losses'].append(pre_wgd_losses)
            n_events_estimate['Post_WGD_Losses'].append(post_wgd_losses)

        return n_events_estimate 
    

    @staticmethod
    def simulate_clone_share(alpha,n_samples):
        alpha = alpha/np.sum(alpha)
        alpha = global_test_vars['prior_weight']*alpha +1

        dirichlet_sample = np.random.dirichlet(alpha,size=n_samples)
        return dirichlet_sample
    def get_density_estimate(self,samples,n_test_points=1000,radius=0.05):
        nn_finder = NearestNeighbors(radius=radius,p=1)
        nn_finder.fit(samples)
        random_mult_indicies = np.random.choice(samples.shape[0],size=n_test_points,replace=False)
        nearest_neighbors = nn_finder.radius_neighbors(samples[random_mult_indicies,:],return_distance=False)
        nn_size = np.array([x.size-1 for x in nearest_neighbors])
        return np.mean(nn_size>0.1),np.mean(nn_size>2.1)
    
    def run_mult_sampling(self,n_snvs,alpha,n_subclones,wgd_timing_distribution,samples_per_run=500,max_samples=5e5):
        timing_store = []
        wgd_timing_store= []
        mult_store =[]

        start_sampling_time = time.perf_counter()
        n_samples = 0
        next_eval_time = None
        eval_count =0
        start_time = time.perf_counter()
        while True:
            if self.wgd_status:
                wgd_timing = np.random.choice(wgd_timing_distribution)
            else:
                wgd_timing =np.nan
            mults,timing = self.sample_mults(wgd_timing,samples_per_run)
            #integrate in the (likely case) that there's only one subclone
            '''if n_subclones ==-1:
                ll_store.append(mult_probabilities.evaluate_likelihood_array_one_subclone(mults))'''
            timing_store.append(timing)

            if n_subclones >0 :
                clone_share = self.simulate_clone_share(alpha,timing.shape[1])
                clonal_share = clone_share[:,0].reshape(-1,1)
                subclone_mults = clone_share[:,1:]
                
                mults = np.concatenate([mults*clonal_share,subclone_mults],axis=1)
            mult_store.append(mults)
            
            
            
           

            wgd_timing_store.extend([wgd_timing]*timing.shape[1])
            if eval_count==100 or (next_eval_time is not None and time.perf_counter()>next_eval_time):
                
                timing_test = np.concatenate(timing_store,axis=1)
                timing_test = np.transpose(timing_test)
         
                if n_subclones >0:
                    subclone_mults = np.concatenate(mult_store,axis=0)[:,-n_subclones:]
                    timing_test = np.concatenate([timing_test,subclone_mults],axis=1)
                sampling_time = time.perf_counter()-start_sampling_time
                
                start_density_time = time.perf_counter()
                density,density_high = self.get_density_estimate(timing_test)
                density_time = time.perf_counter()-start_density_time
                if eval_count >50:
                    next_eval_time = time.perf_counter()+min(max(density_time*5,sampling_time),30.0)
                else:
                    next_eval_time = time.perf_counter()+1.0
                
                start_sampling_time = time.perf_counter()
                if density >= global_test_vars['density_cut_off']:
                    break
            eval_count+=1
            if eval_count*samples_per_run > max_samples:
                density,density_high = self.get_density_estimate(timing_test)
                break
            
        timing_store =np.concatenate(timing_store,axis=1)
        mult_store = np.concatenate(mult_store,axis=0)
        
        wgd_timing_store = np.array(wgd_timing_store)

        return mult_store,timing_store,wgd_timing_store,np.array([density,density_high])

    def save_gz_numpy(self,path,array):
        with open(path,'wb') as f:
            np.save(f,array)
    def load_gz_numpy(self,path):
        with open(path,'rb') as f:
            return np.load(f)
    
    def get_mult_store(self,n_snvs,alpha,n_subclones,wgd_timing_distribution):
        if self.mult_store_dir is not None:
            cn_dir = self.mult_store_dir/f'{self.major_cn}_{self.minor_cn}'
            os.makedirs(cn_dir,exist_ok=True)
            


            #pathlib version
            mult_store_path = cn_dir/f'mult_store_{self.short_id}_wgd_{self.wgd_status}.npy.gz'
            timing_store_path = cn_dir/f'timing_store_{self.short_id}_wgd_{self.wgd_status}.npy.gz'
            wgd_timing_store_path = cn_dir/f'wgd_timing_store_{self.short_id}_wgd_{self.wgd_status}.npy.gz'
            density_store_path = cn_dir/f'density_store_{self.short_id}_wgd_{self.wgd_status}.npy.gz'

            if os.path.exists(mult_store_path):
                mult_store = self.load_gz_numpy(mult_store_path)
                timing_store = self.load_gz_numpy(timing_store_path)
                
                wgd_timing_store = self.load_gz_numpy(wgd_timing_store_path)
                density = self.load_gz_numpy(density_store_path)


                return mult_store,timing_store,wgd_timing_store,density
        mult_store,timing_store,wgd_timing_store,density = self.run_mult_sampling(n_snvs,alpha,n_subclones,wgd_timing_distribution)
        
        if self.mult_store_dir is not None:
            self.save_gz_numpy(mult_store_path,mult_store)
            self.save_gz_numpy(timing_store_path,timing_store)
            self.save_gz_numpy(wgd_timing_store_path,wgd_timing_store)
            self.save_gz_numpy(density_store_path,density)
        
        return mult_store,timing_store,wgd_timing_store,density

    def get_raw_samples_store(self,mult_store,timing_store,wgd_timing_store,ll_store,n_samples=10000):
        
        random_indexes = np.random.randint(0,mult_store.shape[0],size=n_samples)
        
        raw_samples_store = {'Timing':{},'Mult':[],'WGD_Timing':[],'LL':[]}
        for node in self.route_tree.timeable_nodes:
            node_timing = timing_store[self.route_tree.non_phased_node_order.index(node),random_indexes].copy()
            raw_samples_store['Timing'][node] = node_timing
        raw_samples_store['Mult'] = mult_store[random_indexes,:].copy()
        raw_samples_store['WGD_Timing'] = wgd_timing_store[random_indexes].copy()
        raw_samples_store['LL'] = ll_store[random_indexes].copy()
        return raw_samples_store
    
    def run_sampling(self,mult_probabilities,subclone_table,n_snvs,wgd_timing_distribution,phased):

        run_time = time.perf_counter()

        #alpha os subclonal clonal proportions
        alpha = None
        if subclone_table is not None:
            n_clonal_snvs_sample = np.round(subclone_table['N_SNVs'].sum()/subclone_table['Subclone_Fraction'].sum()*(1-subclone_table['Subclone_Fraction'].sum())).astype(int)
            alpha = np.concatenate([[n_clonal_snvs_sample],subclone_table['N_SNVs']])+1
            subclonal_correction_array = mult_probabilities.get_subclonal_correction_array(subclone_table)
            
        
            alpha = alpha/subclonal_correction_array
            
        
        n_subclones = 0 if subclone_table is None else len(subclone_table.index)

        mult_store,timing_store,wgd_timing_store,density = self.get_mult_store(n_snvs,alpha,n_subclones,wgd_timing_distribution)

        ll_store = mult_probabilities.evaluate_likelihood_array(mult_store)
        self.raw_samples = self.get_raw_samples_store(mult_store,timing_store,wgd_timing_store,ll_store)
        
        
        weights = np.exp(ll_store-np.max(ll_store))
        node_timing,wgd_timing_store,mult_store= self.get_weighted_arrays(timing_store,wgd_timing_store,mult_store,weights)
        
        self.n_events_store = self.get_n_events_estimate(node_timing,wgd_timing_store)

        #self.get_n_events_estimate_full(node_timing,wgd_timing_store)
        self.ll_store = ll_store
        self.node_timing = node_timing
        self.wgd_timing_store = np.array(wgd_timing_store)
        self.mult_store = mult_store
        self.run_time = time.perf_counter()-run_time
        self.density = density[0]
        self.density_high = density[1]


class RouteClassifier:
    def __init__(self,major_cn,minor_cn,wgd_status,wgd_trees_status,mult_store_dir):
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.wgd_status = wgd_status
        self.mult_store_dir = mult_store_dir
        self.routes = self.generate_routes(wgd_trees_status)

        self.min_n_events = None

        self.route_probabilities = {}

    def get_best_timing(self):
        best_route_id = max(self.route_probabilities, key=self.route_probabilities.get)
        best_route = self.routes[best_route_id]
        best_timing = []
        for node in best_route.route_tree.timeable_nodes:
            best_timing.append(best_route.get_node_timing(node))
        return np.array(best_timing)
    
    def generate_routes(self,wgd_trees_status):
        possible_routes = {}
        possible_trees = treetools.get_nx_trees(self.major_cn,self.minor_cn,self.wgd_status,wgd_trees_status)
        for tree_id,tree in possible_trees.items():
            route = Route(tree_id,tree,self.major_cn,self.minor_cn,self.wgd_status,self.mult_store_dir)
            possible_routes[tree_id] = route
        return possible_routes
 
    def fit_routes(self,mult_probabilities,subclone_table,n_snvs,wgd_timing_distribution,phased,non_parsimony_penalty=False):
        route_ll_store = []
        route_ids = list(sorted(self.routes.keys()))
        for route_id in route_ids:
            route = self.routes[route_id]
            route.run_sampling(mult_probabilities,subclone_table,n_snvs,wgd_timing_distribution,phased)
            route_ll_store.append(route.ll_store)

        #helps keep the exponentiation under control
        max_point = np.max(np.concatenate(route_ll_store))
        likelihood_store =[]
        
        for route_ll in route_ll_store:
            route_likelihoods = np.exp(route_ll-max_point)
            average_likelihood = np.average(route_likelihoods)
            likelihood_store.append(average_likelihood)
        
        likelihood_store = np.array(likelihood_store)/np.sum(likelihood_store)
        if non_parsimony_penalty:
            route_n_events = np.array([self.routes[route_id].get_average_events('N_Events') for route_id in route_ids])
            likelihood_store = apply_non_parsimony_penalty(likelihood_store,route_n_events)
        
        for i,route_id in enumerate(route_ids):
            self.route_probabilities[route_id] = likelihood_store[i]
    
    
    def get_timing_table(self):
        if len(self.route_probabilities) ==0:
            raise ValueError("routes hasn't been fit yet")
        timing_table_data = {"Route":[],'Probability':[],'Node':[],'Node_Phasing':[],'Timing':[],'Timing_CI_Low':[],'Timing_CI_High':[],'Average_N_Events':[],'Average_Pre_WGD_Losses':[],'Average_Post_WGD_Losses':[],'Time':[],'Density':[],'Density_High':[]}
        for route_id,probability in self.route_probabilities.items():
            route = self.routes[route_id]
            timeable_nodes = route.route_tree.timeable_nodes
            
            if len(timeable_nodes)==0:
                timeable_nodes = [np.nan]
            
            for node in timeable_nodes:
                if not node is np.nan:
                    node_timing = route.get_node_timing(node)
                else:
                    node_timing = np.array([np.nan])
                
                timing_table_data['Route'].append(route.short_id)
                timing_table_data['Probability'].append(probability)
                timing_table_data['Node'].append(node)
                timing_table_data['Node_Phasing'].append(route.route_tree.node_attributes[node]['Phasing'])
                timing_table_data['Timing'].append(np.abs(np.round(np.percentile(node_timing,50),3)))
                timing_table_data['Timing_CI_Low'].append(np.abs(np.round(np.percentile(node_timing,2.5),3)))
                timing_table_data['Timing_CI_High'].append(np.abs(np.round(np.percentile(node_timing,97.5),3)))
                timing_table_data['Average_N_Events'].append(route.get_average_events('N_Events'))
                timing_table_data['Average_Pre_WGD_Losses'].append(route.get_average_events('Pre_WGD_Losses'))
                timing_table_data['Average_Post_WGD_Losses'].append(route.get_average_events('Post_WGD_Losses'))
                timing_table_data['Time'].append(np.round(route.run_time,3))
                timing_table_data['Density'].append(route.density)
                timing_table_data['Density_High'].append(route.density_high)

        return pd.DataFrame(timing_table_data)
    
    def get_timing_dict(self,n_samples = 5000):
        if len(self.route_probabilities) ==0:
            raise ValueError("routes hasn't been fit yet")
        timing_dict = {}
        
        for route_id,probability in self.route_probabilities.items():
            if np.isclose(probability,0):
                continue
            
            route = self.routes[route_id]
            route_samples = {'Timing':{}}
            
            wgd_timing_store =route.wgd_timing_store
            random_indexes = np.random.randint(0,wgd_timing_store.size,size=n_samples)
            route_samples['Timing']['WGD'] = wgd_timing_store[random_indexes]
            
            for node in route.route_tree.timeable_nodes:
                node_timing = route.get_node_timing(node)
                route_samples['Timing'][node] = node_timing[random_indexes]
            route_samples['Mult'] = route.mult_store[random_indexes]
            route_samples['LL'] = route.ll_store[random_indexes]
            route_samples['Raw_Samples'] = route.raw_samples
            timing_dict[route_id] = route_samples
        return timing_dict

    def get_timing_tree_labels(self,route,wgd_info):
        node_labels = {}
        for node in route.route_tree.non_phased_node_order:
            if node in route.route_tree.timeable_nodes:
                timing_dist = route.get_node_timing(node)
                timing = np.abs(np.round(np.percentile(timing_dist,50),3))
                timing_ci_low = np.abs(np.round(np.percentile(timing_dist,5),3))
                timing_ci_high = np.abs(np.round(np.percentile(timing_dist,95),3))
                node_labels[node] = f"{timing} - [{timing_ci_low},{timing_ci_high}]"
            elif node in route.route_tree.wgd_nodes:
                node_labels[node]= f"{wgd_info['Timing']} - [{wgd_info['CI_Low']},{wgd_info['CI_High']}]"
            else:
                node_labels[node] = ''
        return node_labels
    
    def plot_trees(self,plot_output_dir,seg_title,wgd_info):
        os.makedirs(plot_output_dir,exist_ok=True)

        if len(self.route_probabilities.keys())==0:
            raise ValueError("The routes haven't been fit yet")  
        best_route_id = max(self.route_probabilities, key=self.route_probabilities.get)   
        for route_id,route in self.routes.items():
            
            probability = np.round(self.route_probabilities[route_id],4)
            
            route_output_path = f"{plot_output_dir}/route_{route.short_id}.pdf"
            tree_subtitle = f"Route {route.short_id} (Probability = {probability})"
            plot_title = f"{seg_title}\n{tree_subtitle}"
            if route_id==best_route_id:
                plot_title = f"{plot_title} - (Best Fit)"

            node_labels = self.get_timing_tree_labels(route,wgd_info)
            plotting_tree = route.route_tree.main_tree.copy()
            nx.set_node_attributes(plotting_tree,node_labels,'Label')
            
            treetools.plot_tree(plotting_tree,plot_title,output_path=route_output_path)
def check_record_segment_timing(segment,wgd_status,wgd_trees_status):
    if segment.major_cn <=1:
        return False
    return True


def get_wgd_trees_status(segment,wgd_status):
    return 'Default'

def add_wgd_info_to_table(timing_table,wgd_info,wgd_status):
    timing_table = timing_table.copy()
    timing_table['WGD_Timing'] = wgd_info['Timing']
    timing_table['WGD_Timing_CI_Low'] = wgd_info['CI_Low']
    timing_table['WGD_Timing_CI_High'] = wgd_info['CI_High']
    timing_table['WGD_Status'] = wgd_status
    return timing_table 

def write_timing_tables(timing_table,timing_table_path):
    timing_table = timing_table[['Sample_ID','Segment_ID', 'Chromosome', 'Segment_Start', 'Segment_End', 'Major_CN',
       'Minor_CN', 'Total_CN', 'N_Mutations','Mutation_Rate','Route','Probability','Average_N_Events','Average_Pre_WGD_Losses','Average_Post_WGD_Losses','Time','Node','Node_Phasing','Timing', 'Timing_CI_Low', 'Timing_CI_High','WGD_Status','WGD_Timing','WGD_Timing_CI_Low',
       'WGD_Timing_CI_High','Density']]

    if os.path.exists(timing_table_path):
        timing_table_prev = pd.read_csv(timing_table_path,sep='\t')
        timing_table = pd.concat([timing_table_prev,timing_table])
    timing_table.to_csv(timing_table_path,sep="\t",index=False)

def write_timing_dict(timing_dict,dict_dir,segment_id):
    output_path = f'{dict_dir}/{segment_id.replace(":","")}_timing_dict.bz2'
    ofile = bz2.BZ2File(output_path,'wb')
    pickle.dump(timing_dict,ofile)
    ofile.close()

def get_wgd_info(wgd_timing_distribution):
    if wgd_timing_distribution is None:
        wgd_timing = np.nan
        wgd_timing_ci_high = np.nan
        wgd_timing_ci_low = np.nan
    else:
        wgd_timing = np.round(np.percentile(wgd_timing_distribution,50),3)
        wgd_timing_ci_high = np.round(np.percentile(wgd_timing_distribution,95),3)
        wgd_timing_ci_low= np.round(np.percentile(wgd_timing_distribution,5),3)
    return {'Timing':wgd_timing,'CI_High':wgd_timing_ci_high,'CI_Low':wgd_timing_ci_low}


def get_potential_wgd_segments(sample,cn_high):
    if cn_high:
        major_cns = [3,4]
    else:
        major_cns = [2]
    
    potential_wgd_segments = []
    for segment in sample.segments:
        if segment.major_cn in major_cns and segment.n_mutations >= 10 and segment.chromosome in set(map(str,range(1,23))):
            if (segment.major_cn>= segment.minor_cn and not cn_high) or (segment.major_cn>segment.minor_cn and cn_high):
                potential_wgd_segments.append(segment)
    return potential_wgd_segments

def get_combined_distribution(distributions,n_samples=500):
    bins = np.linspace(0,1,201)
    bin_mid_points = (bins[1:]+bins[:(bins.size-1)])/2
    binned_distributions = []
    for distribution in distributions:
        distribution = np.clip(distribution,1e-8,1-1e-8)
        binned_distribution = np.histogram(distribution,bins=bins)[0]
        binned_distribution = binned_distribution/np.sum(binned_distribution)
        binned_distributions.append(binned_distribution)
    binned_distributions = np.array(binned_distributions)
    binned_distributions = np.log(binned_distributions+1e-300)
    
    combined_distribution = np.sum(binned_distributions,axis=0)

    combined_distribution = np.exp(combined_distribution-logsumexp(combined_distribution))
   
    combined_distribution_samples = np.random.choice(bin_mid_points,p=combined_distribution,size=n_samples,replace=True)
    return combined_distribution_samples

def time_wgd_major_cn_2(sample,output_dir,mult_store_dir,timing_dict_dir):

    wgd_timing_table_path = output_dir/f"{sample.sample_id}_gain_timing_table_wgd_segments.tsv"
    potential_wgd_segments = get_potential_wgd_segments(sample,cn_high=False)
 
    segment_ci_store = {}
    segment_width_store = {}
    
    wgd_timing_tables = []
    for segment in potential_wgd_segments:
        print('Timing WGD Segment -',segment)
        #treat minor cn of 2 as 0 as they are equivalent
        pseudo_minor_cn = 0 if segment.minor_cn == 2 else segment.minor_cn

        classifier =  RouteClassifier(segment.major_cn,pseudo_minor_cn,False,'No_WGD',mult_store_dir)
 
        mult_probabilities =segment.multiplicity_probabilities
        #changing minor_cn here here as well
        #change this in a future version
        mult_probabilities.minor_cn = pseudo_minor_cn
        classifier.fit_routes(mult_probabilities,segment.subclone_table,segment.n_snvs,None,segment.phased)
        mult_probabilities.minor_cn = segment.minor_cn

        wgd_timing_table = classifier.get_timing_table()
        segment_dict = segment.get_info_dict()

        for key,val in segment_dict.items():
            wgd_timing_table[key] = val
        wgd_timing_tables.append(wgd_timing_table)
        
        
        classifier_route = list(classifier.routes.values())[0]
        segment_timing = classifier_route.get_node_timing(0)
        if np.isnan(segment_timing).any():
            continue
        segment_timing_ci_low = np.percentile(segment_timing,10/2)
        segment_timing_ci_high = np.percentile(segment_timing,100-10/2)
        
        segment_ci_store[segment.segment_id] = (segment_timing_ci_low,segment_timing_ci_high)
        segment_width_store[segment.segment_id] = segment.width
    
    wgd_timing_table = pd.concat(wgd_timing_tables)
    wgd_timing_table = wgd_timing_table.drop(columns=['Node','Average_Pre_WGD_Losses','Average_Post_WGD_Losses','Probability'])
 
    max_segment_width = sum(segment_width_store.values())
    segments_with_max_overlap,overlap_width,best_overlap_timing = distributiontools.get_ids_with_maximum_overlap(segment_ci_store,segment_width_store)
    
    overlap_proportion = overlap_width/max_segment_width
    
    wgd_timing_table['Best_Overlap_Timing'] = best_overlap_timing
    wgd_timing_table['Overlap_Proportion'] = overlap_proportion
    wgd_timing_table['Intersecting'] = wgd_timing_table['Segment_ID'].isin(segments_with_max_overlap)
    
    wgd_timing_table.to_csv(wgd_timing_table_path,sep="\t",index=False)
    
   
    overlapping_segments = [segment for segment in potential_wgd_segments if segment.segment_id in segments_with_max_overlap]
    non_overlapping_segment_ids = [segment.segment_id for segment in potential_wgd_segments if segment.segment_id not in segments_with_max_overlap]
    
    cn_distributions = get_combined_segment_timing_cn_2(overlapping_segments,sample.subclone_table,sample.purity,mult_store_dir,timing_dict_dir)
    wgd_timing_distribution = get_combined_distribution(cn_distributions)
    
    return wgd_timing_distribution,non_overlapping_segment_ids,overlap_proportion,best_overlap_timing



def get_combined_segment_timing_cn_2(overlapping_segments,subclone_table,sample_purity,mult_store_dir,timing_dict_dir):

    mutation_tables = []
    for segment in overlapping_segments:
        mutation_tables.append(segment.mutation_table)
    combined_mutation_table = pd.concat(mutation_tables)
    assert len(combined_mutation_table['Major_CN'].unique()) ==1
    segment_timing_store = []
    
    for minor_cn,minor_cn_mutation_table in combined_mutation_table.groupby('Minor_CN'):
        minor_cn_mutation_table = minor_cn_mutation_table.copy()
        
        #set irrelevant parameters to merged value
        minor_cn_mutation_table['Segment_ID'] = f'Minor_CN_{minor_cn}'
        minor_cn_mutation_table['Chromosome'] = np.nan
        minor_cn_mutation_table['Segment_Start'] = np.nan
        minor_cn_mutation_table['Segment_End'] = np.nan

        pseudo_minor_cn = 0 if minor_cn == 2 else minor_cn
        #sex is None as it is autosomal segments only
        new_seg = Segment(minor_cn_mutation_table,subclone_table,sample_purity,sex=None)
        classifier =  RouteClassifier(new_seg.major_cn,pseudo_minor_cn,False,'No_WGD',mult_store_dir)
        mult_probabilities =new_seg.multiplicity_probabilities
        mult_probabilities.minor_cn = pseudo_minor_cn
        classifier.fit_routes(mult_probabilities,new_seg.subclone_table,segment.n_snvs,None,new_seg.phased)
        timing_dict= classifier.get_timing_dict()
        write_timing_dict(timing_dict,timing_dict_dir,f'WGD_minor_cn_{minor_cn}')
        mult_probabilities.minor_cn = minor_cn
        segment_timing = classifier.get_best_timing()[0]
        segment_timing_store.append(segment_timing)
    
    return segment_timing_store
def get_combined_segment_timing_cn_high(overlapping_segments,subclone_table,sample_purity,mult_store_dir):

    mutation_tables = []
    for segment in overlapping_segments:
        mutation_tables.append(segment.mutation_table)
    combined_mutation_table = pd.concat(mutation_tables)
    segment_timing_store = []
    
    for segment_id,segment_table in combined_mutation_table.groupby('Segment_ID'):
      
        seg = Segment(segment_table,subclone_table,sample_purity)
        classifier =  RouteClassifier(seg.major_cn,seg.minor_cn,False,'No_WGD',mult_store_dir)
        mult_probabilities =seg.multiplicity_probabilities
        classifier.fit_routes(mult_probabilities,seg.subclone_table,segment.n_snvs,None,seg.phased)
       
        segment_timing = classifier.get_best_timing()[0]
        
        segment_timing_store.append(segment_timing)
    
    return segment_timing_store

def check_permitted_cn_state(major_cn,minor_cn,wgd_status):
    #enforcing no more than 500 routes per cn state
    if minor_cn> major_cn:
        return False
    if minor_cn <0:
        return False
    if major_cn <=0:
        return False
    if major_cn ==1:
        return True
    if major_cn ==2:
        #don't produce timing distributions for WGD segments
        if wgd_status:
            return False
        return True
    if major_cn >=3 and major_cn <=5:
        return True
    if major_cn==6 and minor_cn<=4:
        return True
    if  major_cn ==7 and minor_cn<=3:
        return True
    if major_cn ==8 and minor_cn<=1:
        return True
    return False

def get_timeable_segments(sample,wgd_status,min_mutations=0):

    complex_segments = {}
    for segment in sample.segments:
        if check_permitted_cn_state(segment.major_cn,segment.minor_cn,wgd_status) and segment.n_mutations >= min_mutations:
            segment_cn_state = (segment.major_cn,segment.minor_cn)
            if segment_cn_state not in complex_segments:
                complex_segments[segment_cn_state] = []
            complex_segments[segment_cn_state].append(segment)

    return complex_segments


def write_wgd_calling_info(wgd_timing_distribution,overlap_proportion,best_overlap_timing,major_cn_mode,wgd_status,output_dir,sample_id):
    wgd_info = get_wgd_info(wgd_timing_distribution)
    wgd_info['Major_CN_Mode'] = major_cn_mode
    wgd_info['Overlap_Proportion'] = overlap_proportion
    wgd_info['WGD_Status'] = wgd_status
    wgd_info['Best_Overlap_Timing'] = best_overlap_timing
    wgd_info = pd.DataFrame(wgd_info,index=[0])

    wgd_info.to_csv(output_dir/f"{sample_id}_wgd_calling_info.tsv",index=False,sep="\t")

def process_segments(segments,wgd_timing_distribution,output_dir,mult_store_dir,timing_dict_dir,sample_id,wgd_status,plot_trees,non_parsimony_penalty=False):

    wgd_info = get_wgd_info(wgd_timing_distribution)

    timing_table_path = output_dir/f"{sample_id}_gain_timing_table.tsv"


    for cn_state in segments:
        cn_dir = f'{mult_store_dir}/{cn_state[0]}_{cn_state[1]}'
        for i,segment in enumerate(segments[cn_state]):
            print('Timing gained segment -',segment)
            mult_probabilities =segment.multiplicity_probabilities

            wgd_trees_status = get_wgd_trees_status(segment,wgd_status)
            classifier = RouteClassifier(segment.major_cn,segment.minor_cn,wgd_status,wgd_trees_status,mult_store_dir)
            
            classifier.fit_routes(mult_probabilities,segment.subclone_table,segment.n_snvs,wgd_timing_distribution,segment.phased,non_parsimony_penalty)
            
                
            timing_dict= classifier.get_timing_dict()
            write_timing_dict(timing_dict,timing_dict_dir,segment.segment_id)
            
            if check_record_segment_timing(segment,wgd_status,wgd_trees_status):
                timing_table = classifier.get_timing_table()
                segment_dict = segment.get_info_dict()
                for key,val in segment_dict.items():
                    timing_table[key] = val
                timing_table = add_wgd_info_to_table(timing_table,wgd_info,wgd_status)
                timing_table['Sample_ID'] = sample_id
                
                write_timing_tables(timing_table,timing_table_path)

                if plot_trees:
                    plot_output_dir = f"{output_dir}/{sample_id}_tree_plots/{segment.segment_id}".replace(":","-")
                    classifier.plot_trees(plot_output_dir,str(segment),wgd_info)
        shutil.rmtree(cn_dir)
        
def apply_non_parsimony_penalty(likelihood_store,n_events,l=2.7):
    penalty = np.exp(-l*n_events)
    likelihood_store = likelihood_store*penalty
    likelihood_store = likelihood_store/np.sum(likelihood_store)
    return likelihood_store
    
def process_sample(sample,output_dir,plot_trees=True,min_wgd_overlap=0.6,wgd_override=None,non_parsimony_penalty=False):
    
    output_dir = pathlib.Path(output_dir)
    os.makedirs(output_dir,exist_ok=True)
   
    shutil.rmtree(output_dir)

    
    os.makedirs(output_dir,exist_ok=True)

    timing_dict_dir= output_dir/f"{sample.sample_id}_timing_dicts/"
    os.makedirs(timing_dict_dir,exist_ok=True)

    sample_mutation_table = sample.get_mutation_table()
    sample_mutation_table_path = output_dir/f"{sample.sample_id}_mutation_table.tsv"
    sample_mutation_table.to_csv(sample_mutation_table_path,sep="\t",index=False)
    
    sample_subclone_table = sample.get_subclone_table()
    sample_subclone_table_path = output_dir/f"{sample.sample_id}_subclone_table.tsv"
    if sample_subclone_table is not None:
        sample_subclone_table.to_csv(sample_subclone_table_path,sep="\t",index=False)
    
    mult_store_dir = output_dir/f"{sample.sample_id}_mult_stores_temp"
 
    os.makedirs(mult_store_dir,exist_ok=True)
    major_cn_mode = get_major_cn_mode(sample)
    
    if major_cn_mode >2:
        raise ValueError(f"Major CN mode is {major_cn_mode} for sample {sample.sample_id}, this is currently not supported")

    if wgd_override is not None:
        if major_cn_mode ==1 and wgd_override:
            warnings.warn('Sample was specified as WGD but major CN mode is 1. Please check WGD status of sample')
        if major_cn_mode ==2 and not wgd_override:
            warnings.warn('Sample was specified as non-WGD but major CN mode is 2. Please check WGD status of sample')

        
        if wgd_override:
            wgd_timing_distribution,non_overlapping_segments,overlap_proportion,best_overlap_timing= time_wgd_major_cn_2(sample,output_dir,mult_store_dir,timing_dict_dir)
            wgd_status = True
        else:
            wgd_status = False
            wgd_timing_distribution = None
            non_overlapping_segments = []
            best_overlap_timing = np.nan
            overlap_proportion = np.nan
    else:
        if major_cn_mode ==1:
            wgd_status = False
            wgd_timing_distribution = None
            non_overlapping_segments = []
            best_overlap_timing = np.nan
            overlap_proportion = np.nan
        elif major_cn_mode ==2:
            wgd_timing_distribution,non_overlapping_segments,overlap_proportion,best_overlap_timing= time_wgd_major_cn_2(sample,output_dir,mult_store_dir,timing_dict_dir)
            if overlap_proportion < min_wgd_overlap:
                warnings.warn('The major CN mode is 2, but the overlap proportion is less than 0.6. There are a lot of copy number 2 segments, but they are not overlapping enough to be confident in a WGD call. The sample will be treated as a non-WGD sample. Proceed with caution in downstream analysis for this sample')

                wgd_timing_distribution = None
                non_overlapping_segments = []
                wgd_status = False
            else:
                wgd_status = True

        write_wgd_calling_info(wgd_timing_distribution,overlap_proportion,best_overlap_timing,major_cn_mode,wgd_status,output_dir,sample.sample_id)
    if major_cn_mode <=2:

        timeable_complex_segments = get_timeable_segments(sample,non_overlapping_segments,wgd_status)

    
        timing_table_path = f"{output_dir}/{sample.sample_id}_gain_timing_table.tsv"
            

        if non_parsimony_penalty and not wgd_status:
            warnings.warn('Warning: The non-parsimony penalty is only relevant to WGD samples, but no WGD has been identified. Penalty will not be applied')
            non_parsimony_penalty = False
        process_segments(timeable_complex_segments,wgd_timing_distribution,output_dir,mult_store_dir,timing_dict_dir,sample.sample_id,wgd_status,plot_trees=plot_trees,non_parsimony_penalty=non_parsimony_penalty)

        if os.path.exists(timing_table_path):
            for apply_penalty in [True,False]:
              sample_posterior_table = posteriortablegen.get_sample_posterior_table(timing_table_path,output_dir,sample.sample_id,apply_penalty=apply_penalty)
              sample_posterior_table_summary = posteriortablegen.get_sample_posterior_table_summary(sample_posterior_table)
     
        
          
              posterior_table_path = output_dir/f"{sample.sample_id}_posterior_timing_table_penalty_{apply_penalty}.tsv"
              posterior_table_summary_path = output_dir/f"{sample.sample_id}_posterior_timing_table_summary_{apply_penalty}.tsv"
  
              sample_posterior_table.to_csv(posterior_table_path,sep="\t",index=False)
              sample_posterior_table_summary.to_csv(posterior_table_summary_path,sep="\t",index=False)

    shutil.rmtree(mult_store_dir)
    


    
