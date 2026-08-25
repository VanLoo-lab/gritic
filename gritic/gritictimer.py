import os
import gzip
import hashlib
import json
import logging
import sys
from numbers import Integral, Real

import warnings

import pandas as pd
import numpy as np
import networkx as nx


from scipy.special import logsumexp

from scipy.optimize import nnls
from scipy.linalg import null_space

from gritic.sampletools import Segment, get_major_cn_mode, validate_sample_id
from gritic import dataloader
from gritic.intervaltools import (
    DEFAULT_TIMING_INTERVALS,
    TimingIntervalConfig,
    get_interval_bounds,
)
from gritic.tableschemas import (
    GAIN_TIMING_TABLE_COLUMNS,
    ROUTE_TABLE_COLUMNS,
)
from gritic.timingio import write_timing_archive

from sklearn.neighbors import NearestNeighbors

import time
import shutil

import gritic.distributiontools as distributiontools
import gritic.treetools as treetools
from gritic.route_tree import RouteTree
import gritic.hitandrun as hitandrun

import gritic.posteriortablegen as posteriortablegen

import pathlib


logger = logging.getLogger(__name__)


MIN_WGD_MUTATIONS = 10
ROUTE_CONDITIONAL_SAMPLE_COUNT = 1000
COMBINED_WGD_SAMPLE_COUNT = 500
SUBCLONE_PRIOR_MODES = ('corrected', 'uncorrected')
DEFAULT_SUBCLONE_PRIOR = 'corrected'


def validate_subclone_prior(subclone_prior):
    if subclone_prior not in SUBCLONE_PRIOR_MODES:
        permitted = ', '.join(SUBCLONE_PRIOR_MODES)
        raise ValueError(
            f'subclone_prior must be one of: {permitted}'
        )
    return subclone_prior


def get_sample_clone_fractions(subclone_table):
    if subclone_table is None:
        return np.array([1.0])

    subclone_fractions = np.asarray(
        subclone_table['Subclone_Fraction'],
        dtype=float,
    )
    return np.concatenate((
        [1.0 - np.sum(subclone_fractions)],
        subclone_fractions,
    ))


def get_clone_share_prior_alpha(
    sample_clone_fractions,
    subclonal_correction_array=None,
    subclone_prior=DEFAULT_SUBCLONE_PRIOR,
):
    subclone_prior = validate_subclone_prior(subclone_prior)
    clone_fractions = np.asarray(sample_clone_fractions, dtype=float)
    if clone_fractions.ndim != 1 or clone_fractions.size == 0:
        raise ValueError(
            'sample_clone_fractions must be a non-empty one-dimensional array'
        )
    if (
        not np.isfinite(clone_fractions).all()
        or (clone_fractions < 0).any()
        or not np.isclose(np.sum(clone_fractions), 1.0)
    ):
        raise ValueError(
            'sample_clone_fractions must contain finite non-negative values '
            'that sum to 1'
        )

    if subclone_prior == 'uncorrected':
        prior_fractions = clone_fractions
    else:
        correction = np.asarray(subclonal_correction_array, dtype=float)
        if correction.shape != clone_fractions.shape:
            raise ValueError(
                'subclonal_correction_array must have the same shape as '
                'sample_clone_fractions'
            )
        if not np.isfinite(correction).all() or (correction <= 0).any():
            raise ValueError(
                'subclonal_correction_array must contain finite positive '
                'values'
            )
        # Normalize in log space so a very small but positive detection
        # probability does not overflow the inverse correction.
        log_prior_weights = np.full(clone_fractions.shape, -np.inf)
        positive_fraction = clone_fractions > 0
        log_prior_weights[positive_fraction] = (
            np.log(clone_fractions[positive_fraction])
            - np.log(correction[positive_fraction])
        )
        prior_fractions = np.exp(
            log_prior_weights - logsumexp(log_prior_weights)
        )

    return 1.0 + prior_fractions


def get_segment_cache_identity(segment_id):
    return f'segment:{segment_id}'


def get_pooled_wgd_cache_identity(minor_cn, source_segment_ids):
    source_segment_ids = sorted(map(str, source_segment_ids))
    context = json.dumps(
        {
            'kind': 'pooled_wgd',
            'minor_cn': _json_scalar(minor_cn),
            'source_segment_ids': source_segment_ids,
        },
        sort_keys=True,
        separators=(',', ':'),
    )
    digest = hashlib.sha256(context.encode('utf-8')).hexdigest()
    return f'pooled-wgd:{digest}'


def get_subclone_prior_cache_namespace(
    subclone_prior,
    segment_cache_identity=None,
):
    subclone_prior = validate_subclone_prior(subclone_prior)
    if subclone_prior == 'uncorrected':
        return pathlib.Path('uncorrected')
    if segment_cache_identity is None:
        raise ValueError(
            'segment_cache_identity is required for the corrected '
            'subclone prior'
        )
    identity = str(segment_cache_identity)
    if not identity:
        raise ValueError(
            'segment_cache_identity must not be empty for the corrected '
            'subclone prior'
        )
    digest = hashlib.sha256(identity.encode('utf-8')).hexdigest()
    return pathlib.Path('corrected') / digest


def get_cache_subclone_prior(subclone_prior, subclone_table):
    """Use the shared cache when no subclone prior is sampled."""
    subclone_prior = validate_subclone_prior(subclone_prior)
    if subclone_table is None:
        return 'uncorrected'
    if isinstance(subclone_table, pd.DataFrame) and subclone_table.empty:
        return 'uncorrected'
    return subclone_prior


def _remove_cache_directory_preserving_error(cache_dir):
    """Remove a cache directory without masking an active exception."""
    cache_dir = pathlib.Path(cache_dir)
    if not cache_dir.is_dir():
        return
    active_exception = sys.exc_info()[0] is not None
    try:
        shutil.rmtree(cache_dir)
    except OSError:
        if not active_exception:
            raise
        logger.exception(
            'Failed to remove cache directory %s while preserving the '
            'active processing error',
            cache_dir,
        )


def _remove_cache_file_preserving_error(cache_path):
    """Remove a cache file without masking an active exception."""
    cache_path = pathlib.Path(cache_path)
    if not cache_path.is_file():
        return
    active_exception = sys.exc_info()[0] is not None
    try:
        cache_path.unlink()
    except OSError:
        if not active_exception:
            raise
        logger.exception(
            'Failed to remove cache file %s while preserving the active '
            'processing error',
            cache_path,
        )


def remove_mult_store_cache_namespace(
    mult_store_dir,
    subclone_prior,
    segment_cache_identity=None,
):
    """Remove one resolved cache namespace, if it exists."""
    if mult_store_dir is None:
        return
    namespace = get_subclone_prior_cache_namespace(
        subclone_prior,
        segment_cache_identity,
    )
    cache_dir = pathlib.Path(mult_store_dir, namespace)
    _remove_cache_directory_preserving_error(cache_dir)

    corrected_root = pathlib.Path(mult_store_dir, 'corrected')
    if corrected_root.is_dir():
        try:
            corrected_root.rmdir()
        except OSError:
            pass


def remove_mult_store_route_caches(
    mult_store_dir,
    subclone_prior,
    route_ids,
    wgd_status,
    wgd_timing_distribution,
    segment_cache_identity=None,
):
    """Remove cached arrays for exact routes in one WGD context."""
    if mult_store_dir is None:
        return
    namespace = get_subclone_prior_cache_namespace(
        subclone_prior,
        segment_cache_identity,
    )
    cache_dir = pathlib.Path(mult_store_dir, namespace)
    cache_prefixes = (
        'mult_store',
        'timing_store',
        'wgd_timing_store',
        'density_store',
    )
    for route_id in route_ids:
        cache_key = get_mult_store_cache_key(
            route_id,
            wgd_status,
            wgd_timing_distribution,
        )
        for cache_prefix in cache_prefixes:
            _remove_cache_file_preserving_error(
                cache_dir/f'{cache_prefix}_{cache_key}.npy.gz'
            )
    if cache_dir.is_dir():
        try:
            cache_dir.rmdir()
        except OSError:
            pass


def get_wgd_context_digest(wgd_status, wgd_timing_distribution):
    """Digest the explicit WGD status and exact timing distribution."""
    if not isinstance(wgd_status, (bool, np.bool_)):
        raise ValueError('wgd_status must be a boolean')

    hasher = hashlib.sha256()
    hasher.update(f'wgd-status:{int(wgd_status)}\0'.encode('ascii'))
    if wgd_timing_distribution is None:
        hasher.update(b'no-wgd-timing-distribution')
        return hasher.hexdigest()

    distribution = np.ascontiguousarray(
        np.asarray(wgd_timing_distribution, dtype='<f8')
    )
    hasher.update(b'wgd-timing-distribution\0')
    hasher.update(json.dumps(distribution.shape).encode('ascii'))
    hasher.update(b'\0')
    hasher.update(distribution.tobytes())
    return hasher.hexdigest()


def get_mult_store_cache_key(
    route_id,
    wgd_status,
    wgd_timing_distribution,
):
    """Identify a full route within its explicit WGD sampling context."""
    return (
        f'{route_id}_wgd_'
        f'{get_wgd_context_digest(wgd_status, wgd_timing_distribution)}'
    )


class Route:
    def __init__(
        self,
        route_id,
        tree,
        major_cn,
        minor_cn,
        wgd_status,
        mult_store_dir,
        cache_namespace,
    ):
        self.route_id = route_id
        self.short_id = route_id[:9]
        self.route_tree = RouteTree(tree,major_cn,minor_cn,wgd_status)
        
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.total_cn = self.major_cn + self.minor_cn

        self.wgd_status = wgd_status
        self.mult_store_dir = mult_store_dir
        self.cache_namespace = cache_namespace
        self.ll_store = None
        self.node_timing = None
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
        
    
    def get_weighted_arrays(
        self,
        cumulative_timing,
        wgd_timing_store,
        mult_store,
        weights,
        n_samples=ROUTE_CONDITIONAL_SAMPLE_COUNT,
    ):
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
        
        start_sol = nnls(constraints_matrix, constraints_sum)[0]
        constraints_null = null_space(constraints_matrix)
        
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
        for _ in range(n_samples):
            random_index = np.random.choice(node_timing.shape[1])
            n_events,pre_wgd_losses,post_wgd_losses = self.route_tree.get_n_events(node_timing[:,random_index],wgd_timing[random_index])
            n_events_estimate['N_Events'].append(n_events)
            n_events_estimate['Pre_WGD_Losses'].append(pre_wgd_losses)
            n_events_estimate['Post_WGD_Losses'].append(post_wgd_losses)

        return n_events_estimate 
    

    @staticmethod
    def simulate_clone_share(alpha,n_samples):
        dirichlet_sample = np.random.dirichlet(alpha,size=n_samples)
        return dirichlet_sample
    def get_density_estimate(self,samples,n_test_points=1000,radius=0.05):
        nn_finder = NearestNeighbors(radius=radius,p=1)
        nn_finder.fit(samples)
        random_mult_indicies = np.random.choice(samples.shape[0],size=n_test_points,replace=False)
        nearest_neighbors = nn_finder.radius_neighbors(samples[random_mult_indicies,:],return_distance=False)
        nn_size = np.array([x.size-1 for x in nearest_neighbors])
        return np.mean(nn_size>0.1),np.mean(nn_size>2.1)
    
    def run_mult_sampling(self,alpha,n_subclones,wgd_timing_distribution,samples_per_run=500,max_samples=5e5,density_cut_off=0.9):
        timing_store = []
        wgd_timing_store= []
        mult_store =[]

        start_sampling_time = time.perf_counter()
        next_eval_time = None
        eval_count =0
        
        while True:
            if self.wgd_status:
                wgd_timing = np.random.choice(wgd_timing_distribution)
            else:
                wgd_timing =np.nan
            mults,timing = self.sample_mults(wgd_timing,samples_per_run)
           
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
                if density >= density_cut_off:
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
        with gzip.open(path,'wb') as f:
            np.save(f,array)
    def load_gz_numpy(self,path):
        with gzip.open(path,'rb') as f:
            return np.load(f)
    
    def get_mult_store(self,alpha,n_subclones,wgd_timing_distribution):
        if self.mult_store_dir is not None:
            cache_dir = pathlib.Path(
                self.mult_store_dir,
                self.cache_namespace,
            )
            os.makedirs(cache_dir,exist_ok=True)

            cache_key = get_mult_store_cache_key(
                self.route_id,
                self.wgd_status,
                wgd_timing_distribution,
            )
            cache_suffix = f'{cache_key}.npy.gz'
            mult_store_path = cache_dir/f'mult_store_{cache_suffix}'
            timing_store_path = cache_dir/f'timing_store_{cache_suffix}'
            wgd_timing_store_path = (
                cache_dir/f'wgd_timing_store_{cache_suffix}'
            )
            density_store_path = cache_dir/f'density_store_{cache_suffix}'

            if os.path.exists(mult_store_path):
                mult_store = self.load_gz_numpy(mult_store_path)
                timing_store = self.load_gz_numpy(timing_store_path)
                
                wgd_timing_store = self.load_gz_numpy(wgd_timing_store_path)
                density = self.load_gz_numpy(density_store_path)


                return mult_store,timing_store,wgd_timing_store,density
        mult_store,timing_store,wgd_timing_store,density = self.run_mult_sampling(alpha,n_subclones,wgd_timing_distribution)
        
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
    
    def run_sampling(
        self,
        mult_probabilities,
        subclone_table,
        wgd_timing_distribution,
        subclone_prior=DEFAULT_SUBCLONE_PRIOR,
    ):

        run_time = time.perf_counter()

        subclone_prior = validate_subclone_prior(subclone_prior)
        n_subclones = (
            0
            if subclone_table is None
            else len(subclone_table.index)
        )
        alpha = None
        if n_subclones > 0:
            sample_clone_fractions = get_sample_clone_fractions(
                subclone_table
            )
            subclonal_correction_array = None
            if subclone_prior == 'corrected':
                subclonal_correction_array = (
                    mult_probabilities.get_subclonal_correction_array(
                        subclone_table
                    )
                )
            alpha = get_clone_share_prior_alpha(
                sample_clone_fractions,
                subclonal_correction_array,
                subclone_prior,
            )

        mult_store,timing_store,wgd_timing_store,density = self.get_mult_store(alpha,n_subclones,wgd_timing_distribution)

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
    def __init__(
        self,
        major_cn,
        minor_cn,
        wgd_status,
        wgd_trees_status,
        mult_store_dir,
        *,
        subclone_prior=DEFAULT_SUBCLONE_PRIOR,
        segment_cache_identity=None,
    ):
        self.major_cn = major_cn
        self.minor_cn = minor_cn
        self.wgd_status = wgd_status
        self.mult_store_dir = mult_store_dir
        self.subclone_prior = validate_subclone_prior(subclone_prior)
        self.segment_cache_identity = segment_cache_identity
        self.cache_namespace = None
        self.routes = self.generate_routes(wgd_trees_status)

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
            route = Route(
                tree_id,
                tree,
                self.major_cn,
                self.minor_cn,
                self.wgd_status,
                self.mult_store_dir,
                self.cache_namespace,
            )
            possible_routes[tree_id] = route
        return possible_routes
 
    def fit_routes(self,mult_probabilities,subclone_table,wgd_timing_distribution):
        self.cache_subclone_prior = get_cache_subclone_prior(
            self.subclone_prior,
            subclone_table,
        )
        if getattr(self, 'mult_store_dir', None) is not None:
            self.cache_namespace = get_subclone_prior_cache_namespace(
                self.cache_subclone_prior,
                self.segment_cache_identity,
            )
            for route in self.routes.values():
                route.cache_namespace = self.cache_namespace

        route_ll_store = []
        route_ids = list(sorted(self.routes.keys()))
        for route_id in route_ids:
            route = self.routes[route_id]
            route.run_sampling(
                mult_probabilities,
                subclone_table,
                wgd_timing_distribution,
                self.subclone_prior,
            )
            route_ll_store.append(route.ll_store)

        #helps keep the exponentiation under control
        max_point = np.max(np.concatenate(route_ll_store))
        likelihood_store =[]
        
        for route_ll in route_ll_store:
            route_likelihoods = np.exp(route_ll-max_point)
            average_likelihood = np.average(route_likelihoods)
            likelihood_store.append(average_likelihood)
        
        likelihood_store = np.array(likelihood_store)/np.sum(likelihood_store)

        for i,route_id in enumerate(route_ids):
            self.route_probabilities[route_id] = likelihood_store[i]

    def _get_output_routes(self):
        if len(self.route_probabilities) == 0:
            raise ValueError("routes hasn't been fit yet")

        output_routes = []
        full_ids_by_short_id = {}
        for route_id, probability in self.route_probabilities.items():
            route = self.routes[route_id]
            output_routes.append((route, probability))
            full_ids_by_short_id.setdefault(route.short_id, []).append(route_id)

        collisions = {
            short_id: full_ids
            for short_id, full_ids in full_ids_by_short_id.items()
            if len(full_ids) > 1
        }
        if collisions:
            collision_text = '; '.join(
                f"{short_id}: {', '.join(map(str, full_ids))}"
                for short_id, full_ids in collisions.items()
            )
            raise ValueError(
                'Route identifiers must be unique after shortening to the '
                f'output Route value; collisions: {collision_text}'
            )
        return output_routes

    def get_route_table(self):
        route_table_data = {
            'Route': [],
            'Probability': [],
            'Average_N_Events': [],
            'Average_Pre_WGD_Losses': [],
            'Average_Post_WGD_Losses': [],
            'Time': [],
            'Density': [],
        }
        for route, probability in self._get_output_routes():
            route_table_data['Route'].append(route.short_id)
            route_table_data['Probability'].append(probability)
            route_table_data['Average_N_Events'].append(
                route.get_average_events('N_Events')
            )
            route_table_data['Average_Pre_WGD_Losses'].append(
                route.get_average_events('Pre_WGD_Losses')
            )
            route_table_data['Average_Post_WGD_Losses'].append(
                route.get_average_events('Post_WGD_Losses')
            )
            route_table_data['Time'].append(np.round(route.run_time, 3))
            route_table_data['Density'].append(route.density)

        return pd.DataFrame(route_table_data)

    def get_timing_table(
        self,
        interval=DEFAULT_TIMING_INTERVALS.route_gain,
        rounding_digits=3,
    ):
        timing_table_data = {
            'Route': [],
            'Node': [],
            'Node_Phasing': [],
            'Timing': [],
            'Timing_CI_Low': [],
            'Timing_CI_High': [],
        }
        output_routes = self._get_output_routes()
        for route, _ in output_routes:
            for node in route.route_tree.timeable_nodes:
                node_timing = np.asarray(route.get_node_timing(node))
                timing_ci_low, timing_ci_high = get_interval_bounds(
                    node_timing,
                    interval,
                )
                timing_table_data['Route'].append(route.short_id)
                timing_table_data['Node'].append(node)
                timing_table_data['Node_Phasing'].append(
                    route.route_tree.node_attributes[node]['Phasing']
                )
                timing_table_data['Timing'].append(
                    np.round(np.median(node_timing), rounding_digits)
                )
                timing_table_data['Timing_CI_Low'].append(
                    np.round(timing_ci_low, rounding_digits)
                )
                timing_table_data['Timing_CI_High'].append(
                    np.round(timing_ci_high, rounding_digits)
                )

        timing_table = pd.DataFrame(timing_table_data)
        duplicate_gains = timing_table.duplicated(
            subset=['Route', 'Node'],
            keep=False,
        )
        if duplicate_gains.any():
            duplicate_keys = (
                timing_table.loc[duplicate_gains, ['Route', 'Node']]
                .drop_duplicates()
                .astype(str)
                .agg(':'.join, axis=1)
                .tolist()
            )
            raise ValueError(
                'Each timeable gain must have a unique Route and Node; '
                f"duplicates: {', '.join(duplicate_keys)}"
            )

        route_ids = {route.short_id for route, _ in output_routes}
        unexpected_routes = set(timing_table['Route']) - route_ids
        if unexpected_routes:
            raise ValueError(
                'Every timing-table Route must occur in the route table; '
                f"unexpected routes: {', '.join(map(str, unexpected_routes))}"
            )
        return timing_table
    
    def get_timing_dict(self,n_samples = 5000):
        timing_dict = {}
        
        for route, _ in self._get_output_routes():
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
            timing_dict[route.short_id] = route_samples
        return timing_dict

    def get_timing_tree_labels(
        self,
        route,
        wgd_info,
        gain_interval=DEFAULT_TIMING_INTERVALS.tree_gain,
        rounding_digits=3,
    ):
        node_labels = {}
        for node in route.route_tree.non_phased_node_order:
            if node in route.route_tree.timeable_nodes:
                timing_dist = np.asarray(route.get_node_timing(node))
                timing = np.round(np.median(timing_dist), rounding_digits)
                timing_ci_low, timing_ci_high = get_interval_bounds(
                    timing_dist,
                    gain_interval,
                )
                timing_ci_low = np.round(timing_ci_low, rounding_digits)
                timing_ci_high = np.round(timing_ci_high, rounding_digits)
                node_labels[node] = f"{timing} - [{timing_ci_low},{timing_ci_high}]"
            elif node in route.route_tree.wgd_nodes:
                node_labels[node] = (
                    f"{wgd_info['WGD_Timing']} - "
                    f"[{wgd_info['WGD_Timing_CI_Low']},"
                    f"{wgd_info['WGD_Timing_CI_High']}]"
                )
            else:
                node_labels[node] = ''
        return node_labels
    
    def plot_trees(
        self,
        plot_output_dir,
        seg_title,
        wgd_info,
        gain_interval=DEFAULT_TIMING_INTERVALS.tree_gain,
    ):
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

            node_labels = self.get_timing_tree_labels(
                route,
                wgd_info,
                gain_interval=gain_interval,
            )
            plotting_tree = route.route_tree.main_tree.copy()
            nx.set_node_attributes(plotting_tree,node_labels,'Label')
            
            treetools.plot_tree(plotting_tree,plot_title,output_path=route_output_path)
def add_wgd_info_to_route_table(route_table,wgd_info,wgd_status:bool):
    route_table = route_table.copy()
    route_table['WGD_Timing'] = wgd_info['WGD_Timing']
    route_table['WGD_Timing_CI_Low'] = wgd_info['WGD_Timing_CI_Low']
    route_table['WGD_Timing_CI_High'] = wgd_info['WGD_Timing_CI_High']
    route_table['WGD_Status'] = wgd_status
    return route_table


def _append_table(table, table_path, columns):
    table = table[columns]
    if os.path.exists(table_path):
        string_converters = {
            column: str
            for column in ('Sample_ID', 'Segment_ID', 'Chromosome', 'Route')
            if column in columns
        }
        previous_table = pd.read_csv(
            table_path,
            sep='\t',
            converters=string_converters,
        )
        table = pd.concat([previous_table, table], ignore_index=True)
    _write_tsv(table, table_path)


def _write_tsv(table, table_path):
    """Write a public TSV."""
    table.to_csv(table_path, sep='\t', index=False)


def write_route_table(route_table, route_table_path):
    _append_table(route_table, route_table_path, ROUTE_TABLE_COLUMNS)


def write_gain_timing_table(timing_table, timing_table_path):
    _append_table(
        timing_table,
        timing_table_path,
        GAIN_TIMING_TABLE_COLUMNS,
    )


def _initialize_table(table_path, columns):
    if not os.path.exists(table_path):
        _write_tsv(
            pd.DataFrame(columns=columns),
            table_path,
        )

def get_wgd_info(
    wgd_timing_distribution,
    interval=DEFAULT_TIMING_INTERVALS.sample_wgd,
    rounding_digits=3,
):
    if wgd_timing_distribution is None:
        wgd_timing = np.nan
        wgd_timing_ci_high = np.nan
        wgd_timing_ci_low = np.nan
    else:
        wgd_timing = np.round(
            np.median(wgd_timing_distribution),
            rounding_digits,
        )
        wgd_timing_ci_low, wgd_timing_ci_high = get_interval_bounds(
            wgd_timing_distribution,
            interval,
        )
        wgd_timing_ci_low = np.round(
            wgd_timing_ci_low,
            rounding_digits,
        )
        wgd_timing_ci_high = np.round(
            wgd_timing_ci_high,
            rounding_digits,
        )
    return {
        'WGD_Timing': wgd_timing,
        'WGD_Timing_CI_Low': wgd_timing_ci_low,
        'WGD_Timing_CI_High': wgd_timing_ci_high,
    }


def get_potential_wgd_segments(
    sample,
    min_wgd_mutations: int = MIN_WGD_MUTATIONS,
):
    valid_autosomes = set(sample.autosomes)
    potential_wgd_segments = []
    for segment in sample.segments:
        if segment.major_cn ==2 and segment.n_mutations >= min_wgd_mutations and segment.chromosome in valid_autosomes:
            potential_wgd_segments.append(segment)
    return potential_wgd_segments

def get_combined_distribution(
    distributions,
    n_samples=COMBINED_WGD_SAMPLE_COUNT,
    eps=1e-300,
):
    bins = np.linspace(0,1,201)
    bin_mid_points = (bins[1:]+bins[:(bins.size-1)])/2
    binned_distributions = []
    for distribution in distributions:
        distribution = np.clip(distribution,eps,1.0-eps)
        binned_distribution = np.histogram(distribution,bins=bins)[0]
        binned_distribution = binned_distribution/np.sum(binned_distribution)
        binned_distributions.append(binned_distribution)
    binned_distributions = np.array(binned_distributions)
    binned_distributions = np.log(binned_distributions+eps)
    
    combined_distribution = np.sum(binned_distributions,axis=0)

    combined_distribution = np.exp(combined_distribution-logsumexp(combined_distribution))
   
    combined_distribution_samples = np.random.choice(bin_mid_points,p=combined_distribution,size=n_samples,replace=True)
    return combined_distribution_samples


def _time_wgd_segment(
    segment,
    mult_store_dir,
    interval_config,
    subclone_prior,
):
    segment_cache_identity = get_segment_cache_identity(
        getattr(segment, 'segment_id', str(segment))
    )
    subclone_table = getattr(segment, 'subclone_table', None)
    cache_subclone_prior = get_cache_subclone_prior(
        subclone_prior,
        subclone_table,
    )
    try:
        pseudo_minor_cn = 0 if segment.minor_cn == 2 else segment.minor_cn
        classifier = RouteClassifier(
            segment.major_cn,
            pseudo_minor_cn,
            False,
            'No_WGD',
            mult_store_dir,
            subclone_prior=subclone_prior,
            segment_cache_identity=segment_cache_identity,
        )

        mult_probabilities = segment.multiplicity_probabilities
        original_minor_cn = mult_probabilities.minor_cn
        try:
            mult_probabilities.minor_cn = pseudo_minor_cn
            classifier.fit_routes(
                mult_probabilities,
                subclone_table,
                None,
            )
        finally:
            mult_probabilities.minor_cn = original_minor_cn

        wgd_route_table = classifier.get_route_table()
        wgd_timing_table = classifier.get_timing_table(
            interval=interval_config.route_gain,
        )
        wgd_timing_table = pd.merge(
            wgd_route_table,
            wgd_timing_table,
            on='Route',
            how='inner',
            validate='one_to_many',
        )
        for key, val in segment.get_info_dict().items():
            wgd_timing_table[key] = val

        classifier_route = list(classifier.routes.values())[0]
        segment_timing = classifier_route.get_node_timing(0)
        return wgd_timing_table, segment_timing
    finally:
        if cache_subclone_prior == 'corrected':
            remove_mult_store_cache_namespace(
                mult_store_dir,
                cache_subclone_prior,
                segment_cache_identity,
            )


def time_wgd_major_cn_2(
    sample,
    output_dir,
    mult_store_dir,
    timing_dict_dir,
    interval_config=DEFAULT_TIMING_INTERVALS,
    subclone_prior=DEFAULT_SUBCLONE_PRIOR,
):
    subclone_prior = validate_subclone_prior(subclone_prior)

    wgd_timing_table_path = output_dir/f"{sample.sample_id}_gain_timing_table_wgd_segments.tsv"
    potential_wgd_segments = get_potential_wgd_segments(sample)
    if not potential_wgd_segments:
        raise ValueError(
            "No eligible segments with enough mutations for WGD timing in "
            f"sample {sample.sample_id}: at least one autosomal major-copy-"
            f"number-2 segment with {MIN_WGD_MUTATIONS} or more retained "
            "mutations is required."
        )
 
    segment_ci_store = {}
    segment_width_store = {}
    
    wgd_timing_tables = []
    for segment in potential_wgd_segments:
        logger.info('Timing WGD segment: %s', segment)
        wgd_timing_table, segment_timing = _time_wgd_segment(
            segment,
            mult_store_dir,
            interval_config,
            subclone_prior,
        )
        wgd_timing_tables.append(wgd_timing_table)

        if not np.isfinite(segment_timing).all():
            continue
        segment_timing_ci_low, segment_timing_ci_high = (
            get_interval_bounds(
                segment_timing,
                interval_config.wgd_overlap,
            )
        )
        segment_ci_store[segment.segment_id] = (
            segment_timing_ci_low,
            segment_timing_ci_high,
        )
        segment_width_store[segment.segment_id] = segment.width

    if not segment_ci_store:
        raise ValueError(
            "No eligible major-copy-number-2 segment produced a finite WGD "
            f"timing interval for sample {sample.sample_id}."
        )
    
    wgd_timing_table = pd.concat(wgd_timing_tables)
    wgd_timing_table = wgd_timing_table.drop(columns=['Node','Average_Pre_WGD_Losses','Average_Post_WGD_Losses','Probability'])
 
    max_segment_width = sum(segment_width_store.values())
    segments_with_max_overlap,overlap_width,best_overlap_timing = distributiontools.get_ids_with_maximum_overlap(segment_ci_store,segment_width_store)
    
    overlap_proportion = overlap_width/max_segment_width
    
    wgd_timing_table['Best_Overlap_Timing'] = best_overlap_timing
    wgd_timing_table['Overlap_Proportion'] = overlap_proportion
    wgd_timing_table['Intersecting'] = wgd_timing_table['Segment_ID'].isin(segments_with_max_overlap)
    
    _write_tsv(wgd_timing_table, wgd_timing_table_path)
    
   
    overlapping_segments = [segment for segment in potential_wgd_segments if segment.segment_id in segments_with_max_overlap]
    non_overlapping_segment_ids = [segment.segment_id for segment in potential_wgd_segments if segment.segment_id not in segments_with_max_overlap]
    
    cn_distributions = get_combined_segment_timing_cn_2(
        overlapping_segments,
        sample.subclone_table,
        sample.purity,
        mult_store_dir,
        timing_dict_dir,
        subclone_prior=subclone_prior,
    )
    wgd_timing_distribution = get_combined_distribution(cn_distributions)
    
    return wgd_timing_distribution,non_overlapping_segment_ids,overlap_proportion,best_overlap_timing


def _time_combined_wgd_segment(
    minor_cn,
    mutation_table,
    combined_width,
    source_segment_ids,
    subclone_table,
    sample_purity,
    apply_reads_correction,
    min_mutation_alt_count,
    coverage_vaf_quantile,
    mult_store_dir,
    timing_dict_dir,
    subclone_prior,
):
    mutation_table = mutation_table.copy()
    mutation_table['Segment_ID'] = f'Minor_CN_{minor_cn}'
    mutation_table['Chromosome'] = np.nan
    mutation_table['Segment_Start'] = 0
    mutation_table['Segment_End'] = combined_width

    new_seg = Segment(
        mutation_table,
        subclone_table,
        sample_purity,
        sex=None,
        apply_reads_correction=apply_reads_correction,
        min_mutation_alt_count=min_mutation_alt_count,
        coverage_vaf_quantile=coverage_vaf_quantile,
    )
    segment_cache_identity = get_pooled_wgd_cache_identity(
        minor_cn,
        source_segment_ids,
    )
    cache_subclone_prior = get_cache_subclone_prior(
        subclone_prior,
        new_seg.subclone_table,
    )
    try:
        pseudo_minor_cn = 0 if minor_cn == 2 else minor_cn
        classifier = RouteClassifier(
            new_seg.major_cn,
            pseudo_minor_cn,
            False,
            'No_WGD',
            mult_store_dir,
            subclone_prior=subclone_prior,
            segment_cache_identity=segment_cache_identity,
        )
        mult_probabilities = new_seg.multiplicity_probabilities
        original_minor_cn = mult_probabilities.minor_cn
        try:
            mult_probabilities.minor_cn = pseudo_minor_cn
            classifier.fit_routes(
                mult_probabilities,
                new_seg.subclone_table,
                None,
            )
        finally:
            mult_probabilities.minor_cn = original_minor_cn

        timing_dict = classifier.get_timing_dict()
        write_timing_archive(
            timing_dict,
            timing_dict_dir,
            f'WGD_minor_cn_{minor_cn}',
        )
        return classifier.get_best_timing()[0]
    finally:
        if cache_subclone_prior == 'corrected':
            remove_mult_store_cache_namespace(
                mult_store_dir,
                cache_subclone_prior,
                segment_cache_identity,
            )



def get_combined_segment_timing_cn_2(
    overlapping_segments,
    subclone_table,
    sample_purity,
    mult_store_dir,
    timing_dict_dir,
    subclone_prior=DEFAULT_SUBCLONE_PRIOR,
):
    subclone_prior = validate_subclone_prior(subclone_prior)

    mutation_tables = []
    combined_width_by_minor_cn = {}
    source_segment_ids_by_minor_cn = {}
    min_mutation_alt_counts = {
        segment.min_mutation_alt_count
        for segment in overlapping_segments
    }
    if len(min_mutation_alt_counts) != 1:
        raise ValueError(
            'Combined WGD segments must use one minimum mutation '
            'alternate-read count'
        )
    min_mutation_alt_count = min_mutation_alt_counts.pop()
    coverage_vaf_quantiles = {
        segment.coverage_vaf_quantile
        for segment in overlapping_segments
    }
    if len(coverage_vaf_quantiles) != 1:
        raise ValueError(
            'Combined WGD segments must use one coverage VAF quantile'
        )
    coverage_vaf_quantile = coverage_vaf_quantiles.pop()
    apply_reads_correction_values = {
        segment.apply_reads_correction
        for segment in overlapping_segments
    }
    if len(apply_reads_correction_values) != 1:
        raise ValueError(
            'Combined WGD segments must use one reads-correction setting'
        )
    apply_reads_correction = apply_reads_correction_values.pop()
    for segment in overlapping_segments:
        mutation_tables.append(segment.mutation_table)
        combined_width_by_minor_cn[segment.minor_cn] = (
            combined_width_by_minor_cn.get(segment.minor_cn, 0)
            + int(segment.width)
        )
        source_segment_ids_by_minor_cn.setdefault(
            segment.minor_cn,
            [],
        ).append(getattr(segment, 'segment_id', str(segment)))
    combined_mutation_table = pd.concat(mutation_tables)
    assert len(combined_mutation_table['Major_CN'].unique()) ==1
    segment_timing_store = []
    grouped_mutations = combined_mutation_table.groupby('Minor_CN')
    for minor_cn, minor_cn_mutation_table in grouped_mutations:
        segment_timing_store.append(
            _time_combined_wgd_segment(
                minor_cn,
                minor_cn_mutation_table,
                combined_width_by_minor_cn[minor_cn],
                source_segment_ids_by_minor_cn[minor_cn],
                subclone_table,
                sample_purity,
                apply_reads_correction,
                min_mutation_alt_count,
                coverage_vaf_quantile,
                mult_store_dir,
                timing_dict_dir,
                subclone_prior,
            )
        )
    
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
        #we don't produce timing distributions for WGD segments
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

def get_timeable_segments(
    sample,
    *,
    wgd_status,
    excluded_segment_ids=None,
    min_mutations=0,
):

    complex_segments = {}
    excluded_segment_ids = (
        set()
        if excluded_segment_ids is None
        else set(excluded_segment_ids)
    )
    for segment in sample.segments:
        if segment.segment_id in excluded_segment_ids:
            continue
        if check_permitted_cn_state(segment.major_cn,segment.minor_cn,wgd_status) and segment.n_mutations >= min_mutations:
            segment_cn_state = (segment.major_cn,segment.minor_cn)
            if segment_cn_state not in complex_segments:
                complex_segments[segment_cn_state] = []
            complex_segments[segment_cn_state].append(segment)

    return complex_segments


def _json_scalar(value):
    if isinstance(value, np.generic):
        value = value.item()
    if (
        value is None
        or pd.isna(value)
        or (isinstance(value, float) and not np.isfinite(value))
    ):
        return None
    return value


def write_wgd_calling_info(
    wgd_info,
    overlap_proportion,
    best_overlap_timing,
    major_cn_mode,
    wgd_status,
    output_dir,
    sample_id,
):
    calling_info = {
        **wgd_info,
        'Major_CN_Mode': major_cn_mode,
        'Overlap_Proportion': overlap_proportion,
        'WGD_Status': wgd_status,
        'Best_Overlap_Timing': best_overlap_timing,
    }
    calling_info = {
        key: _json_scalar(value)
        for key, value in calling_info.items()
    }
    output_path = output_dir / f'{sample_id}_wgd_calling_info.json'
    with output_path.open('w', encoding='utf-8') as output_file:
        json.dump(calling_info, output_file, indent=2, allow_nan=False)
        output_file.write('\n')


def _process_segment(
    segment,
    wgd_timing_distribution,
    output_dir,
    mult_store_dir,
    timing_dict_dir,
    sample_id,
    wgd_status,
    plot_trees,
    wgd_info,
    interval_config,
    subclone_prior,
    route_table_path,
    timing_table_path,
):
    segment_cache_identity = get_segment_cache_identity(
        getattr(segment, 'segment_id', str(segment))
    )
    subclone_table = getattr(segment, 'subclone_table', None)
    cache_subclone_prior = get_cache_subclone_prior(
        subclone_prior,
        subclone_table,
    )
    try:
        classifier = RouteClassifier(
            segment.major_cn,
            segment.minor_cn,
            wgd_status,
            'Default',
            mult_store_dir,
            subclone_prior=subclone_prior,
            segment_cache_identity=segment_cache_identity,
        )
        classifier.fit_routes(
            segment.multiplicity_probabilities,
            subclone_table,
            wgd_timing_distribution,
        )

        segment_route_table = classifier.get_route_table()
        segment_timing_table = classifier.get_timing_table(
            interval=interval_config.route_gain,
        )
        timing_dict = classifier.get_timing_dict()
        write_timing_archive(
            timing_dict,
            timing_dict_dir,
            segment.segment_id,
        )

        for key, val in segment.get_info_dict().items():
            segment_route_table[key] = val
        segment_route_table = add_wgd_info_to_route_table(
            segment_route_table,
            wgd_info,
            wgd_status,
        )
        segment_route_table['Sample_ID'] = sample_id
        segment_route_table = posteriortablegen.add_penalized_probability(
            segment_route_table
        )

        segment_timing_table['Segment_ID'] = segment.segment_id
        segment_timing_table['Sample_ID'] = sample_id

        write_route_table(segment_route_table, route_table_path)
        write_gain_timing_table(segment_timing_table, timing_table_path)

        if plot_trees:
            plot_output_dir = (
                f'{output_dir}/{sample_id}_tree_plots/{segment.segment_id}'
            )
            classifier.plot_trees(
                plot_output_dir,
                str(segment),
                wgd_info,
                gain_interval=interval_config.tree_gain,
            )
    finally:
        if cache_subclone_prior == 'corrected':
            remove_mult_store_cache_namespace(
                mult_store_dir,
                cache_subclone_prior,
                segment_cache_identity,
            )


def process_segments(
    segments,
    wgd_timing_distribution,
    output_dir,
    mult_store_dir,
    timing_dict_dir,
    sample_id,
    wgd_status,
    plot_trees,
    wgd_info,
    interval_config=DEFAULT_TIMING_INTERVALS,
    subclone_prior=DEFAULT_SUBCLONE_PRIOR,
):
    subclone_prior = validate_subclone_prior(subclone_prior)

    route_table_path = output_dir/f"{sample_id}_route_table.tsv"
    timing_table_path = output_dir/f"{sample_id}_gain_timing_table.tsv"

    _initialize_table(route_table_path, ROUTE_TABLE_COLUMNS)
    _initialize_table(timing_table_path, GAIN_TIMING_TABLE_COLUMNS)

    for cn_state in segments:
        route_ids = tuple(treetools.get_nx_trees(
            cn_state[0],
            cn_state[1],
            wgd_status,
            'Default',
        ))
        try:
            for segment in segments[cn_state]:
                logger.info('Timing gained segment: %s', segment)
                _process_segment(
                    segment,
                    wgd_timing_distribution,
                    output_dir,
                    mult_store_dir,
                    timing_dict_dir,
                    sample_id,
                    wgd_status,
                    plot_trees,
                    wgd_info,
                    interval_config,
                    subclone_prior,
                    route_table_path,
                    timing_table_path,
                )
        finally:
            remove_mult_store_route_caches(
                mult_store_dir,
                'uncorrected',
                route_ids,
                wgd_status,
                wgd_timing_distribution,
            )

def _validate_wgd_count(wgd_count):
    if wgd_count is None:
        return None
    if isinstance(wgd_count, bool) or not isinstance(wgd_count, Integral):
        raise ValueError('wgd_count must be None or an integer equal to 0 or 1')
    if wgd_count not in (0, 1):
        raise ValueError('wgd_count must be None or an integer equal to 0 or 1')
    return int(wgd_count)


def _validate_min_wgd_overlap(min_wgd_overlap):
    if (
        isinstance(min_wgd_overlap, (bool, np.bool_))
        or not isinstance(min_wgd_overlap, Real)
        or not np.isfinite(min_wgd_overlap)
        or not 0 <= min_wgd_overlap <= 1
    ):
        raise ValueError(
            'min_wgd_overlap must be a finite number between 0 and 1'
        )
    return float(min_wgd_overlap)


def _process_sample_with_mult_store(
    sample,
    output_dir,
    timing_dict_dir,
    mult_store_dir,
    plot_trees,
    min_wgd_overlap,
    wgd_count,
    interval_config,
    subclone_prior,
    major_cn_mode,
):
    if wgd_count is not None:
        if major_cn_mode == 1 and wgd_count == 1:
            warnings.warn(
                'Sample was specified with WGD count 1 but major CN mode is '
                '1. Please check the WGD count of the sample'
            )
        if major_cn_mode == 2 and wgd_count == 0:
            warnings.warn(
                'Sample was specified with WGD count 0 but major CN mode is '
                '2. Please check the WGD count of the sample'
            )

        if wgd_count == 1:
            (
                wgd_timing_distribution,
                non_overlapping_segment_ids,
                overlap_proportion,
                best_overlap_timing,
            ) = time_wgd_major_cn_2(
                sample,
                output_dir,
                mult_store_dir,
                timing_dict_dir,
                interval_config=interval_config,
                subclone_prior=subclone_prior,
            )
            wgd_status = True
        else:
            wgd_status = False
            wgd_timing_distribution = None
            non_overlapping_segment_ids = []
            best_overlap_timing = np.nan
            overlap_proportion = np.nan
    elif major_cn_mode == 1:
        wgd_status = False
        wgd_timing_distribution = None
        non_overlapping_segment_ids = []
        best_overlap_timing = np.nan
        overlap_proportion = np.nan
    else:
        (
            wgd_timing_distribution,
            non_overlapping_segment_ids,
            overlap_proportion,
            best_overlap_timing,
        ) = time_wgd_major_cn_2(
            sample,
            output_dir,
            mult_store_dir,
            timing_dict_dir,
            interval_config=interval_config,
            subclone_prior=subclone_prior,
        )
        if overlap_proportion < min_wgd_overlap:
            warnings.warn(
                'The major CN mode is 2, but the overlap proportion is '
                f'less than {min_wgd_overlap}. There are a lot of copy '
                'number 2 segments, but they are not overlapping enough '
                'to be confident in a WGD call. The sample will be '
                'treated as a non-WGD sample. Proceed with caution in '
                'downstream analysis for this sample'
            )
            wgd_timing_distribution = None
            non_overlapping_segment_ids = []
            wgd_status = False
        else:
            wgd_status = True

    wgd_info = get_wgd_info(
        wgd_timing_distribution,
        interval=interval_config.sample_wgd,
    )
    write_wgd_calling_info(
        wgd_info,
        overlap_proportion,
        best_overlap_timing,
        major_cn_mode,
        wgd_status,
        output_dir,
        sample.sample_id,
    )
    timeable_complex_segments = get_timeable_segments(
        sample,
        wgd_status=wgd_status,
        excluded_segment_ids=non_overlapping_segment_ids,
    )

    route_table_path = output_dir/f'{sample.sample_id}_route_table.tsv'
    timing_table_path = output_dir/f'{sample.sample_id}_gain_timing_table.tsv'
    process_segments(
        timeable_complex_segments,
        wgd_timing_distribution,
        output_dir,
        mult_store_dir,
        timing_dict_dir,
        sample.sample_id,
        wgd_status,
        plot_trees=plot_trees,
        wgd_info=wgd_info,
        interval_config=interval_config,
        subclone_prior=subclone_prior,
    )

    if os.path.exists(route_table_path):
        for apply_penalty in (False, True):
            (
                sample_posterior_table,
                sample_posterior_route_draw_table,
            ) = posteriortablegen.get_sample_posterior_tables(
                route_table_path,
                timing_table_path,
                output_dir,
                sample.sample_id,
                apply_penalty=apply_penalty,
            )
            sample_posterior_table_summary = (
                posteriortablegen.get_sample_posterior_table_summary(
                    sample_posterior_table,
                    sample_posterior_route_draw_table,
                    interval=interval_config.posterior_summary,
                )
            )
            posterior_table_summary_path = (
                output_dir
                / f'{sample.sample_id}_posterior_timing_table_summary_'
                f'penalty_{apply_penalty}.tsv'
            )
            _write_tsv(
                sample_posterior_table_summary,
                posterior_table_summary_path,
            )


def process_sample(
    sample,
    output_dir,
    plot_trees=False,
    min_wgd_overlap=0.6,
    wgd_count=None,
    interval_config=DEFAULT_TIMING_INTERVALS,
    subclone_prior=DEFAULT_SUBCLONE_PRIOR,
):
    if not isinstance(interval_config, TimingIntervalConfig):
        raise TypeError('interval_config must be a TimingIntervalConfig')
    subclone_prior = validate_subclone_prior(subclone_prior)
    wgd_count = _validate_wgd_count(wgd_count)
    min_wgd_overlap = _validate_min_wgd_overlap(min_wgd_overlap)
    validate_sample_id(sample.sample_id)
    major_cn_mode = get_major_cn_mode(sample)
    if major_cn_mode not in (1, 2):
        raise ValueError(
            f"Major CN mode is {major_cn_mode} for sample {sample.sample_id}; "
            "GRITIC currently supports only modal major copy numbers 1 and 2"
        )

    output_dir = pathlib.Path(output_dir, pathlib.Path(sample.sample_id))
    if output_dir.exists():
        if not output_dir.is_dir():
            raise FileExistsError(
                f'Sample output path is not a directory: {output_dir}'
            )
        if any(output_dir.iterdir()):
            raise FileExistsError(
                'Sample output directory must be absent or empty before a '
                f'GRITIC run: {output_dir}'
            )
    else:
        os.makedirs(output_dir, exist_ok=False)


    timing_dict_dir= output_dir/f"{sample.sample_id}_timing_dicts/"
    os.makedirs(timing_dict_dir,exist_ok=True)

    sample_mutation_table = sample.get_mutation_table()
    sample_mutation_table_path = output_dir/f"{sample.sample_id}_mutation_table.tsv"
    _write_tsv(sample_mutation_table, sample_mutation_table_path)

    sample_subclone_table = sample.get_subclone_table()
    sample_subclone_table_path = output_dir/f"{sample.sample_id}_subclone_table.tsv"
    if sample_subclone_table is None:
        sample_subclone_table = pd.DataFrame(
            columns=dataloader.SUBCLONE_OUTPUT_COLUMNS,
        )
    _write_tsv(sample_subclone_table, sample_subclone_table_path)

    mult_store_dir = output_dir/f"{sample.sample_id}_mult_stores_temp"

    os.makedirs(mult_store_dir,exist_ok=True)
    try:
        _process_sample_with_mult_store(
            sample,
            output_dir,
            timing_dict_dir,
            mult_store_dir,
            plot_trees,
            min_wgd_overlap,
            wgd_count,
            interval_config,
            subclone_prior,
            major_cn_mode,
        )
    finally:
        _remove_cache_directory_preserving_error(mult_store_dir)
