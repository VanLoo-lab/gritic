import numpy as np
from scipy.linalg import null_space
from numba import njit

@njit
def get_random_direction(dimension):
    x = np.random.normal(0,1,size=dimension)
    return x/np.linalg.norm(x)

#enforcing non-negativity limits
@njit
def get_direction_limit(direction,current_position):
    negative_direction = direction[direction<-1e-10]
    negative_position = current_position[direction<-1e-10]
    lower_limit = -np.max(np.divide(negative_position,negative_direction))
    return lower_limit
@njit
def get_direction_range(direction,current_position):
    positive_limit = get_direction_limit(direction,current_position)
    negative_limit = -get_direction_limit(-direction,current_position)
    return negative_limit,positive_limit
@njit
def get_new_position(current_position,null_dimension,A_null):
    direction_components = get_random_direction(null_dimension)
    direction = np.sum(np.multiply(direction_components,A_null),axis=1)
    valid_direction_length_range = get_direction_range(direction,current_position)
    valid_direction_length = np.random.uniform(valid_direction_length_range[0],valid_direction_length_range[1])
    new_position = current_position+direction*valid_direction_length

    assert not (new_position<-1e-10).any()
    assert not (new_position>1+1e-10).any()

    new_position = np.clip(new_position,0.0,1.0)
    return new_position

@njit
def run_chain(current_position,null_dimension,A_null,timing_state,burn_in=25,skips=5,n_samples=100):
    n_samples_actual = skips*n_samples+burn_in
    hit_and_run_store = np.zeros((n_samples_actual,timing_state.size))
    
    for i in range(n_samples_actual):
        new_position = get_new_position(current_position,null_dimension,A_null)
        
        hit_and_run_store[i,:] = new_position
        current_position = new_position
    return hit_and_run_store[burn_in::skips]



@njit()
def hit_and_run(A_null,timing_state,n_samples=500,burn_in=25,skips=5):
 
    if A_null.shape[1] ==0:
        hit_and_run_store = np.zeros((n_samples,timing_state.size))
        for i in range(n_samples):
            hit_and_run_store[i,:] = timing_state
        return hit_and_run_store
    #note that this can cause problems hitting the 1 limit
    #but I think this is ok
    if A_null.shape[1] ==1:

        
        direction = A_null[:,0]
        valid_length_range = get_direction_range(direction,timing_state)
        valid_length_samples = np.linspace(valid_length_range[0],valid_length_range[1],n_samples)
        np.random.shuffle(valid_length_samples)
        extra_vectors = np.outer(valid_length_samples,direction)
        hit_and_run_store = extra_vectors+timing_state

        return hit_and_run_store
    
    null_dimension = A_null.shape[1]
    first_position = timing_state
    current_position = timing_state
    
    process_count =0
    while np.sum(np.abs(current_position-first_position))<1e-2:
        current_position = get_new_position(current_position,null_dimension,A_null)
        process_count+=1
        if process_count>50000:
            
            print("Warning: process count exceeded 50000")
            print('final pos',current_position)
            print('initial pos',first_position)
            print('='*50)
            break

    return run_chain(current_position,null_dimension,A_null,timing_state,burn_in=burn_in,skips=skips,n_samples=n_samples)
