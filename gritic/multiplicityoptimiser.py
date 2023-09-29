import numpy as np
from scipy.optimize import minimize,LinearConstraint

def unconstrained_mult_optimisation(mult_probabilities,n_subclones):
    n_dim = mult_probabilities.major_cn + n_subclones
    start_point = get_point_on_simplex(n_dim)
    #start_point = np.ones(n_dim)/n_dim
    bounds = [(0,1)]*n_dim
    #force a sum to one
    constraints = LinearConstraint(np.ones((1,n_dim)),1,1)
    sol = minimize(unconstrained_mult_likelihood,start_point, bounds=bounds,args=(n_subclones,mult_probabilities),constraints=constraints)
    if not sol.success:
        return None
    return sol.x
def unconstrained_mult_likelihood(x,n_subclones,mult_probabilities):
    
    likelihood = mult_probabilities.evaluate_likelihood(x)
    return -likelihood

def mult_likelihood(x,n_subclones,mult_probabilities,timing_matrix):
    timing_state =x[:(x.size-n_subclones)]
    subclonal_mult_state = x[(x.size-n_subclones):]
    
    clonal_mult_state = np.dot(timing_state,timing_matrix)

    likelihood = mult_probabilities.evaluate_likelihood(clonal_mult_state,subclonal_mult_state)
    return -likelihood

def get_point_on_simplex(n_dim):
    samples = np.random.uniform(0, 1, size=(n_dim - 1))
    samples = np.sort(samples)
    samples = np.insert(samples, 0, 0)
    samples = np.insert(samples, n_dim, 1)
    simplex_point = np.diff(samples)
    return simplex_point
#to do build a likelihood function where we sum over subclonal states and multiplicity likelihoods
def solve_timing(mult_probabilities,n_subclones,timing_matrix,full_constraints,full_constraints_sum):
    
    start_state = get_point_on_simplex(timing_matrix.shape[0])
    bounds = [(0, 1)]*start_state.size
    
    if n_subclones> 0:
        #these don't form part of the constraints
        full_constraints = np.append(full_constraints,np.zeros((full_constraints.shape[0],n_subclones)),1)
        subclone_state = np.random.uniform(0.01,0.4,size=n_subclones)
        start_state = np.concatenate([start_state,subclone_state])
        bounds.extend([(0,None)]*subclone_state.size)

    lin_cons=LinearConstraint(full_constraints,full_constraints_sum,full_constraints_sum)

    minimize_attempts = 0

    sol=minimize(mult_likelihood,start_state,constraints=lin_cons, bounds=bounds,args=(n_subclones,mult_probabilities,timing_matrix))
    
    log_likelihood = -sol.fun
    solution = sol.x
    timing_state_solution = solution[:(solution.size-n_subclones)]
    
    final_clonal_mult_state = np.matmul(timing_state_solution,timing_matrix)
    final_clonal_mult_state = final_clonal_mult_state[:mult_probabilities.major_cn]
    final_clonal_mult_state = final_clonal_mult_state/np.sum(final_clonal_mult_state)
    final_subclonal_mult_state = solution[(solution.size-n_subclones):]

    final_mult_state = np.concatenate((final_clonal_mult_state,final_subclonal_mult_state))
    final_mult_state = final_mult_state/np.sum(final_mult_state)
    return final_mult_state,timing_state_solution,log_likelihood,sol.success