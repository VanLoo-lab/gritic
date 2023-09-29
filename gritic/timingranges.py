import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
def get_equalities(mult_state,timing_matrix,wgd_constraints_matrix,wgd_timing,total_cn):
    normalising_constant = total_cn/np.sum(np.multiply(np.arange(mult_state.size)+1,mult_state))
    short_mult_state = mult_state[1:]
    
    equality_matrix = timing_matrix.T
    equality_sums = normalising_constant*short_mult_state
    if wgd_constraints_matrix is not None:
        equality_matrix = np.vstack([equality_matrix,wgd_constraints_matrix])
        equality_sums = np.concatenate([equality_sums,np.ones(wgd_constraints_matrix.shape[0])*wgd_timing])

    return equality_matrix,equality_sums

def get_range(mult_state,n_nodes,timing_matrix,sum_constraints_matrix, wgd_constraints_matrix,wgd_timing,total_cn,tol=1e-7):
    equality_matrix,equality_sums = get_equalities(mult_state,timing_matrix,wgd_constraints_matrix,wgd_timing,total_cn)
    
    timing_range = []
    for i in range(n_nodes):
        timing_range.append({})
        for optim_type in ['Min','Max']:
            c = np.zeros(n_nodes)
            c[i] = 1
            if optim_type=='Max':
                c[i]=-1
            bounds = [(-tol,1+tol)]*n_nodes
            #tol options
            #options={'dual_feasibility_tolerance':1e-5,'primal_feasibility_tolerance':1e-5,'ipm_optimality_tolerance':1e-5}
            res = linprog(c,A_eq=equality_matrix,b_eq=equality_sums,A_ub=sum_constraints_matrix,b_ub=np.ones(sum_constraints_matrix.shape[0])+tol,bounds=bounds,method='highs')
            #res = linprog(c,A_eq=equality_matrix,b_eq=equality_sums,A_ub=None,b_ub=None,bounds=bounds,method='highs',options={'dual_feasibility_tolerance':1e-5,'primal_feasibility_tolerance':1e-5,'ipm_optimality_tolerance':1e-5})

            if not res.success:
                print(equality_matrix)
                print(equality_sums)
                print(sum_constraints_matrix)
                print(res)
                print(n_nodes)
                '''_, inds = sympy.Matrix(equality_matrix).T.rref()
                print(inds)'''
                raise ValueError('No success')
            
            timing_range[i][optim_type] = res.x
    return timing_range