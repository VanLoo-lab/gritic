import numpy as np
from scipy.optimize import LinearConstraint, minimize


def unconstrained_mult_optimisation(mult_probabilities, n_subclones):
    n_dim = mult_probabilities.major_cn + n_subclones
    start_point = get_point_on_simplex(n_dim)
    bounds = [(0, 1)] * n_dim
    constraints = LinearConstraint(np.ones((1, n_dim)), 1, 1)
    sol = minimize(
        unconstrained_mult_likelihood,
        start_point,
        bounds=bounds,
        args=(mult_probabilities,),
        constraints=constraints,
    )
    if not sol.success:
        return None
    return sol.x


def unconstrained_mult_likelihood(x, mult_probabilities):
    return -mult_probabilities.evaluate_likelihood(x)


def get_point_on_simplex(n_dim):
    return np.random.dirichlet(np.ones(n_dim))
