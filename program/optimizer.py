from grid_solver import *
from scipy.optimize import Bounds
from scipy.optimize import differential_evolution, direct, dual_annealing, shgo


def txs_cost_function(pos, grid):
    for i in range(n):
        if 0.0 <= pos[2*i] <= 5.0 and 2.90 <= pos[2*i+1] <= 3.10:
            # to prevent the algorithm from placing the emitter inside this wall
            return 0.0
    grid.fill_power(pos.astype(np.float32))
    return -grid.get_rms()


def get_bounds(n):
    # filling bounds within where the algorithm can place emitters
    bottom_left = [0.0 for _ in range(2*n)]
    upper_right = []
    for _ in range(n):
        upper_right.append(15.0)
        upper_right.append(8.0)
    return Bounds(bottom_left, upper_right)


@measure_execution_time
def txs_optimize(grid):
    bounds = get_bounds(grid.n)
    result = differential_evolution(txs_cost_function, bounds, args=[grid], strategy='best1bin', popsize=40)
    return result
