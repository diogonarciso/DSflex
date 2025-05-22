#%% Cell 1
# Paths definition and parameters for algorithm execution

root_dir_str            = ''                           # Root directory
code_dir_str            = root_dir_str + '/code'       # files with code
prob_dir_str            = root_dir_str + '/examples'   # files with definition of feasible spaces
resl_dir_str            = root_dir_str + '/results'    # files with results from algorithms execution

problem_str             = ''

import_path             = prob_dir_str + '/' + problem_str + '.xlsx'
export_path_min         = resl_dir_str + '/' + problem_str + '_min.xlsx'
export_path_max         = resl_dir_str + '/' + problem_str + '_max.xlsx'

paths                   = {'import':import_path, 'export_min':export_path_min, 'export_max':export_path_max}

# Insert path directories to locate necessary files
import sys
sys.path.insert(1, code_dir_str)

# Import parameters
import par
alg_par                     = par.alg_par
graph_display_settings_2D   = par.graph_display_settings_2D
unique_colours              = par.unique_colours

#%%############################################################################
# Cell 2
# Solve in a single call
###############################################################################
from composite_funcs import full_prob_calc
full_prob_calc(paths, alg_par, graph_display_settings_2D, unique_colours)

#%%############################################################################
# Cell 3
# Or solve step-by-step and check results
###############################################################################
# Import libraries 
from aux_funcs import feas_space_load_and_initial_proc
from single_constr_funcs import single_constr_proc
from n_constr_funcs import n_constr_proc

#%% Cell 4
# Load feasible space, and initial steps of processing
feas_space, dim, bnds_feas_space, empty_feas_space_bool, feas_space_constr_redund_bool, full_dim_feas_space_bool, feas_space_corners, all_corners_square_sums, alg_par, graph_display_settings_2D = feas_space_load_and_initial_proc(paths, alg_par, graph_display_settings_2D, unique_colours)

#%% Cell 5
# Check results
print('Feasible space:')
print(feas_space)
#%%
print('\nDimensions:')
print(dim)
#%%
print('\nBounds for LP:')
print(bnds_feas_space)
#%%
print('\nEmpty feasible space bool:')
print(empty_feas_space_bool)
#%%
print('\nRedundant bounds of the feasible space:')
print(feas_space_constr_redund_bool)
#%%
print('\nFull dimensional feasible space:')
print(full_dim_feas_space_bool)
#%%
print('\nCorners of the feasible space (theta):')
print(feas_space_corners)
#%%
print('\nNorms of all corners coordinates:')
print(all_corners_square_sums ** (1/2))
#%%
print('\nParameters for algorithm execution:')
print(alg_par)
#%%
print('\nParameters for graphical display:')
print(graph_display_settings_2D)

#%% Cell 6
# Obtain the minimum distance functions per bound of the feasible space and critical regions (including informaton
# on redundant constraints). For 2D problems, the maps of critical regions and the average distance is also displayed.
min_optimizer_funcs, min_dist_funcs, CRs_min_min_bounds, CRs_max_min_bounds, avg_min_dist_func = single_constr_proc(feas_space, dim, feas_space_corners, feas_space_constr_redund_bool, bnds_feas_space, alg_par, graph_display_settings_2D, paths)

#%% Cell 7
# Check results
print('Optimizer functions ([theta | 1]):')
print(min_optimizer_funcs)
#%%
print('\nMinimum distance functions ([theta | 1]):')
print(min_dist_funcs)
#%%
print('\nAverage minimum distance function ([theta | 1]):')
print(avg_min_dist_func)
#%%
print('\nCR bounds: min-min ([theta | -1] <= 0))')
print(CRs_min_min_bounds)
#%%
print('\nTotal number of CR bounds: min-min')
print(CRs_min_min_bounds.shape[0]*CRs_min_min_bounds.shape[1])
#%%
print('\nCR bounds: max-min ([theta | -1] <= 0))')
print(CRs_max_min_bounds)
#%%
print('\nTotal number of CR bounds: max-min')
if (CRs_max_min_bounds != None):
    print(CRs_max_min_bounds.shape[0]*CRs_max_min_bounds.shape[1])
else:
    print(None)
#%% Cell 8
# Obtain the maximum distance functions per corner of the feasible space and critical regions (including informaton
# on redundant constraints). For 2D problems, the maps of critical regions and the average distance is also displayed.
max_optimizers, max_dist_funcs, CRs_max_max_bounds, CRs_min_max_bounds, avg_max_dist_func = n_constr_proc(feas_space, dim, bnds_feas_space, feas_space_corners, all_corners_square_sums, alg_par, graph_display_settings_2D, paths)

#%% Cell 9
# Check results
print('Optimizers/corners (theta):')
print(max_optimizers)
#%%
print('\nMaximum distance functions (excluding second order sensitivities; [theta | 1]):')
print(max_dist_funcs)
#%%
print('\nAverage minimum distance function (excluding second order sensitivities; [theta | 1]):')
print(avg_max_dist_func)
#%%
print('\nCR bounds: max-max ([theta | -1] <= 0))')
print(CRs_max_max_bounds)
#%%
print('\nTotal number of CR bounds: max-max')
print(CRs_max_max_bounds.shape[0]*CRs_max_max_bounds.shape[1])
#%%
print('\nCR bounds: min-max ([theta | -1] <= 0))')
print(CRs_min_max_bounds)
#%%
print('\nTotal number of CR bounds: min-max')
if (CRs_min_max_bounds != None):
    print(CRs_min_max_bounds.shape[0]*CRs_min_max_bounds.shape[1])
else:
    print(None)
