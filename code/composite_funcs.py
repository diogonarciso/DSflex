from aux_funcs import feas_space_load_and_initial_proc
from single_constr_funcs import single_constr_proc
from n_constr_funcs import n_constr_proc

#%% Function delivers the calculation of all min, max metrics
def full_prob_calc(paths, alg_par, graph_display_settings_2D, unique_colours):

    # Load feasible space, and initial steps of processing
    feas_space, dim, bnds_feas_space, empty_feas_space_bool, feas_space_constr_redund_bool, full_dim_feas_space_bool, feas_space_corners, all_corners_square_sums, alg_par, graph_display_settings_2D = feas_space_load_and_initial_proc(paths, alg_par, graph_display_settings_2D, unique_colours)

    # Obtain the minimum distance functions per bound of the feasible space and critical regions (including informaton
    # on redundant constraints). For 2D problems, the maps of critical regions and the average distance is also displayed.
    min_optimizer_funcs, min_dist_funcs, CRs_min_min_bounds, CRs_max_min_bounds, avg_min_dist_func = single_constr_proc(feas_space, dim, feas_space_corners, feas_space_constr_redund_bool, bnds_feas_space, alg_par, graph_display_settings_2D, paths)

    # Obtain the maximum distance functions per corner of the feasible space and critical regions (including informaton
    # on redundant constraints). For 2D problems, the maps of critical regions and the average distance is also displayed.
    max_optimizers, max_dist_funcs, CRs_max_max_bounds, CRs_min_max_bounds, avg_max_dist_func = n_constr_proc(feas_space, dim, bnds_feas_space, feas_space_corners, all_corners_square_sums, alg_par, graph_display_settings_2D, paths)

    return
