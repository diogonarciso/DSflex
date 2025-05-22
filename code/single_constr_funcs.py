# Import libraries
from numpy import append, arctan2, argsort, array, average, dot, identity, linspace, meshgrid, min, max, ones, reshape, zeros
from numpy.linalg import inv
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from aux_funcs import convex_hull_corners_calc, LP_ineq_eq, polygon_mask_2D, non_redund_constr_calc
from report_funcs import df_calc, save_xls

#%%############################################################################
# Composite function for single constraints (calculation of all maps).
###############################################################################
def single_constr_proc(feas_space, dim, feas_space_corners, feas_space_constr_redund_bool, bnds_feas_space, alg_par, graph_display_settings_2D, paths):

    # Calculate projection infomation for all m inequality constraints defining the feasible space. Depending on the
    # calculation mode, this matrix may include information on distance functions only, or also optimizers.
    opt_min = all_proj(feas_space, dim, alg_par)

    if (alg_par['min_calc_mode'] == 'dist_only'):
        min_optimizer_funcs = None
        min_dist_funcs      = opt_min
    else:
        min_optimizer_funcs = opt_min[:,0:dim['n'],:]
        min_dist_funcs      = opt_min[:,dim['n'],:]

    ###########################################################################
    # Minimum-minimum distance
    ###########################################################################

    # Calculation of all bounds for all critical regions
    CRs_min_min_bounds = all_CR_min_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_bool, alg_par)

    ###########################################################################
    # Maximum-minimum distance
    ###########################################################################

    # Calculation of all bounds for all critical regions
    if (alg_par['max_min_calc']):
        CRs_max_min_bounds = all_CR_max_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_bool, alg_par)
    else:
        CRs_max_min_bounds = None

    ###########################################################################
    # Average-minimum distance
    ###########################################################################

    # Calculation of average sensitivities
    avg_min_dist_func = avg_dist_sensitivities_calc(min_dist_funcs)

    ###########################################################################
    # 2D plots
    ###########################################################################
    if ((dim['n'] == 2) and alg_par['display_2D_bool']):

        # Calculation of redundant bounds on the definition of all critical regions (min-min)
        CRs_min_min_redund = all_CR_min_redund_constr_calc(dim, CRs_min_min_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par)

        # Plot critical regions (minimum-minimum distance to the bounds of the feasible space problem).
        min_dist_CR_plot_2D(dim, feas_space_constr_redund_bool, CRs_min_min_bounds, CRs_min_min_redund, bnds_feas_space, alg_par, graph_display_settings_2D, 'CR: min-min')

        if (alg_par['max_min_calc']):
            # Calculation of redundant bounds on the definition of all critical regions (max-min)
            CRs_max_min_redund = all_CR_max_redund_constr_calc(dim, CRs_max_min_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par)
    
            # Plot critical regions (maximum-minimum distance to the bounds of the feasible space problem).
            min_dist_CR_plot_2D(dim, feas_space_constr_redund_bool, CRs_max_min_bounds, CRs_max_min_redund, bnds_feas_space, alg_par, graph_display_settings_2D, 'CR: max-min')

        # Plot minimum-minimum distance contours for all points in the feasible space to nearest constraints.
        single_constr_contour_plot_2D(feas_space, dim, min_dist_funcs, None, CRs_min_min_bounds, alg_par, graph_display_settings_2D, 'min_dist', 'Flexibility map: min-min')

        # Plot average-minimum distance contours for all points in the feasible space to the m constraints defining it.
        single_constr_contour_plot_2D(feas_space, dim, None, avg_min_dist_func, None, alg_par, graph_display_settings_2D, 'avg_dist', 'Flexibility map: avg-min')

        if (alg_par['max_min_calc']):
            # Plot maximum-minimum distance contour for all points in the feasible space to nearest constraints.
            single_constr_contour_plot_2D(feas_space, dim, min_dist_funcs, None, CRs_max_min_bounds, alg_par, graph_display_settings_2D, 'max_dist', 'Flexibility map: max-min')

    ###########################################################################
    # Save results to Excel
    ###########################################################################
    if (alg_par['save_results']):

        # Calculation of DataFrames from (matrix) results
        min_optimizer_funcs_df, min_dist_funcs_df, CRs_min_min_bounds_df, CRs_max_min_bounds_df, avg_min_dist_func_df = single_constr_df_calc(dim, alg_par, min_optimizer_funcs, min_dist_funcs, CRs_min_min_bounds, CRs_max_min_bounds, avg_min_dist_func)

        # Saving results to Excel
        sheet_list = ['min_optimizer_funcs', 'min_dist_funcs', 'CRs_min_min_bounds', 'CRs_max_min_bounds', 'avg_min_dist_func']
        save_xls([min_optimizer_funcs_df, min_dist_funcs_df, CRs_min_min_bounds_df, CRs_max_min_bounds_df, avg_min_dist_func_df], sheet_list, paths['export_min'])

    return min_optimizer_funcs, min_dist_funcs, CRs_min_min_bounds, CRs_max_min_bounds, avg_min_dist_func

#%%############################################################################
# Projection / minimum distance functions to all bounds of the feasible space
###############################################################################
#%% Function "single_constr_proj" is applied to all m single constraints for the calculation of minimum distances/projections; results are delivered
# in an output matrix of suitable dimensons.
def all_proj(feas_space, dim, alg_par):

    # Outputs initialisation
    if (alg_par['min_calc_mode'] == 'dist_only'):

        # In this case, the output matrix includes m rows for each of the corresponding sensitivities for the distance functions
        output  = zeros((dim['m'], dim['n']+1))

    elif (alg_par['min_calc_mode'] == 'projection'):

        # In this case, the output matrix includes m matrices for each of the corresponding sensitivities for the optimizer and distance functions
        output  = zeros((dim['m'], dim['n']+1, dim['n']+1))

    # The results for each of the m single constraints is calculated individually and the output matrix filled accordingly
    for i in range(dim['m']):
        output[i] = single_proj(feas_space, dim, alg_par, i)

    return output

#%% Given an affine bound in the n-dim space, this function returns a vector for the calculation of the
# minimum distance between a point theta and the given bound. Two modes are defined: in the first mode
# only the output vector for minimum distance calculation is returned; in the second mode, the projection
# from theta into the given bound is also calculated (the optimizer, as a function of theta).
def single_proj(feas_space, dim, alg_par, input_index):

    # Calculate norm of normal vector (pointing towards the interior region of the feasible space)
    norm = (sum(feas_space['A'][input_index,:] ** 2)) ** (1/2)

    if (alg_par['min_calc_mode'] == 'dist_only'):

        # Initialisation - in this case, includes only the minimum distance information (sensitivities to theta
        # plus the constant independent term in the last position)
        output              = zeros(dim['n']+1)

        # Distance vector calculated trivially from constraint and norm
        output[0:dim['n']]  = -feas_space['A'][input_index,:] / norm
        output[dim['n']]    = feas_space['b'][input_index] / norm

    elif (alg_par['min_calc_mode'] == 'projection'):

        # Initialisation - projection problem statement via a linear system
        proj_mat = identity(dim['n']+1)

        # Fill last row and last column of projection matrix
        proj_mat[dim['n'], 0:dim['n']]  = -feas_space['A'][input_index,:]
        proj_mat[0:dim['n'], dim['n']]  = -feas_space['A'][input_index,:]
        proj_mat[dim['n'], dim['n']]    = 0

        # Calculate inverse matrix - output is in this case sensitivities for the optimizer and distance functions
        output  = inv(proj_mat)

        # In the original form, the last row of the inverse matrix delivers the sensitivities for the calculation of a multiplicative factor (l)
        # with respect to the given normal vector; multiplying this by the norm, the output matrix is updated to return the sensitivities for
        # the calculation of the minimum distance between a given theta and given bound (equivalent to first mode)
        output[dim['n'],:]  = output[dim['n'],:] * norm

        # The last row of the inverse matrix is multiplied by "-b" so that all calculations are consistent with the multiplication of this matrix
        # with vector [theta | 1]
        output[:,dim['n']]  = output[:,dim['n']] * (-feas_space['b'][input_index])

    return output

#%%############################################################################
# Minimum distance of all minimum distances to the bounds of the feasible space
# (multiparametric programming problem): calculation of critical regions.
###############################################################################
#%% This function delivers the bounds of critical regions for all single constraint (active sets).
def all_CR_min_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_bool, alg_par):

    # Initialisation of vector of CRs bounds
    CR_bounds = zeros((dim['m'], dim['m'], dim['n']+1))

    # Calculations of bounds per each of the m single constraint (active sets)
    for i in range(dim['m']):
        if not(feas_space_constr_redund_bool[i]):
            CR_bounds[i,:,:] = single_CR_min_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_bool, i)

    return CR_bounds

#%% For a given single constraint (active set), this function obtains the set of m-1 optimality constraints, which are
# obtained by comparing the corresponding minimum distance with the distances from the remaining m-1 constraints (criterium
# for optimality). A single feasibility constraint (the original feasible space constraint is also added).
def single_CR_min_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_check, input_index):

    # Initialisation: x-dependent terms in the first n positions, independent term in last position
    CR_bounds = zeros((dim['m'], dim['n']+1))

    # Define relevant bounds (m-1 optimality constriants and a single feasibility constraint)
    for j in range(dim['m']):

        # This relates to the comparison of distance between the reference distance for the i-th constraint and any
        # other constraint j (not the same as i), for which it makes sense to compare the distances and obtain the
        # corresponding optimality constraint
        if (j != input_index):
            if not(feas_space_constr_redund_check[j]):
                CR_bounds[j,:] = single_opt_cond(min_dist_funcs, 'min', input_index, j)

        # In this case, a comparison of the same distance function would not make sense. Here, instead of an
        # optimality condition, a feasibility condition is defined which is the original i-th constraint of the
        # feasible space
        else:
            CR_bounds[j,0:dim['n']] = feas_space['A'][input_index,:]
            CR_bounds[j,dim['n']]   = -feas_space['b'][input_index]

    return CR_bounds

#%% This function check redundancy of all constraints in the definition of all critical regions (single constraint).
def all_CR_min_redund_constr_calc(dim, all_CR_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par):

    # Initialialisation of matrix for redundancy classification
    redund_bool = zeros((dim['m'], dim['m'])).astype(bool)

    # For all single constraints, obtains redundancy information on the definition of their critical regions; for those
    # single constraints which have been identified as redundant, an all-False vector is defined (does not require 
    # calculation from the "single" function).
    for i in range(dim['m']):
        if (feas_space_constr_redund_bool[i]):
            redund_bool[i,:] = ones(dim['m'])
        else:
            redund_bool[i,:] = single_CR_min_redund_constr_calc(dim, all_CR_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par, i)

    return redund_bool

#%% This function check redundancy of all constraints in the definition of a given critical region (single constraint). It is 
# adapted from function "feas_space_redund_constr_calc", where some additional considerations are required to perform this
# assessment.
def single_CR_min_redund_constr_calc(dim, all_CR_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par, input_index):

    # Initialialisation of vector for redundancy classification
    redund_bool = zeros(dim['m']).astype(bool)

    # Tests individually all constraints via LP in the context of all bounds of the convex hull; at position "i = index_input"
    # is the feasibility constraint (defining the feasible space), which does not require assessment. A refined logic is
    # implemented to dismiss the comparison with redundant constraints (on the definition of the feasible space), which are
    # not relevant for critical region construction.
    for i in range(dim['m']):
        if not(i == input_index):
            if (feas_space_constr_redund_bool[i]):
                redund_bool[i] = 1
            else:
                status, x_opt = LP_ineq_eq(alg_par['cost'], all_CR_bounds[input_index,:,0:dim['n']], -all_CR_bounds[input_index,:,dim['n']], all_CR_bounds[input_index,i,0:dim['n']].reshape(1,dim['n']), -array(all_CR_bounds[input_index,i,dim['n']]), bnds_feas_space)
                if (status != 0):
                    redund_bool[i] = 1

    return redund_bool

#%%############################################################################
# Maximum distance of all minimum distances to all bounds of the feasible space
# (multiparametric programming problem): calculation of critical regions.
###############################################################################
#%% This function delivers the bounds of critical regions for all single constraint (active sets).
def all_CR_max_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_bool, alg_par):

    # Initialisation of vector of CRs bounds
    CR_bounds = zeros((dim['m'], 2*(dim['m']-1), dim['n']+1))

    # Calculations of bounds per each of the m single constraint (active sets)
    for i in range(dim['m']):
        if not(feas_space_constr_redund_bool[i]):
            CR_bounds[i,:,:] = single_CR_max_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_bool, i)

    return CR_bounds

#%% For a given single constraint (active set), this function obtains the set of m-1 optimality constraints, which are
# obtained by comparing the corresponding minimum distance with the distances from the remaining m-1 constraints (criterium
# for optimality). m+1 feasibility constraint (all constraints of the feasible space, except the original constraint).
def single_CR_max_bounds_calc(feas_space, dim, min_dist_funcs, feas_space_constr_redund_check, input_index):

    # Initialisation: x-dependent terms in the first n positions, independent term in last position
    CR_bounds   = zeros((2*(dim['m']-1), dim['n']+1))
    counter     = 0

    # Define relevant bounds (first set: m-1 optimality constraints)
    for j in range(dim['m']):

        # This relates to the comparison of distance between the reference distance for the i-th constraint and any
        # other constraint j (not the same as i), for which it makes sense to compare the distances and obtain the
        # corresponding optimality constraint
        if (j != input_index):
            if not(feas_space_constr_redund_check[j]):
                CR_bounds[counter,:] = single_opt_cond(min_dist_funcs, 'max', input_index, j)

            # Increment
            counter += 1

    # Second set: all bounds of the feasible space except at "index_input" must also be considered (m-1 feasibility constraints).
    # The first step consists of identifying which bounds of the feasible space are included: all but the "index_input"^th bound
    select_bool = ones(dim['m']).astype(bool)
    select_bool[input_index] = 0

    # Adding the corresponding bounds to the results matrix
    CR_bounds[dim['m']-1:2*(dim['m']-1),0:dim['n']] = feas_space['A'][select_bool,:]
    CR_bounds[dim['m']-1:2*(dim['m']-1),dim['n']]   = -feas_space['b'][select_bool]

    return CR_bounds

#%% This function check redundancy of all constraints in the definition of all critical regions (single constraint).
def all_CR_max_redund_constr_calc(dim, all_CR_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par):

    # Initialialisation of matrix for redundancy classification
    redund_bool = zeros((dim['m'], 2*(dim['m']-1))).astype(bool)

    # For all single constraints, obtains redundancy information on the definition of their critical regions; for those
    # single constraints which have been identified as redundant, an all-False vector is defined (does not require 
    # calculation from the "single" function).
    for i in range(dim['m']):
        if (feas_space_constr_redund_bool[i]):
            redund_bool[i,:] = ones(2*(dim['m']-1))
        else:
            redund_bool[i,:] = single_CR_max_redund_constr_calc(dim, all_CR_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par, i)

    return redund_bool

#%% This function check redundancy of all constraints in the definition of a given critical region (single constraint). It is 
# adapted from function "feas_space_redund_constr_calc", where some additional considerations are required to perform this
# assessment.
def single_CR_max_redund_constr_calc(dim, all_CR_max_bounds, feas_space_constr_redund_bool, bnds_feas_space, alg_par, input_index):

    # Initialialisation of vector for redundancy classification
    redund_bool = zeros(2*(dim['m']-1)).astype(bool)

    # Tests individually all constraints via LP in the context of all bounds of the convex hull; at position "i = index_input"
    # is the feasibility constraint (defining the feasible space), which does not require assessment. A refined logic is
    # implemented to dismiss the comparison with redundant constraints (on the definition of the feasible space), which are
    # not relevant for critical region construction.
    for i in range(dim['m']-1):
        if (feas_space_constr_redund_bool[i]):
            redund_bool[i] = 1
        else:
            status, x_opt = LP_ineq_eq(alg_par['cost'], all_CR_max_bounds[input_index,:,0:dim['n']], -all_CR_max_bounds[input_index,:,dim['n']], all_CR_max_bounds[input_index,i,0:dim['n']].reshape(1,dim['n']), -array(all_CR_max_bounds[input_index,i,dim['n']]), bnds_feas_space)
            if (status != 0):
                redund_bool[i] = 1

    for i in range(dim['m']-1):
        status, x_opt = LP_ineq_eq(alg_par['cost'], all_CR_max_bounds[input_index,:,0:dim['n']], -all_CR_max_bounds[input_index,:,dim['n']], all_CR_max_bounds[input_index,dim['m']-1+i,0:dim['n']].reshape(1,dim['n']), -array(all_CR_max_bounds[input_index,dim['m']-1+i,dim['n']]), bnds_feas_space)
        if (status != 0):
            redund_bool[dim['m']-1+i] = 1

    return redund_bool

#%%############################################################################
# Optimality conditions: minimum and maximum (of minimum distance).
###############################################################################
#%% This function delivers an optimality condition for the i-th single active constraint, by comparing the corresponding
# distance vector, with the distance vector obtained for the j-th single active constraint: distance(i) <= distance(j).
# A vector is obtained as a result of this calculation, which is specific for the critical region of the i-th constraint
# (or active set).
def single_opt_cond(min_dist_funcs, mode, i, j):

    # Subtraction of the two corresponding vectors (from dist_i <= dist_j or dist_i >= dist_j)
    if (mode == 'min'):
        opt_cond = min_dist_funcs[i,:] - min_dist_funcs[j,:]
    elif (mode == 'max'):
        opt_cond = -(min_dist_funcs[i,:] - min_dist_funcs[j,:])

    return opt_cond

#%%############################################################################
# Average distance of all minimum distances to the bounds of the feasible space:
# calculation of average sensitivities.
###############################################################################
def avg_dist_sensitivities_calc(min_dist_funcs):

    # Sensitivity of average distance with respect to all coordinates of x (and the constant, final position)
    avg_min_dist_func = average(min_dist_funcs[:,:], axis=0)

    return avg_min_dist_func

#%%############################################################################
# Plotting critical regions: min-min and max-min cases.
###############################################################################
#%% Plot all critical regions, taking the bounds as calculated in previous steps. To this end, the corners of each of
# the critical regions must be calculated, from which a set of coloured polygons are then generated.
def min_dist_CR_plot_2D(dim, feas_space_redund, CRs_bounds, CRs_redund, bnds_feas_space, alg_par, graph_display_settings_2D, title_str, print_bool=False, theta=None, x_proj=None):

    # Initialisation on information of critical regions defining full-dim volumes
    CRs_full_dim = zeros(dim['m']).astype(bool)

    # Initialise plot
    plt.figure(figsize=(graph_display_settings_2D['fig_width'], graph_display_settings_2D['fig_height']))
    plt.title(title_str, fontsize=graph_display_settings_2D['font_size_title'])
    plt.xlabel(xlabel=r'$x_1$', fontsize=graph_display_settings_2D['font_size_axis'])
    plt.ylabel(ylabel=r'$x_2$', fontsize=graph_display_settings_2D['font_size_axis'])
    plt.xticks(fontsize=graph_display_settings_2D['font_size_axis'])
    plt.yticks(fontsize=graph_display_settings_2D['font_size_axis'])

    plt.xlim([graph_display_settings_2D['min_x1'], graph_display_settings_2D['max_x1']])
    plt.ylim([graph_display_settings_2D['min_x2'], graph_display_settings_2D['max_x2']])

    # Graphical representation of critical regions
    for i in range(dim['m']):

        # Only non redundant constraints of the feasible space yield critical regions
        if not(feas_space_redund[i]):

            # Making use of the vector of non redundant bounds (critical region definition), obtain the corresponding
            # vector of index positions
            non_redund_constr_indices = non_redund_constr_calc(CRs_redund[i])

            # Obtain the complete set of corners included in the critical region
            CRs_full_dim[i], unique_corners = convex_hull_corners_calc(dim, CRs_bounds[i], non_redund_constr_indices, bnds_feas_space, alg_par, 'CR', True)

            # Create a polygon for each critical region, using the corresponding vertices
            if (CRs_full_dim[i]):
                plt.fill(unique_corners[:,0], unique_corners[:,1], graph_display_settings_2D['colours'][i])

    # Print a single point - typically a projection on a bound
    if (print_bool):
        plt.scatter(x=theta[0], y=theta[1], s=graph_display_settings_2D['scatter_size'], c='black')
        plt.scatter(x=x_proj[0], y=x_proj[1], s=graph_display_settings_2D['scatter_size'], c='black')

    # Plot image
    plt.show()

    return

#%%############################################################################
# This function delivers a heatmap, based on the average distance functions, or some distance calculated 
# from a set of critical regions and the corresponding distance functions (for a 2D grid).
###############################################################################
def single_constr_contour_plot_2D(feas_space, dim, min_dist_funcs, avg_min_dist_func, bounds, alg_par, graph_display_settings_2D, mode, title_str):

    ###########################################################################
    # 2D arrays for x1, x2 and the average distance at all of their coordinates
    ###########################################################################

    # Auxiliary arrays of x and y to generate 2D plot
    linspace_x1 = linspace(graph_display_settings_2D['min_x1'], graph_display_settings_2D['max_x1'], graph_display_settings_2D['N']) 
    linspace_x2 = linspace(graph_display_settings_2D['min_x2'], graph_display_settings_2D['max_x2'], graph_display_settings_2D['N']) 

    # 2D meshgrids
    [X1, X2] = meshgrid(linspace_x1, linspace_x2)

    # Initialisation of array of distances
    dist_array      = zeros((graph_display_settings_2D['N'], graph_display_settings_2D['N']))
    dist_non_zero   = zeros(0)

    # Calculation of distance for all points in the meshgrid. Some of these points will not actually be shown
    # in the final plot, since the meshgrid typically extends beyond the feasible space
    if (mode == 'avg_dist'):
        for x in range(graph_display_settings_2D['N']):
            for y in range(graph_display_settings_2D['N']):
                prod_min = max(dot(feas_space['A'], array([X1[x,y], X2[x,y]])) - feas_space['b'])
                if (prod_min <= 0):
                    dist_temp = dot(avg_min_dist_func, array([X1[x,y], X2[x,y], 1]))
                    dist_array[x,y] = dist_temp
                    dist_non_zero = append(dist_non_zero, dist_temp)
    else: # if (mode == 'max_dist' or mode == 'min_dist'):
        for x in range(graph_display_settings_2D['N']):
            for y in range(graph_display_settings_2D['N']):
                for i in range(bounds.shape[0]):
                    prod_min = max(dot(bounds[i,:,:], array([X1[x,y], X2[x,y], 1])))
                    if (prod_min <= 0):
                        dist_temp = dot(min_dist_funcs[i,:], array([X1[x,y], X2[x,y], 1]))
                        dist_array[x,y] = dist_temp
                        dist_non_zero = append(dist_non_zero, dist_temp)
                        break

    ###########################################################################
    # Calculation of min and max values for colour bar: these are set to the minimum and maximum distance values found
    # in all corners of the feasible space
    ###########################################################################

    min_dist = min(dist_non_zero)
    max_dist = max(dist_non_zero)

    # Updating distances with value 0 to minimum distance to improve colorbar display
    for x in range(graph_display_settings_2D['N']):
        for y in range(graph_display_settings_2D['N']):
            if(dist_array[x,y] == 0):
                dist_array[x,y] = min_dist

    ###########################################################################
    # Figure set-up
    ###########################################################################

    # Initialisation
    plt.figure(figsize=(graph_display_settings_2D['fig_width'], graph_display_settings_2D['fig_height']))

    # Data to figure
    CS = plt.contourf(X1, X2, dist_array)

    # Title and axes
    plt.title(title_str, fontsize=graph_display_settings_2D['font_size_title'])
    plt.xlabel(xlabel=r'$x_1$', fontsize=graph_display_settings_2D['font_size_axis'])
    plt.ylabel(ylabel=r'$x_2$', fontsize=graph_display_settings_2D['font_size_axis'])
    plt.xticks(fontsize=graph_display_settings_2D['font_size_axis'])
    plt.yticks(fontsize=graph_display_settings_2D['font_size_axis'])
    plt.xlim([graph_display_settings_2D['min_x1'], graph_display_settings_2D['max_x1']])
    plt.ylim([graph_display_settings_2D['min_x2'], graph_display_settings_2D['max_x2']])

    # Colorbar
    fmt = lambda x, pos: '{:.2f}'.format(x)
    colorbar_ticks = linspace(min_dist, max_dist, 10, endpoint=True)
    cbar = plt.colorbar(CS, ticks=colorbar_ticks, format=FuncFormatter(fmt))
    cbar.ax.set_ylim(min_dist, max_dist)
    cbar.ax.tick_params(labelsize=graph_display_settings_2D['label_size'])

    # Mask polygons - set of polygons marked in white to emphasize feasible space
    for i in range(dim['m']):
        mask_corners = polygon_mask_2D(feas_space, dim, alg_par, graph_display_settings_2D, i)
        order = argsort(arctan2(mask_corners[:,0] - average(mask_corners[:,0]), mask_corners[:,1] - average(mask_corners[:,1])))
        plt.fill(mask_corners[order,0], mask_corners[order,1], "w")

    # Plot
    plt.show()

    return

#%% Conversion of matrix result outputs in an equivalent DataFrame format to be saved in an Excel file.
def single_constr_df_calc(dim, alg_par, min_optimizer_funcs, min_dist_funcs, CRs_min_min_bounds, CRs_max_min_bounds, avg_min_dist_func):

    # Optimizer functions
    if (alg_par['min_calc_mode'] == 'dist_only'):
        min_optimizer_funcs_df  = df_calc(min_optimizer_funcs)
    else: # if (alg_par['min_calc_mode'] == 'projection'):
        min_optimizer_funcs_df  = df_calc(reshape(min_optimizer_funcs, ((dim['n'])*dim['m'], dim['n']+1)))

    # Minimum distance functions
    min_dist_funcs_df       = df_calc(min_dist_funcs)

    # Information on critical cregions for min-min distance
    CRs_min_min_bounds_df   = df_calc(reshape(CRs_min_min_bounds, (dim['m']**2, dim['n']+1)))

    # Information on critical cregions for max-min distance
    if (alg_par['max_min_calc']):
        CRs_max_min_bounds_df   = df_calc(reshape(CRs_max_min_bounds, (dim['m']*(dim['m']-1)*2, dim['n']+1)))
    else:
        CRs_max_min_bounds_df   = df_calc(CRs_max_min_bounds)

    # Average minimum distance function
    avg_min_dist_func_df    = df_calc(avg_min_dist_func)

    return min_optimizer_funcs_df, min_dist_funcs_df, CRs_min_min_bounds_df, CRs_max_min_bounds_df, avg_min_dist_func_df
