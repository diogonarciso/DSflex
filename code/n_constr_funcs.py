# Import libraries
from numpy import append, arctan2, argsort, array, average, dot, linspace, meshgrid, min, max, sum, reshape, zeros
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from aux_funcs import convex_hull_corners_calc, LP_ineq, LP_ineq_eq, non_redund_constr_calc, polygon_mask_2D
from report_funcs import df_calc, save_xls

#%%############################################################################
# Composite function for n constraints  (calculation of all maps)
###############################################################################
def n_constr_proc(feas_space, dim, bnds_feas_space, feas_space_corners, corners_square_sums, alg_par, graph_display_settings_2D, paths):

    # Calculate information for all distances to the t corners. In this case, optimizers are the corners themselves
    # (no additional calculation required - since this calculation was performed as an initial step).
    # Concerning the maximum distance functions to each of the t corners, this is also obtained trivially from
    # the corners coordinates matrix (C) and their norms: (d^max)^2 = theta'*theta * [1, 1, ..., 1] -2C*theta + ||C||
    # All information was already calculated in a preliminary step, and is now returned in the form [-2C | ||C|| ]
    max_dist_funcs                 = zeros((dim['t'], dim['n']+1))
    max_dist_funcs[:,0:dim['n']]   = -2 * feas_space_corners
    max_dist_funcs[:,dim['n']]     = corners_square_sums

    ###########################################################################
    # Maximum-maximum distance
    ###########################################################################

    # Calculation of all bounds for all critical regions
    CRs_max_max_bounds = all_CR_bounds_calc(feas_space, dim, feas_space_corners, max_dist_funcs, 'max')

    ###########################################################################
    # Minimum-maximum distance
    ###########################################################################

    # Calculation of all bounds for all critical regions
    if (alg_par['min_max_calc']):
        CRs_min_max_bounds = all_CR_bounds_calc(feas_space, dim, feas_space_corners, max_dist_funcs, 'min')
    else:
        CRs_min_max_bounds = None

    ###########################################################################
    # Average-maximum distance
    ###########################################################################

    # Calculation of average distance sensitivities (orders 1 and 0)
    avg_max_dist_func = avg_dist_sensitivities_calc(dim, max_dist_funcs)

    ###########################################################################
    # 2D plots
    ###########################################################################
    if ((dim['n'] == 2) and alg_par['display_2D_bool']):

        # Checking if critical regions define empty volumes, and determine redundant bounds in their definition.
        CRs_max_max_empty, CRs_max_max_redund = all_CR_redund_constr_calc(dim, CRs_max_max_bounds, bnds_feas_space, alg_par)

        # Critical regions for max-max distance problem.
        max_dist_CR_plot_2D(dim, feas_space_corners, CRs_max_max_bounds, CRs_max_max_empty, CRs_max_max_redund, bnds_feas_space, alg_par, graph_display_settings_2D, 'CR: max-max')

        if (alg_par['min_max_calc']):
            # Checking if critical regions define empty volumes, and determine redundant bounds in their definition.
            CRs_min_max_empty, CRs_min_max_redund = all_CR_redund_constr_calc(dim, CRs_min_max_bounds, bnds_feas_space, alg_par)
    
            # Critical regions for min-max distance problem.
            max_dist_CR_plot_2D(dim, feas_space_corners, CRs_min_max_bounds, CRs_min_max_empty, CRs_min_max_redund, bnds_feas_space, alg_par, graph_display_settings_2D, 'CR: min-max')

        # Plot maximum-maximum distance contours for all points in the feasible space to nearest constraints.
        n_constr_contour_plot_2D(feas_space, dim, max_dist_funcs, None, CRs_max_max_bounds, alg_par, graph_display_settings_2D, 'max_dist', 'Flexibility map: max-max')

        # Plot average-maximum distance contours for all points in the feasible space to the m constraints defining it.
        n_constr_contour_plot_2D(feas_space, dim, None, avg_max_dist_func, None, alg_par, graph_display_settings_2D, 'avg_dist', 'Flexibility map: avg-max')

        if (alg_par['min_max_calc']):
            # Plot minimum-maximum distance contours for all points in the feasible space to nearest constraints.
            n_constr_contour_plot_2D(feas_space, dim, max_dist_funcs, None, CRs_min_max_bounds, alg_par, graph_display_settings_2D, 'min_dist', 'Flexibility map: min-max')

    ###########################################################################
    # Save results to Excel
    ###########################################################################
    if (alg_par['save_results']):

        # Calculation of DataFrames from (matrix) results
        max_optimizers_df, max_dist_funcs_df, CRs_max_max_bounds_df, CRs_min_max_bounds_df, avg_max_dist_func_df = n_constr_df_calc(dim, alg_par, feas_space_corners, max_dist_funcs, CRs_max_max_bounds, CRs_min_max_bounds, avg_max_dist_func)

        # Saving results to Excel
        sheet_list = ['max_optimizers', 'max_dist_funcs', 'CRs_max_max_bounds', 'CRs_min_max_bounds', 'avg_max_dist_func']
        save_xls([max_optimizers_df, max_dist_funcs_df, CRs_max_max_bounds_df, CRs_min_max_bounds_df, avg_max_dist_func_df], sheet_list, paths['export_max'])

    return feas_space_corners, max_dist_funcs, CRs_max_max_bounds, CRs_min_max_bounds, avg_max_dist_func

#%%############################################################################
# Calculation of critical regions for max-max and min-max problems.
###############################################################################
#%% Calculation of all bounds for all critical regions.
def all_CR_bounds_calc(feas_space, dim, feas_space_corners, max_dist_funcs, mode):

    # Initialisation
    CR_bounds = zeros((dim['t'], dim['t']-1+dim['m'], dim['n']+1))

    # Calculation of bounds for all critical regions
    for i in range(dim['t']):
        CR_bounds[i,:,:] = single_CR_bounds_calc(feas_space, dim, feas_space_corners, max_dist_funcs, i, mode)

    return CR_bounds

#%% Calculation of the bounds for a critical region associated with a corner. This includes a set of optimality
# conditions from the distance comparison against other corners, and also the feasibility constraints which bound
# the corner of interest.
def single_CR_bounds_calc(feas_space, dim, feas_space_corners, max_dist_funcs, input_index, mode):

    # Initialisation 
    CR_bounds   = zeros((dim['t']-1+dim['m'], dim['n']+1))
    counter     = 0

    # First set of constraints: optimality conditions - from distance comparison against other corners of the feasible space
    for j in range(dim['t']):

        # Comparison of distance against other vertices only
        if (input_index != j):
            CR_bounds[counter] = single_opt_cond(max_dist_funcs, mode, input_index, j)
            counter += 1

    # Second set of constraints: all constraints of the feasible space
    CR_bounds[dim['t']-1:dim['t']-1+dim['m'], 0:dim['n']]   = feas_space['A'][:,:]
    CR_bounds[dim['t']-1:dim['t']-1+dim['m'], dim['n']]     = -feas_space['b'][:]

    return CR_bounds

#%% Given a corner at the i-th position and a corner at the j-th position of the complete set of corners included in the
# feasible space, this function returns the optimality condition for the definition of the critical region associated with
# the i-th corner (by comparing against the distance function for the j-th corner). This optimality is valid for the i-th
# corner (multiple by -1 to obtain the result for the j-th corner, against the i-th corner).
def single_opt_cond(max_dist_funcs, mode, i, j):

    if (mode == 'max'):
        opt_cond = -(max_dist_funcs[i,:] - max_dist_funcs[j,:])
    elif (mode == 'min'):
        opt_cond = max_dist_funcs[i,:] - max_dist_funcs[j,:]

    return opt_cond

#%%############################################################################
# Average distance of all maximum distances to the bounds of the feasible space:
# calculation of average sensitivities.
###############################################################################
def avg_dist_sensitivities_calc(dim, max_dist_funcs):

    # Coefficients for average distance calculation
    avg_max_dist_func = sum(max_dist_funcs[:,:], axis=0) / dim['t']

    return avg_max_dist_func

#%%############################################################################
# CR checks: empty space and redundant constraints.
###############################################################################
#%% Execute "single_CR_redund_constr_calc" function for all corners of the feasible space.
def all_CR_redund_constr_calc(dim, CR_bounds, bnds_feas_space, alg_par):

    # Initialialisation of matrix for empty CR check and redundancy classification
    empty_CR_bool   = zeros(dim['t']).astype(bool)
    redund_bool     = zeros((dim['t'], dim['t']-1+dim['m'])).astype(bool)

    # Calculation of empty CR checks and redundancy for all corners of the feasible space
    for i in range(dim['t']):
        empty_CR_bool[i], redund_bool[i,:] = single_CR_redund_constr_calc(dim, CR_bounds, bnds_feas_space, alg_par, i)

    return empty_CR_bool, redund_bool

#%% Checks if given set of bounds defines a non-empty critical region, and if so, then checks which of them are redundant
# in the definition of the critical region.
def single_CR_redund_constr_calc(dim, CR_bounds, bnds_feas_space, alg_par, input_index):

    # Initialisation of vector for redundancy classification
    redund_bool = zeros(dim['t']-1+dim['m']).astype(bool)

    # Check if the critical region defined via the set of given constraints includes at least one point
    status = LP_ineq(alg_par['cost'], CR_bounds[input_index,:,0:dim['n']], -CR_bounds[input_index,:,dim['n']], bnds_feas_space)
    if (status != 0):
        empty_CR_bool = True
    else:
        empty_CR_bool = False

    if not(empty_CR_bool):

        # Tests individually all constraints
        for i in range(dim['t']-1+dim['m']):

            status, x_opt = LP_ineq_eq(alg_par['cost'], CR_bounds[input_index,:,0:dim['n']], -CR_bounds[input_index,:,dim['n']], CR_bounds[input_index,i,0:dim['n']].reshape(1,dim['n']), -array(CR_bounds[input_index,i,dim['n']]), bnds_feas_space)
            if (status != 0):
                redund_bool[i] = 1

    return empty_CR_bool, redund_bool

#%%############################################################################
# Plotting critical regions: max-max and min-max cases.
###############################################################################
#%% Plot all critical regions, taking the bounds as calculated in previous steps. To this end, the corners of each of
# the critical regions must be calculated, from which a set of coloured polygons are then generated.
def max_dist_CR_plot_2D(dim, feas_space_corners, CR_bounds, empty_CR_bool, redund_bool, bnds_feas_space, alg_par, graph_display_settings_2D, title_str):

    # Initialisation on information of critical regions defining full-dim volumes
    full_dim_CR_bool = zeros(dim['t']).astype(bool)

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
    for i in range(dim['t']):

        # Print critical region (if it includes at least at least one point in the feasible space, and is in fact full-dim)
        if not(empty_CR_bool[i]):

            # Indexes calculation
            non_redund_constr_indices = non_redund_constr_calc(redund_bool[i,:])

            # Obtain the complete set of corners included in the critical region
            full_dim_CR_bool[i], unique_corners = convex_hull_corners_calc(dim, CR_bounds[i], non_redund_constr_indices, bnds_feas_space, alg_par, 'CR', True)

            # Create a polygon for each critical region, using the corresponding vertices
            if (full_dim_CR_bool[i]):
                plt.fill(unique_corners[:,0], unique_corners[:,1], graph_display_settings_2D['colours'][i])

    # Optional graphical representation of the corresponding corners
    if (graph_display_settings_2D['display_corners']):
        
        for i in range(dim['t']):
    
            # Print critical region (if it includes at least at least one point in the feasible space, and is in fact full-dim)
            if not(empty_CR_bool[i]):
    
                # Indexes calculation
                non_redund_constr_indices = non_redund_constr_calc(redund_bool[i,:])
    
                # Obtain the complete set of corners included in the critical region
                full_dim_CR_bool[i], unique_corners = convex_hull_corners_calc(dim, CR_bounds[i], non_redund_constr_indices, bnds_feas_space, alg_par, 'CR', True)
    
                # Create a polygon for each critical region, using the corresponding vertices
                if (full_dim_CR_bool[i]):
                    plt.scatter(x=feas_space_corners[i,0], y=feas_space_corners[i,1], s=graph_display_settings_2D['scatter_size'], c=graph_display_settings_2D['colours'][i])

    # Plot image
    plt.show()

    return

#%%############################################################################
# This function delivers a heatmap, based on the average distance functions, or some distance calculated 
# from a set of critical regions and the corresponding distance functions (for a 2D grid).
###############################################################################
def n_constr_contour_plot_2D(feas_space, dim, max_dist_funcs, avg_max_dist_func, bounds, alg_par, graph_display_settings_2D, mode, title_str):

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
                    dist_temp = corners_dist_calc(dim, avg_max_dist_func, array([X1[x,y], X2[x,y]]))
                    dist_array[x,y] = dist_temp
                    dist_non_zero = append(dist_non_zero, dist_temp)
    else: # if (mode == 'max_dist' or mode == 'min_dist'):
        for x in range(graph_display_settings_2D['N']):
            for y in range(graph_display_settings_2D['N']):
                for i in range(bounds.shape[0]):
                    prod_min = max(dot(bounds[i,:,:], array([X1[x,y], X2[x,y], 1])))
                    if (prod_min <= 0):
                        dist_temp = corners_dist_calc(dim, max_dist_funcs[i,:], array([X1[x,y], X2[x,y]]))
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

#%% Average maximum distance (to all corners of the feasible space).
def corners_dist_calc(dim, max_dist_func, x):

    dist = (dot(x, x) + (dot(max_dist_func[0:dim['n']], x) + max_dist_func[dim['n']])) ** (1/2)

    return dist

#%% Conversion of matrix result outputs in an equivalent DataFrame format to be saved in an Excel file.
def n_constr_df_calc(dim, alg_par, feas_space_corners, max_dist_funcs, CRs_max_max_bounds, CRs_min_max_bounds, avg_max_dist_func):

    # Optimizer functions
    max_optimizers_df       = df_calc(feas_space_corners)

    # Minimum distance functions
    max_dist_funcs_df       = df_calc(max_dist_funcs)

    # Information on critical cregions for max-max distance
    CRs_max_max_bounds_df   = df_calc(reshape(CRs_max_max_bounds, ((dim['m']+dim['t']-1)*dim['t'], dim['n']+1)))

    # Information on critical cregions for min-max distance
    if (alg_par['max_min_calc']):
        CRs_min_max_bounds_df   = df_calc(reshape(CRs_min_max_bounds, ((dim['m']+dim['t']-1)*dim['t'], dim['n']+1)))
    else:
        CRs_min_max_bounds_df   = df_calc(CRs_min_max_bounds)

    # Average minimum distance function
    avg_max_dist_func_df    = df_calc(avg_max_dist_func)

    return max_optimizers_df, max_dist_funcs_df, CRs_max_max_bounds_df, CRs_min_max_bounds_df, avg_max_dist_func_df

#%%############################################################################
# Maximum distance calculation to the corners of the feasible space
###############################################################################
#%% Calculation of distance between a given x and a given corner of the feasible space. This expression appears to be
# the simplest most general way of representing these functions, where the dependence on x (x^2 + x^1 + x^0) is employed
# and where the corresponding constants are calculated. In the case of n constraints, the structure of functions is
# modified with respect to a single constraint: the constant elements for its representation are (implicitly) the identity
# matrix for x^2, the vector of corner cordinates for x^1, and the sum of squares for x^0 - all previously calculated.
# General distance representations for all corners in this context is not particularly relevant.
def dist_calc(feas_space_corners, all_corners_square_sums, x, input_index):

    dist = (dot(x, x) + dot(-2*feas_space_corners[input_index], x) + all_corners_square_sums[input_index]) ** (1/2)

    return dist
