# Import libraries
from pandas import read_excel
from numpy import append, arctan2, argsort, array, average, dot, min, max, ones, unique, zeros
from numpy.linalg import inv, matrix_rank
from scipy.optimize import linprog
from math import ceil, comb
from itertools import combinations

#%%############################################################################
# Composite function for initital processing of feasible space
###############################################################################
def feas_space_load_and_initial_proc(paths, alg_par, graph_display_settings_2D, unique_colours):

    # Load Ax < b matrices from Excel file, and configure a set of auxiliary parameters for algorithm execution
    feas_space, dim, bnds_feas_space, alg_par = problem_load_set_up(paths['import'], alg_par)

    # Check if constraints define an empty space, and redundancy assessment for all constraints defining the feasible space.
    empty_feas_space_bool, feas_space_constr_redund_bool = feas_space_redund_constr_calc(feas_space, dim, bnds_feas_space, alg_par)

    # Calculations proceed only if feasible space is not empty
    if not(empty_feas_space_bool):

        # Calculate indices of non redundant constraints in the definition of the feasible space.
        feas_space_non_redund_constr_indices = non_redund_constr_calc(feas_space_constr_redund_bool)

        # Calculate all corners of the feasible space, and the indices of constraints active at their coordinates.
        full_dim_feas_space_bool, feas_space_corners = convex_hull_corners_calc(dim, feas_space, feas_space_non_redund_constr_indices, bnds_feas_space, alg_par, 'feas_space', True)

        # Update dim dictionary (information on the number of corners of the feasible space)
        dim, graph_display_settings_2D = dicts_update(dim, feas_space_corners, alg_par, graph_display_settings_2D, unique_colours)

        # Sum of quare sums for all corners - auxiliary for distance calculations
        all_corners_square_sums = corners_coords_square_sum(dim, feas_space_corners)

    else:

        # If space is empty, no corners exist
        feas_space_corners = None

    return feas_space, dim, bnds_feas_space, empty_feas_space_bool, feas_space_constr_redund_bool, full_dim_feas_space_bool, feas_space_corners, all_corners_square_sums, alg_par, graph_display_settings_2D

#%%############################################################################
# Feasible space functions
###############################################################################
#%% Read information on feasible space from Excel file, and create a set of dictionaries to be used as arguments for
# developed functions.
def problem_load_set_up(import_path, alg_par):

    # Read feasible space information from path
    df = read_excel(import_path, sheet_name='feas_space', index_col=0)

    # Create arrays A and b, and dictionary for feasible space
    b = df.iloc[:,0].to_numpy()
    A = df.iloc[:,1:].to_numpy()
    feas_space = {'A':A, 'b':b}

    # Calculate problem dimensions and disctionary
    m = A.shape[0]
    n = A.shape[1]
    dim = {'n':n, 'm':m}

    # Update dictionary of parameters
    alg_par['cost'] = ones(dim['n'])

    # Bounds per variable and initialisation - for optimization calculations
    bnds_single      = (-alg_par['large_pos_number'], alg_par['large_pos_number'])
    bnds_feas_space  = []

    # Create list of bounds for n dimensions
    for i in range(dim['n']):
        bnds_feas_space.append(bnds_single)

    return feas_space, dim, bnds_feas_space, alg_par

#%% Check redundancy of all bounds in the definition of a convex hull defined via a set of bounds of the form Ax <= b.
def feas_space_redund_constr_calc(feas_space, dim, bnds_feas_space, alg_par):

    # Check if feasible space is not empty (in case a feasible solution is found: status=0). In this case, LP is solved
    # only with the m constraints defining the feasible space.
    status = LP_ineq(alg_par['cost'], feas_space['A'], feas_space['b'], bnds_feas_space)

    if (status == 0):
        empty_feas_space_bool = False
    else:
        empty_feas_space_bool = True

    # Initialialisation of vector for redundancy classification
    constr_redund_bool = zeros(dim['m']).astype(bool)

    if not(empty_feas_space_bool):

        # Tests individually all constraints via LP in the context of all bounds of the convex hull
        for i in range(dim['m']):
            status, x_opt = LP_ineq_eq(alg_par['cost'], feas_space['A'], feas_space['b'], feas_space['A'][i,:].reshape(1,dim['n']), array([feas_space['b'][i]]), bnds_feas_space)
            if (status != 0):
                constr_redund_bool[i] = 1

    return empty_feas_space_bool, constr_redund_bool

#%% Update dicts, based on the corners of the feasible space.
def dicts_update(dim, feas_space_corners, alg_par, graph_display_settings_2D, unique_colours):

    # dim update
    dim['t'] = feas_space_corners.shape[0]

    # graph_display_settings_2D update - only if 2D plots are to be used
    if ((dim['n'] == 2) and alg_par['display_2D_bool']):

        # Bounds of 2D graph display
        min_x1_corners = min(feas_space_corners[:,0]) - graph_display_settings_2D['edge_gap']
        max_x1_corners = max(feas_space_corners[:,0]) + graph_display_settings_2D['edge_gap']
        min_x2_corners = min(feas_space_corners[:,1]) - graph_display_settings_2D['edge_gap']
        max_x2_corners = max(feas_space_corners[:,1]) + graph_display_settings_2D['edge_gap']

        # graph_display_settings_2D update - bounds
        graph_display_settings_2D['min_x1'] = min_x1_corners
        graph_display_settings_2D['max_x1'] = max_x1_corners
        graph_display_settings_2D['min_x2'] = min_x2_corners
        graph_display_settings_2D['max_x2'] = max_x2_corners

        # graph_display_settings_2D update - corners
        corners = array([[min_x1_corners, min_x2_corners],
                         [min_x1_corners, max_x2_corners],
                         [max_x1_corners, min_x2_corners],
                         [max_x1_corners, max_x2_corners]])

        graph_display_settings_2D['corners'] = corners

        # Colours initialisation
        colours = []

        # Repeat list of colours if solution includes more than 20 colours
        min_c_d = ceil(min(array([dim['m'], dim['t']])/len(unique_colours)))

        for i in range(min_c_d):
            colours = colours + unique_colours

        graph_display_settings_2D['colours'] = colours

    return dim, graph_display_settings_2D

#%%############################################################################
# General functions (used in multiple calculations)
###############################################################################
#%% This function converts the boolean vector containing the information on redundancy into an array listing the positions
# of all non redundant constraints.
def non_redund_constr_calc(constr_redund_bool):

    # Total number of constraints - for bool input vector
    total_constr = constr_redund_bool.shape[0]

    # Initialisation
    total_non_redund_constr     = total_constr - sum(constr_redund_bool[:])
    non_redund_constr_indices   = zeros(total_non_redund_constr).astype(int)
    counter                     = 0

    # Conversion from boolean vector to an equivalent index positions vector
    for i in range(total_constr):
        if not(constr_redund_bool[i]):
            non_redund_constr_indices[counter] = i
            counter += 1

    return non_redund_constr_indices

#%% This function delivers the full list of corners from a convex hull defined via a set of bounds of the form Ax <= b.
# To this end, it takes all combinations of n non-redundant bounds defining the convex hull (minimum requirement to define
# a corner) as a set of n equality constraints and the full set of constraints defining the convex hull (as inequality
# constraints) to check via LP which indeed define vertices. To achieve this, a vector listing the index of all non
# redundant constraints is calculated beforehand and used as a function argument. This function is applicable to the
# feasible space and to critical regions.
def convex_hull_corners_calc(dim, input_bounds, non_redundant_constr_indices, bnds_feas_space, alg_par, mode, order_bool):

    # Total number of combinations of constraints potentially defining vertices.
    total_combs = comb(non_redundant_constr_indices.shape[0], dim['n'])

    # Initialisation: 
    statuses        = zeros(total_combs)
    corners_raw     = zeros((total_combs, dim['n']))    # not all results stored in this array are necessarily a corner - a final processing step is required to exclude irrelevant combinations
    total_corners   = 0
    counter         = 0                                 # auxiliary counter for cycle below

    # Test all candidate combinations possibly defining a corner in the feasible space (via LP)
    if (mode == 'feas_space'):

        for i in combinations(non_redundant_constr_indices, dim['n']):

            # Test single combination of constraints
            statuses[counter], corners_raw[counter,:] = LP_ineq_eq(
                alg_par['cost'], input_bounds['A'][non_redundant_constr_indices,:],
                input_bounds['b'][non_redundant_constr_indices], input_bounds['A'][array(i),:],
                input_bounds['b'][array(i)], bnds_feas_space)

            # Increment on statuses and corners index
            counter += 1

    elif (mode == 'CR'):

        for i in combinations(non_redundant_constr_indices, dim['n']):

            # Test single combination of constraints
            statuses[counter], corners_raw[counter,:] = LP_ineq_eq(
                alg_par['cost'], input_bounds[non_redundant_constr_indices,0:dim['n']],
                -input_bounds[non_redundant_constr_indices,dim['n']], input_bounds[array(i),0:dim['n']],
                -input_bounds[array(i),dim['n']], bnds_feas_space)

            # Increment on statuses and corners index
            counter += 1

    # Process results from first test: keep only results from feasible LP calculations (where a corner is found)
    for i in range(total_combs):
        if (statuses[i] == 0):
            total_corners += 1

    # Initialisation of full list of corners included in the convex hull
    corners = zeros((total_corners, dim['n']))
    counter = 0

    # Final list of corners (a final pre-processing step is still required after this to remove duplicates)
    for i in range(total_combs):
        if (statuses[i] == 0):
            corners[counter]            = corners_raw[i]
            counter += 1

    # Remove duplicates (if any), check they define a full-dimensional region, and update their order
    # for correct display (in 2D only)
    full_dim_CR_bool, unique_corners = corners_pre_proc(dim, corners, alg_par, order_bool)

    return full_dim_CR_bool, unique_corners

#%% This function takes a set of corners for critical region representation, checks if they define a full-dimensional
# critical region, and if so adjusts their order to ensure the critical region is printed as a full polygon.
def corners_pre_proc(dim, corners, alg_par, order_bool):

    # Eliminate duplicates from raw set of corners (directly from all combinations of interceptions)
    unique_corners          = unique(corners, axis=0)
    total_unique_corners    = unique_corners.shape[0]

    # Check if the unique set of corners defines a full-dimensional critical region
    if (total_unique_corners > dim['n']):

        # Initialisation: polygon directions
        polygon_direct = zeros((total_unique_corners-1, dim['n']))

        # Directions calculated taking a unique corner as reference (difference to remaining corners)
        for i in range(total_unique_corners-1):
            polygon_direct[i,:] = unique_corners[i+1,:] - unique_corners[0,:]

        # Rank
        rank = matrix_rank(polygon_direct)

        # Conclusion on full dimensionality
        if (rank == dim['n']):
            full_dim_CR_bool = True
        else:
            full_dim_CR_bool = False

    else:
        full_dim_CR_bool = False

    # Calculate clockwise order of unique vertices - to ensure they are correctly displayed as a full polygon
    # (this function is applicable only to 2D problems)
    if ((dim['n'] == 2) and alg_par['display_2D_bool'] and full_dim_CR_bool and order_bool):
        order = argsort(arctan2(unique_corners[:,0] - average(unique_corners[:,0]), unique_corners[:,1] - average(unique_corners[:,1])))
        return full_dim_CR_bool, unique_corners[order]

    else:
        return full_dim_CR_bool, unique_corners

#%% This is an au auxiliary function, which facilitates at a later stage to calculate the distance between any x and
# any of the corners of the feasible space.
def corners_coords_square_sum(dim, feas_space_corners):

    # Initialisation
    all_corners_square_sums = zeros(dim['t'])

    for i in range(dim['t']):
        all_corners_square_sums[i] = dot(feas_space_corners[i], feas_space_corners[i])

    return all_corners_square_sums

#%% This function takes a single constraint of the feasible space as the key argument, and calculates the coordinates of
# all corners defined via the interception of the "exterior" halfspace and the bounds of the display (applicable only
# to 2D). This information can then be used to mask the display of plots, to emphasize the feasible space only.
def polygon_mask_2D(feas_space, dim, alg_par, graph_display_settings_2D, input_index):

    # Initialisation: mask corners (empty) and inequality checks per constraint
    raw_mask_corners        = zeros((0,2))
    display_corners_check   = zeros(4)      # for corners 'c1', 'c2', 'c3' and 'c4'

    # Corners testing
    for i in range(4):
        display_corners_check[i] = dot(-feas_space['A'][input_index,:], graph_display_settings_2D['corners'][i]) + feas_space['b'][input_index]
        if (display_corners_check[i] <= 0):
            raw_mask_corners = append(raw_mask_corners, graph_display_settings_2D['corners'][i].reshape([1,dim['n']]), axis=0)

    ###########################################################################
    # Check min and max x1 bounds (x1 = lb, x1 = ub)
    ###########################################################################
    # Both cases require a similar treatment, since the same matrix is used in them. The candidate points in both cases
    # are calculated as the interception of the given feasible space constraint and [1, 0] (for x1 >= lb, x1 <= ub)
    interception_matrix = array([feas_space['A'][input_index,:], array([1,0])])
    rank                = matrix_rank(interception_matrix)

    # Only non-parallel lines yield a point - calculations proceed only if this criterium is satisfied.
    if (rank == dim['n']):

        # Inverse matrix
        matrix_inv = inv(interception_matrix)

        # Check at x1 = lb: calculate point and check if it is in the bounds of the display; if so add point 
        interception_vector = array([feas_space['b'][input_index], graph_display_settings_2D['min_x1']])
        candidate_corner    = dot(matrix_inv, interception_vector)
        raw_mask_corners    = mask_corners_test_update(dim, candidate_corner, raw_mask_corners, graph_display_settings_2D)

        # Check at x1 = ub: calculate point and check if it is in the bounds of the display; if so add point 
        interception_vector = array([feas_space['b'][input_index], graph_display_settings_2D['max_x1']])
        candidate_corner    = dot(matrix_inv, interception_vector)
        raw_mask_corners    = mask_corners_test_update(dim, candidate_corner, raw_mask_corners, graph_display_settings_2D)

    ###########################################################################
    # Check min and max x2 bounds (x2 = lb, x2 = ub)
    ###########################################################################
    # Both cases require a similar treatment, since the same matrix is used in them. The candidate points in both cases
    # are calculated as the interception of the given feasible space constraint and [0, 1] (for x2 = lb, x2 = ub)
    interception_matrix = array([feas_space['A'][input_index,:], array([0,1])])
    rank                = matrix_rank(interception_matrix)

    # Only non-parallel lines yield a point - calculations proceed only if this criterium is satisfied.
    if (rank == dim['n']):

        # Inverse matrix
        matrix_inv = inv(interception_matrix)

        # Check at x2 = lb: calculate point and check if it is in the bounds of the display; if so add point 
        interception_vector = array([feas_space['b'][input_index], graph_display_settings_2D['min_x2']])
        candidate_corner    = dot(matrix_inv, interception_vector)
        raw_mask_corners    = mask_corners_test_update(dim, candidate_corner, raw_mask_corners, graph_display_settings_2D)

        # Check at x2 = ub: calculate point and check if it is in the bounds of the display; if so add point 
        interception_vector = array([feas_space['b'][input_index], graph_display_settings_2D['max_x2']])
        candidate_corner    = dot(matrix_inv, interception_vector)
        raw_mask_corners    = mask_corners_test_update(dim, candidate_corner, raw_mask_corners, graph_display_settings_2D)

    # Remove duplicates (if any), check they define a full-dimensional region, and update their order
    # for correct display (in 2D only)
    full_dim_CR_bool, mask_corners = corners_pre_proc(dim, raw_mask_corners, alg_par, True)

    return mask_corners

#%% This function takes an array of coordinates for polygon representation, and a candidate vector of coordinates that,
# if included in the bounds of the 2D display is also added to this array.
def mask_corners_test_update(dim, point, mask_corners, graph_display_settings_2D):

    if (point[0] >= graph_display_settings_2D['min_x1'] and point[0] <= graph_display_settings_2D['max_x1'] and point[1] >= graph_display_settings_2D['min_x2'] and point[1] <= graph_display_settings_2D['max_x2']):
        mask_corners = append(mask_corners, point.reshape([1,dim['n']]), axis=0)

    return mask_corners

#%%############################################################################
# Optimization functions
###############################################################################
#%% Linear programming (LP), including inequality constraints only. Check if a set of bounds of the form Ax <= b
# return a feasible solution.
def LP_ineq(cost, A_ub, b_ub, bounds):

    status = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs').status

    return status

#%% LP, including inequality and equality constraints. Check if a set of bounds of the form Ax <= b and a subset
# of the matching equality constraints Ax = b return a feasible solution (feasibility and optimal solution).
def LP_ineq_eq(cost, A_ub, b_ub, A_eq, b_eq, bounds):

    sol = linprog(cost, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

    status  = sol.status
    x_opt   = sol.x

    return status, x_opt
