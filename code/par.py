# Parameters for algorithm execution
min_calc_mode           = 'dist_only'   # 'dist_only' mode delivers only the distance functions; 'projection' mode delivers also the optimizer functions
max_min_calc            = False         # min-max flexibility calculated only if this bool is set to True
min_max_calc            = False         # max-min flexibility calculated only if this bool is set to True
display_2D_bool         = True          # If True, displays critical regions and heatmaps
save_results            = False         # If True, saves results to Excel
large_pos_number        = 1e10          # Lower and upper bounds for LPs are set from this parameter

alg_par = {'min_calc_mode':min_calc_mode, 'max_min_calc':max_min_calc, 'min_max_calc':min_max_calc,
           'display_2D_bool':display_2D_bool, 'save_results':save_results, 'large_pos_number':large_pos_number}

# Parameters for 2D graphs
N_display_points    = 100   # Number of points per axis for delivering the distance heatmaps
edge_gap            = 0.5   # Controls the gap between the minimum/maximum values of x1 and x2 in the design space and the corresponding display bounds
fig_width           = 4     # Options for figures
fig_height          = 4
font_size_title     = 20
font_size_axis      = 15    # Axes: ticks and labels
scatter_size        = 200   # Size of corners (if enabled)
label_size          = 12    # Colorbar (font sizes of scale)
display_corners     = True  # Highlights corners for the corresponding critical regions in the same colour if set to True

# Colours used in the representation of critical regions (20 colours); depending on the number of critical regions, some colours may be repeated
unique_colours      = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black',
                       'rosybrown', 'coral', 'sienna', 'burlywood', 'navajowhite', 'gold', 'honeydew', 'turquoise', 'lightsteelblue']

graph_display_settings_2D = {'N':N_display_points, 'edge_gap':edge_gap, 'fig_width':fig_width, 'fig_height':fig_height, 'scatter_size':scatter_size, 
                             'font_size_title':font_size_title, 'font_size_axis':font_size_axis, 'label_size':label_size, 'display_corners':display_corners}