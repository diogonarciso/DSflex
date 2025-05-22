# DSflex
Introduction:

This project contains the Python files delivering the functionality for flexibility assessment, as presented in XX.
This page provides a brief introduction to setting up the system, configuring examples, and obtaining and interpreting results.
For any comments, please contact diogo.narciso@tecnico.ulisboa.pt.

System set-up:
1) Download and install Python. Spyder was used to develop and test the code files in this project; any other preferred interface can deliver the desired functionality.
2) Download the code files in this project to a subfolder "code" in the selected root directory.
3) Open the "main.py" file and edit line 4 to identify the root variable (as in step 2). Note that forward slashes should be used in paths definition (e.g. 'C:/Users/User1/...')
4) Create a subfolder "examples" in the root folder, where all the examples to be solved should be included. Some examples can be found in the article' Supplementary Materials.
5) Create a subfolder "results" in the root folder, where all results will be saved.
6) This project makes use of libraries "pandas", "numpy", "scipy", "matplotlib", which may be installed by default in the Python working environment.
7) If any of these libraries is missing, please install them (an error message will be displayed when attempting to run the main functions if this is the case).
8) Incompatibilities between versions of these libraries are not expected, but for reference, the versions used in code development are the following: pandas (2.2.2), numpy (1.26.4), scipy (1.13.1), matplotlib (3.8.4).

Examples configuration:
1) The bounds of design spaces are defined as Ax <= b, where x denotes the vector of design variables and A and b are a matrix and vector, respectively, defining the inequalities bounding the design space.
2) All examples should be defined via an Excel file, with information on A and b included in a sheet called "feas_space".
3) The first row of this worksheet should list "b" and the columns of A (e.g., "A_1"—weights for the first design variable, "A_2"—weights for the second design variable, ...) in sequence from the second row of the worksheet.
4) The inequalities should then be listed from the second row of the worksheet, with their index in the first column and the corresponding weights in the remaining columns.
5) For reference, check the example files provided (Step 4 in system set-up).

Using the code:
1) Seven code files enable the full functionality of the 5 algorithms in the reference article.
2) All code files except the "main.py" and "par.py" files should not be edited unless some additional development is required. Executing the functions in the first file delivers flexibility results, and some options for running them may be set in the second file.
3) Given an example in the form Ax <= b, it must be saved as an Excel file in the "examples" subfolder and its name specified accordingly in line 9 of the "main.py" file (e.g. "ex1"; without the ".xlsx" file extension).
4) Then, two options are available to obtain results: (i) execute the first 2 cells of the "main.py" file, or (ii) all cells in sequence except cell 2.
5) In the first option, the relevant paths are first set-up (cell 1), and then all flexibility results are calculated (cell 2).
6) Critical regions and flexibility maps are plotted, and 2 Excel files are (optionally) saved in the "results" subfolder.
7) In the second option, the user must execute all cells in sequence (except cell 2). This enables a step-by-step execution of Algorithms 1-3 (cell 4), Algorithm 4 (cell 6), and Algorithm 5 (cell 8) and selectively inspecting the solutions obtained (cells 5, 7, and 9).

Excel results:
1) If this feature is enabled, two Excel files are created per example solved.
2) The first file includes the results relating to the partial minimum distances ("file_name_min.xlsx"), and the second the results relating to the partial maximum distances ("file_name_max.xlsx").
3) All results files include several sheets, each presenting information on one of the flexibility metrics created (e.g. minimum distance functions). In all sheets, except "avg_min_dist_func" and "avg_max_dist_func", the first n columns represent the dependence on theta, and the last row is the independent term.
4) In the sheets identified in Step 3, the dependence on theta of average distance functions is represented in the first n rows and the independent term is defined in the last row.
5) In sheets "min_dist_funcs" and "max_dist_funcs", a single distance funtion per the corresponding bound/corner are listed.
6) In the remaining sheets, results from numpy are 3D arrays, which are transformed into 2D arrays, showing the results per bound/corner in sequence. Functions "single_constr_df_calc" and "n_constr_df_calc" control how these arrays are reshaped.

Options:
1) The "par.py" file allows enabling/disabling several options and editing several parameters related to the display of 2D plots.
2) By default, parameter "min_calc_mode" is set to "dist_only". If set to "projection", not only the minimum distance to each bound is delivered, but also explicit functions for the projection (optimizer) on the bounds are calculated.
3) As discussed in the article, a preliminary assessment shows that the max-min and min-max distances are not particularly relevant for flexibility assessment, and their calculation is disbled by default. To check results for these metrics, it suffices to set parameters "max_min_calc" and/or "min_max_calc" to True.
4) The display of critical regions and flexibility maps (2D only) is enabled/disabled by parameter "display_2D_bool", and saving results to Excel may be enabled via parameter "save_results".
5) Several other options for displaying 2D plots are also available.
