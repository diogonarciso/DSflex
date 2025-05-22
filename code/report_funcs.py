from pandas import DataFrame, ExcelWriter

#%%############################################################################
# Create DataFrame
###############################################################################
#%% DataFrame from matrix.
def df_calc(mat_input):

    df = DataFrame(mat_input)

    return df

#%%############################################################################
# Export to Excel
###############################################################################
#%% Auxiliary function - save multiple dataframes to Excel.
def save_xls(list_dfs, sheet_list, xls_path):
    with ExcelWriter(xls_path) as writer:
        for n, df in enumerate(list_dfs):
            df.to_excel(excel_writer=writer, sheet_name=sheet_list[n])
