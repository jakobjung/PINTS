# script that combines p-values from all the input files using Fisher's method
#
# import libraries
import numpy as np
import pandas as pd
import sys
from scipy.stats import chi2
from functools import reduce
from scipy import stats

# read in the csv files
# csv_files = sys.argv[1:]
csv_files = ["../data/MAPS_test_intarna_pvalues.csv", "../data/MAPS_test.csv"]


# check how many files are there
print("Number of files: ", len(csv_files))

# read in the csv files
dfs = [pd.read_csv(file) for file in csv_files]

# combine the input files by a columns called "gene_name"
# use reduce to merge the dataframes. I use reduce because we don't know how many dataframes we have
df_final = reduce(lambda left, right: pd.merge(left, right, on='gene_name', how='outer'), dfs)

# combine the p-values using Fisher's method
# get the p-values & gene_name columns
p_values_df = df_final.filter(regex='p_value_')
# add gene name column
p_values_df['gene_name'] = df_final['gene_name']
# for NA in p_value_intaRNA, add 1
p_values_df.fillna(1, inplace=True)



# combine the p-values using Fisher's method (do it for each row)
p_values_df["combined_p_value"] = p_values_df.filter(regex='p_value_').apply(lambda x: stats.combine_pvalues(x, method='fisher')[1], axis=1)


# write it as a csv file
p_values_df.to_csv("../data/MAPS_test_noncombined_pvalues.csv")







