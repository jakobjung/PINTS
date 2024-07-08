# Script that checks the csv file for the correct format
# Usage: python check_csv.py <csv_file>
# Example: python check_csv.py data.csv

# import libraries
import numpy as np
import pandas as pd
import sys

# set errors to false
errors = False

# read in the csv file
csv_file = sys.argv[1]

# read in the csv file
df = pd.read_csv(csv_file)

# print column names
print("Column names:")
print(df.columns)

# convert column names to lower case
df.columns = [x.lower() for x in df.columns]

# retain on





# check whether genome	type	start	end	strand	p_value	gene_name are present in the csv file.
if 'genome' in df.columns and 'type' in df.columns and 'start' in df.columns and 'end' in df.columns and\
        'strand' in df.columns and 'p_value' in df.columns and 'gene_name' in df.columns:
    print("All required columns are present")
    # change p_value column name to p_value_[filename]
    # get filename without .csv and path
    csv_filename = csv_file.split('/')[-1].replace('.csv', '')
    print(csv_filename)
    df = df.rename(columns={'p_value': 'p_value_' + csv_filename})
    # retain only "CDS" and "sRNA" in the "type" column
    df = df[df['type'].isin(["CDS", "sRNA"])]
    # print all types
    print("All types:")
    print(df['type'].unique())
    # save the csv file
    df.to_csv(csv_file, index=False)
else:
    # write error message with name of the csv file
    error_message = "Error: Missing required columns in the csv file: " + csv_file + "\n" + \
                    "Required columns are: genome, type, start, end, strand, p_value, gene_name"
    errors = True

    # raise an error
    raise ValueError(error_message)







