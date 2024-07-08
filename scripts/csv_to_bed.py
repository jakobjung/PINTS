# script that converts the csv file to a bed file.
# Usage: python csv_to_bed.py <csv_file>

# import libraries
import numpy as np
import pandas as pd
import sys

# read in the csv file
csv_file = sys.argv[1]

# read in the csv file
df = pd.read_csv(csv_file)

# convert column names to lower case
df.columns = [x.lower() for x in df.columns]

# put the columns in the correct order for bed file
df = df[['genome', 'start', 'end', 'gene_name', 'type', 'strand']]
print(df)

# if type is "CDS" or "gene" and the strand is "+" , change start to start-50 and end to start+50
df.loc[(df['type'] == 'CDS') & (df['strand'] == '+'), 'start'] = df['start'] - 50
df.loc[(df['type'] == 'CDS') & (df['strand'] == '+'), 'end'] = df['start'] + 100

# if type is "CDS" or "gene" and the strand is "-" , change start to end-50 and end to end+50
df.loc[(df['type'] == 'CDS') & (df['strand'] == '-'), 'start'] = df['end'] - 51
df.loc[(df['type'] == 'CDS') & (df['strand'] == '-'), 'end'] = df['start'] + 100

print (df)

# save the bed file without header
filename_bed = csv_file.replace('.csv', '.bed')
df.to_csv(filename_bed, sep='\t', header=False, index=False)
