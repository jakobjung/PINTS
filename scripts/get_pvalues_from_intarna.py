# script that creates p-values from the intarna results

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gumbel_l
import seaborn as sns
import sys


# make file1 and file2 as arguments (added to python call)
file_intarna = sys.argv[1]
file_intarna_shuff = sys.argv[2]

# read in the intarna results. delim is ;
intarna_results = pd.read_csv(file_intarna, sep=";")
# read in shuffled intarna results
intarna_results_shuff = pd.read_csv(file_intarna_shuff, sep=";")

# # read in the intarna results. delim is ;
# intarna_results = pd.read_csv("../data/MAPS_test_intarna.csv", sep=";")
#
# # read in shuffled intarna results
# intarna_results_shuff = pd.read_csv("../data/MAPS_test_shuffled1_intarna.csv", sep=";")


# save plot as png
# initiate new plot
plt.figure()
# plot distribution of MFE values.use sns
sns.histplot(intarna_results["E"], bins=50)
plt.savefig("./analysis/IntaRNA_MFE_dist.png")
# close plot
plt.close()

# do same with shuffled
plt.figure()
sns.histplot(intarna_results_shuff["E"], bins=50)
plt.savefig("./analysis/IntaRNA_MFE_dist_shuffled.png")
plt.close()


# try doing same but fit a gumbel distribution to the shuffled MFE values
# fit the gumbel distribution
shuff_mfe = intarna_results_shuff["E"]
params = gumbel_l.fit(shuff_mfe)
loc, scale = params


# plot the fitted gumbel distribution
plt.figure()
sns.histplot(shuff_mfe, bins=50,  stat="density")
#plt.hist(shuff_mfe, bins=50, density=True)
x = np.linspace(min(shuff_mfe), max(shuff_mfe), 1000)
# plot the fitted gumbel distribution. color is orange. use sns
plt.plot(x, gumbel_l.pdf(x, loc=loc, scale=scale), color="darkorange")
# change x axis title
plt.xlabel("Minimum free energy (MFE) in ΔG")
# remove y axis title
plt.ylabel("")
# add legend
plt.legend(["Gumbel fit", "Density"])
plt.savefig("./analysis/IntaRNA_MFE_dist_gumbel.png")
plt.close()

# calculate the cdf
cdf = gumbel_l.cdf(intarna_results["E"], loc=loc, scale=scale)

# calculate the p-values
intarna_results["p_value_IntaRNA"] = cdf

# plot the distribution of p-values
plt.figure()
sns.histplot(intarna_results["p_value_IntaRNA"], bins=30)
plt.savefig("./analysis/IntaRNA_pvalue_dist.png")
plt.close()


# extract the gene name from intarna results
intarna_results["gene_name"] = intarna_results["id1"].str.split("::", expand=True)[0]

# only pick columns gene_name, E and p_value
intarna_results = intarna_results[["gene_name", "E", "p_value_IntaRNA"]]

# rename E to MFE(Gibbs free energy)
intarna_results = intarna_results.rename(columns={"E": "MFE(Gibbs free energy)"})


# save the results
intarna_results.to_csv("./data/MAPS_test_intarna_pvalues.csv", index=False)





