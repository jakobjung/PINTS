#calculate pvalues of MFEs
import pandas as pd
from read_intarna import read_IntaRNA_out_summarize

#function for calculating the MFE pvalue empirically
def get_stat_pval(input_shuff, input_org):
    # Merge the two input tables
    temp = pd.merge(input_shuff, input_org, 
                                  how="left", on = 'geneID', 
                                  suffixes = ('_shuff', '_org'))
    
    # Get the stats of the shuffled MFE values for each gene
    MFE_org_shuff_stats = temp.groupby(['geneID', 'MFE_org']).MFE_shuff.describe().reset_index()
    
    # Get the p-value for the observed MFE-org of each geneID
    MFE_org_shuff_pval = temp.groupby('geneID').apply(lambda x: len(x[x['MFE_shuff']<=x['MFE_org']])/len(x)).reset_index(name='pval')
    
    # Now put the stats and the p values in one table
    MFE_stats_pval = pd.merge(MFE_org_shuff_stats, MFE_org_shuff_pval, how="left", on = 'geneID')
    
    # Delete what we don't need
    MFE_org_shuff_stats = None
    MFE_org_shuff_pval = None
    
    
    return MFE_stats_pval
    

def calc_pval(PinT_SPI1_MFE):
	SPI1_MFE_shuff1 = read_IntaRNA_out_summarize('../new_saves/SPI1_PinT_shuffled_256_out.csv')
	print (SPI1_MFE_shuff1.head(1))
	SPI1_MFE_shuff1["geneID"] = SPI1_MFE_shuff1["geneID"].replace("-shuffled-\d+", "", regex=True)
	# Again let's see what the MFE values look like
	print (SPI1_MFE_shuff1.describe())
	
	SPI1_MFE_shuff1 = pd.concat([SPI1_MFE_shuff1, PinT_SPI1_MFE])
	print (SPI1_MFE_shuff1.shape)
	
	# Use the get_stat_pval to summarize all the info regrading the first 256 shuffles
	SPI1_MFE_256 = get_stat_pval(SPI1_MFE_shuff1, PinT_SPI1_MFE)
	print (SPI1_MFE_256.head())
	
	
	
#NOTE:
#repeat for all the n-shuffles, ie, 256x32 shuffled sequences
#concatenate them, summarize, calculate pvalues

#use list.dir to get all shuffled intarna files
#also have to save each of the files with p values
