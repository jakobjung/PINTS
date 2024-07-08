#reading IntaRNA outputs

import pandas as pd

def read_IntaRNA_out_summarize(file):
    # Read the output of intaRNA
    temp = pd.read_csv(file, sep = ";", dtype = {'E': 'float32'})
    
    # Make a new column that uses column id1 to get the geneID only
    temp['geneID'] = temp['id1'].str.split('::').str[0]
    
     # Pick columns of interest
    temp = temp[['geneID', 'E']]
    
    # Rename E to MFE
    temp = temp.rename(columns = {'E': 'MFE'})
    #temp = temp.rename(columns ={'id1':'Gene name'})
    
    # Change data type to make smaller
    #temp=temp.astype({'MFE': 'float32'})
    
    return temp
    
def called():
	# Read and summarize the MFE values for all of the genes 
	print ("reading intarna outs and testing...")
	PinT_SPI1_MFE = read_IntaRNA_out_summarize('SPI1_PinToutMFEs.csv')
	print (PinT_SPI1_MFE.min())
	return PinT_SPI1_MFE
