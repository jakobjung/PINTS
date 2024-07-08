import pandas as pd
import numpy as np
from Bio import SeqIO 
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import seaborn as sns
from scipy import stats
from statsmodels.stats import multitest as mt 

#hardcoded links - may not be a problem if the files are already created
#through the previous steps
#
#TODO: plotting MFEs/p-vals vs shuffles for a bunch of genes
#can users choose?

#doesn't contain user-chosen criteria for comparing p-values

#loading Pulse
def load_Pulse():
	#OR pass file by value when calling the fn
	#the filename can be replaced w/ an object
	Pulse_FC_SPI1 = pd.read_csv('SPI1_Pulse_Tidy.csv')
	Pulse_FC_SPI1 = Pulse_FC_SPI1.rename(columns={'logFC':'Pulse_logFC','PValue':'Pulse_pValue', 'FDR':'Pulse_FDR'})
	print(f'{Pulse_FC_SPI1.head()}')
	return Pulse_FC_SPI1
	
def load_MAPS():
	#the filename can be replaced w/ an object
	MAPS_FC_SPI1_2sided = pd.read_csv('SPI1_MAPS_Tidy.csv')
	MAPS_FC_SPI1_2sided = MAPS_FC_SPI1_2sided[['Gene_name','Locus_tag','MS2-PinT_SPI-1_vs_PinT_SPI-1_log2FoldChange','MS2-PinT_SPI-1_vs_PinT_SPI-1_pvalue','MS2-PinT_SPI-1_vs_PinT_SPI-1_padj']]
	MAPS_FC_SPI1_2sided = MAPS_FC_SPI1_2sided.rename(columns={'MS2-PinT_SPI-1_vs_PinT_SPI-1_log2FoldChange':'MAPS_logFC','MS2-PinT_SPI-1_vs_PinT_SPI-1_pvalue':'MAPS_pValue', 'MS2-PinT_SPI-1_vs_PinT_SPI-1_padj':'MAPS_FDR'})
	#convert p-val to one-tailed
	MAPS_FC_SPI1 = MAPS_FC_SPI1_2sided.copy()
	MAPS_FC_SPI1['MAPS_pValue'] = MAPS_FC_SPI1['MAPS_pValue']*0.5
	print (MAPS_FC_SPI1.head())
	return MAPS_FC_SPI1
	
def merge_tables(MAPS_FC_SPI1, Pulse_FC_SPI1):
	MAPS_Pulse_SPI1_commongenes = pd.merge(MAPS_FC_SPI1, Pulse_FC_SPI1, 
                       how='inner',
                       on=['Gene_name'])
	return MAPS_Pulse_SPI1_commongenes      

def load_MFE(MAPS_Pulse_SPI1_commongenes):
	#, Pulse_FC_SPI1, MAPS_FC_SPI1):
	#the filename can be replaced w/ an object
	
	MFE_SPI1_commongenes= pd.read_csv('SPI1_MFE_pval_tidy_8192.csv')
	MFE_SPI1_commongenes=MFE_SPI1_commongenes.rename(columns={'geneID':'Gene_name','MFE_org':'MFE','pval':'MFE_pValue'})
	#
	#unnecessary
	#{
	#flag = 0 #tracer to check mismatches
	#length of Pulse dataset
	#len_pulse = len(Pulse_FC_SPI1)
	#number of genes
	#genes_pulse = Pulse_FC_SPI1['Gene_name'].nunique()		
	
	#length of MAPS dataset
	#len_maps = len(MAPS_FC_SPI1)
	#number of genes
	#genes_maps = MAPS_FC_SPI1['Gene_name'].nunique()
	
	#length of MFE dataset
	#len_MFE = len(MFE_SPI1_commongenes)
	#number of genes
	#genes_MFE = MFE_SPI1_commongenes['Gene_name'].nunique()
	#if not len_MFE == genes_MFE:
	#	flag = 1

	#length of MAPS+Pulse common genes dataset
	#this can be used to check mismatches in lengths and gene counts
	#len_MP = len(MAPS_Pulse_SPI1_commongenes)
	#genes_MP = MAPS_Pulse_SPI1_commongenes['Gene_name'].nunique()
	#if not len_MP == genes_MP:
	#	flag = 1
	#}
	#
	
	#if flag = 1, then merge the datasets based on common genes
	#unnecessary, though
	#since the datasets need to be merged anyway. we can simply get rid of duplicates
	SPI1_combined_V2 = pd.merge(MAPS_Pulse_SPI1_commongenes, MFE_SPI1_commongenes, how='inner', on='Gene_name')

	SPI1_combined_V2.drop_duplicates(subset ='Gene_name',
                     keep = False, inplace = True)
	print('Length of combined_df: ', len(SPI1_combined_V2))
	print('Number of genes in combined_df: ', SPI1_combined_V2['Gene_name'].nunique())
	return (SPI1_combined_V2)

def plot_pval_dist(SPI1_combined_V2):
	# Create the figure for distribution of p-values
	fig, ax = plt.subplots(ncols=3, figsize=(20,5))


	# Create the titles for each subplot
	ax[0].set_xlabel('pVal of log2FC_Pulse data')
	ax[1].set_xlabel('one-sided pVal of log2FC_MAPS data')
	ax[2].set_xlabel('one-sided pVal of MFE values')
	
	# Plot
	SPI1_combined_V2['Pulse_pValue'].hist(bins =100, density = False, ax=ax[0])
	SPI1_combined_V2['MAPS_pValue'].hist(bins =100, density = False, ax=ax[1])
	SPI1_combined_V2['MFE_pValue'].hist(bins =100, density = False, ax=ax[2])
	plt.savefig("../plots/pvalue_MFEs.png")#of course can be changed to vector graphics

def load_shuffles():
	stats_256 = pd.read_csv('new_saves/SPI1_MFE_256.csv')
	stats_512 = pd.read_csv('new_saves/SPI1_MFE_512.csv')
	stats_1024 = pd.read_csv('new_saves/SPI1_MFE_1024.csv')
	stats_2048 = pd.read_csv('new_saves/SPI1_MFE_2048.csv')
	Stats = pd.concat([stats_256, stats_512, 
                   stats_1024, stats_2048])
#TODO: plotting MFEs/p-vals vs shuffles for a bunch of genes
#can users choose?

def combine_pval(SPI1_combined_V2):
	SPI1_combined_V2['combi_3pvals']=SPI1_combined_V2.apply(lambda x: stats.combine_pvalues([x['MAPS_pValue'], x['MFE_pValue'], x['Pulse_pValue']], 
                                                                      method='fisher', 
                                                                      weights=None)[1], axis =1)

	##p-val histograms##
	fig, ax = plt.subplots( figsize=(10,5))
	ax.set_xlabel('Combination of all 3 p-vals by Fisher method')
	SPI1_combined_V2['combi_3pvals'].hist(bins =100, density = False, ax=ax)
	plt.savefig("../plots/pvalue_fisher.png")
	
	##multiple testing correction##
	SPI1_combined_V3 = SPI1_combined_V2.sort_values(by = ['combi_3pvals'])
	SPI1_combined_V3['FDR3'] = mt.fdrcorrection(SPI1_combined_V3['combi_3pvals'], alpha=0.05)[1]
	#SPI1_combined_V3.to_csv('../new_saves/SPI1_pvals.csv', index=False)
	
	##plot##
	#describe
	SPI1_combined_V3.drop(SPI1_combined_V3.loc[SPI1_combined_V3['Gene_name']=='STnc440'].index, inplace=True)
	#MAPS vs Pulse p-vals
	fig, ax = plt.subplots(figsize=(20,10))
	SPI1_combined_V3.plot.scatter(x = 'MAPS_pValue', y = 'Pulse_pValue', s = 25, c = 'combi_3pvals', colormap = 'jet', ax =ax)
	plt.savefig("../plots/pvalues_maps_vs_pulse_allgenes.png")
	
	##significant genes p<0.05##
	temp =SPI1_combined_V3[SPI1_combined_V3['combi_3pvals']<0.05]

	fig, ax = plt.subplots(figsize=(20,10))

	SPI1_combined_V3.plot.scatter(x = 'MAPS_pValue', y = 'Pulse_pValue', s = 25, c = 'grey', ax =ax)
	temp.plot.scatter(x = 'MAPS_pValue', y = 'Pulse_pValue', s = 25, c = 'combi_3pvals', colormap = 'Reds_r', ax =ax)
	plt.savefig("../plots/pvalues_maps_vs_pulse_lt0.5.png")
	
	return SPI1_combined_V3

def identify_genes(SPI1_combined_V3):
	#combined p-val < 0.05
	#log2FC > 0 (MAPS)
	SPI1_filtered_2 = SPI1_combined_V3[(SPI1_combined_V3['FDR3']<0.05) & (SPI1_combined_V3['MAPS_logFC']>0)]
	print('Number of genes that have combi_3pvals < 0.05 and MAPS_logFC > 0 :', len(SPI1_filtered_2))
	
	#Filter 2
	#where do these genes come from?
	True_gene_names = ['SL1344_4251', 'hilA', 'sopE', 'grxA','rpsV','fliC','sopE2','ecnB','ugtL','rtsA','grxA','SteC']
	True_targets=SPI1_combined_V3[SPI1_combined_V3['Gene_name'].isin(True_gene_names)]
	#Filter 3
	Hits_in_both=SPI1_filtered_2[SPI1_filtered_2['Gene_name'].isin(True_gene_names)]
	
	##plot p-values##
	##Pulse, MFE##
	
	fig, ax = plt.subplots(figsize=(15,10))
	ax.set_yscale('log')
	ax.invert_yaxis()

	SPI1_combined_V3.plot.scatter(x = 'MFE_pValue', y = 'Pulse_pValue', s = 20, c = 'gainsboro', label = 'all', ax =ax)
	SPI1_filtered_2.plot.scatter(x = 'MFE_pValue', y = 'Pulse_pValue', s = 20, c = 'darkgrey', label = 'SPI1_filtered_2', ax =ax)

	# To further color dots representing genes with certain properties
	temp=None
	temp = SPI1_filtered_2[SPI1_filtered_2['FDR3']<1e-7]
	temp.plot.scatter(x = 'MFE_pValue', y = 'Pulse_pValue', s = 20, c = 'mediumturquoise', label = 'combined_FDR<1e-7', ax =ax)

	True_targets.plot.scatter(x = 'MFE_pValue', y = 'Pulse_pValue', s = 20, c = 'forestgreen', label = 'True_targets', ax =ax)

	temp2=Hits_in_both[Hits_in_both['FDR3']<1e-7]
	temp2.plot.scatter(x = 'MFE_pValue', y = 'Pulse_pValue', s = 20, c = 'tomato', label = 'Hits in both',  ax =ax)

	#adding gene names
	True_targets[['MFE_pValue','Pulse_pValue','Gene_name']].apply(lambda row: ax.text(*row, fontsize=13, fontweight='normal'),axis=1)
	plt.savefig("../plots/Pulse_MFE_sig.png")
	
	##MAPS, MFE##
	fig, ax = plt.subplots(figsize=(15,10))
	ax.set_yscale('log')
	ax.invert_yaxis()

	SPI1_combined_V3.plot.scatter(x = 'MFE_pValue', y = 'MAPS_pValue', s = 20, c = 'gainsboro', label = 'all', ax =ax)
	SPI1_filtered_2.plot.scatter(x = 'MFE_pValue', y = 'MAPS_pValue', s = 20, c = 'darkgrey', label = 'SPI1_filtered_2', ax =ax)


	# To further color dots representing genes with certain properties
	temp=None
	temp = SPI1_filtered_2[abs(SPI1_filtered_2['FDR3'])<1e-7]
	temp.plot.scatter(x = 'MFE_pValue', y = 'MAPS_pValue', s = 20, c = 'mediumturquoise', label = 'combined_FDR<1e-7', ax =ax)

	True_targets.plot.scatter(x = 'MFE_pValue', y = 'MAPS_pValue', s = 20, c = 'forestgreen', label = 'True targets', ax =ax)

	temp2=None
	temp2=Hits_in_both[abs(Hits_in_both['FDR3'])<1e-7]
	temp2.plot.scatter(x = 'MFE_pValue', y = 'MAPS_pValue', s = 20, c = 'tomato', label = 'Hits in both',  ax =ax)

	True_targets[['MFE_pValue','MAPS_pValue','Gene_name']].apply(lambda row: ax.text(*row, fontsize=13, fontweight='normal'),axis=1)
	plt.savefig("../plots/MAPS_MFE_sig.png")
	
	
	##MAPS, Pulse##
	
	fig, ax = plt.subplots(figsize=(15,10))
	ax.set_yscale('log')
	ax.invert_yaxis()

	SPI1_combined_V3.plot.scatter(x = 'Pulse_pValue', y = 'MAPS_pValue', s = 20, c = 'gainsboro', label = 'all', ax =ax)
	SPI1_filtered_2.plot.scatter(x = 'Pulse_pValue', y = 'MAPS_pValue', s = 20, c = 'darkgrey', label = 'SPI1_filtered_2', ax =ax)

	# To further color dots representing genes with certain properties
	temp=None
	temp = SPI1_filtered_2[SPI1_filtered_2['FDR3']<1e-7]
	temp.plot.scatter(x = 'Pulse_pValue', y = 'MAPS_pValue', s = 20, c = 'mediumturquoise', label = 'combined_FDR<1e-7', ax =ax)

	True_targets.plot.scatter(x = 'Pulse_pValue', y = 'MAPS_pValue', s = 20, c = 'forestgreen', label = 'True targets', ax =ax)
	temp2=None
	temp2=Hits_in_both[Hits_in_both['FDR3']<1e-7]
	temp2.plot.scatter(x = 'Pulse_pValue', y = 'MAPS_pValue', s = 20, c = 'tomato', label = 'Hits in both',  ax =ax)

	True_targets[['Pulse_pValue','MAPS_pValue','Gene_name']].apply(lambda row: ax.text(*row, fontsize=13, fontweight='normal'),axis=1)
	plt.savefig("../plots/MAPS_Pulse_sig.png")	

def start():
	MAPS_FC_SPI1 = load_MAPS() #either let the fn load the file or pass by value
	Pulse_FC_SPI1 = load_Pulse()#ditto
	MAPS_Pulse_SPI1_commongenes = merge_tables(MAPS_FC_SPI1, Pulse_FC_SPI1)
	#get combined datasets w/ MFEs
	SPI1_combined_V2 = load_MFE(MAPS_Pulse_SPI1_commongenes)
	plot_pval_dist(SPI1_combined_V2)
	#loading shuffles
	load_shuffles()
	#combine p-values, apply correction
	SPI1_combined_V3 = combine_pval(SPI1_combined_V2)
	#identifying significant genes
	identify_genes(SPI1_combined_V3)
	
