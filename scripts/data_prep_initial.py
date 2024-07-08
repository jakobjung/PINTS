#load data from MAPS and Pulse expression tables
#initial script

import pandas as pd
import numpy as np
from Bio import SeqIO 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats import multitest as mt 
import os

from read_intarna import called
from shuffle_seqs import shuffle_nts
from subprocess import call
from calc_pval import calc_pval


#MAPS
def load_MAPS(maps):
#the param passed is a list variable
#if there are multiple datasets it can be looped through
#and each dataset can be loaded
#
#for now just the first element is retrieved
#it's in start2.py, only the first element is passed on
	print ("load MAPS")
	#MAPS_SPI1_PinT= pd.read_csv('../DataFiles_PinT/MAPS/MAPS_SPI1_PinT.csv')
	MAPS_SPI1_PinT= pd.read_csv(maps)
	print (f"{MAPS_SPI1_PinT.head()}")

	# Pick only iteresting columns
	MAPS_SPI1_PinT_v1 = MAPS_SPI1_PinT.copy()
	MAPS_SPI1_PinT_v1 = MAPS_SPI1_PinT_v1[['Genome', 'Type', 'Start', 'End', 'Strand', 'Locus_tag',
       'Name', 'MS2-PinT_SPI-1_vs_PinT_SPI-1_log2FoldChange', 'MS2-PinT_SPI-1_vs_PinT_SPI-1_pvalue',
       'MS2-PinT_SPI-1_vs_PinT_SPI-1_padj', 'MS2-PinT_SPI-1_vs_MS2_SPI-1_log2FoldChange', 'MS2-PinT_SPI-1_vs_MS2_SPI-1_pvalue',
       'MS2-PinT_SPI-1_vs_MS2_SPI-1_padj']]
       
	#We only focus on CDS rows
	MAPS_SPI1_PinT_v2 = MAPS_SPI1_PinT_v1.copy()
	MAPS_SPI1_PinT_v2=MAPS_SPI1_PinT_v1[(MAPS_SPI1_PinT_v2['Type']!='3UTR')  &  (MAPS_SPI1_PinT_v1['Type']!='5UTR')]
	MAPS_SPI1_PinT_v2=MAPS_SPI1_PinT_v2.rename(columns={'Name':'Gene_name',})
	return MAPS_SPI1_PinT_v2


#Pulse
def load_Pulse(pulse):
	print ("load Pulse")
	Pulse_SPI1_PinT = pd.read_csv(pulse)
	# Check how many unique Gene name are contained in this table
	print (f"Number of unique genes: {Pulse_SPI1_PinT['Gene name'].nunique()}")

	# select interesting columns and rename some of them
	Pulse_SPI1_PinT_v1=Pulse_SPI1_PinT.copy()
	Pulse_SPI1_PinT_v1=Pulse_SPI1_PinT_v1[['logFC','PValue','Gene name']]
	Pulse_SPI1_PinT_v1=Pulse_SPI1_PinT_v1.rename(columns={'Gene name':'Gene_name','logFC':'pulse_logFC','PValue':'pulse_PValue',})

	#print (Pulse_SPI1_PinT_v1.head())
	return Pulse_SPI1_PinT_v1


#Merge the two tables in one, keep only common genes 
#add sRNA sequences regardless of window sizes
def merge_tables(Pulse_SPI1_PinT_v1, MAPS_SPI1_PinT_v2):
	SPI1_MAPS_Pulse_commongenes = pd.merge(Pulse_SPI1_PinT_v1, MAPS_SPI1_PinT_v2, 
                       how='inner',
                       on=['Gene_name'])
	#print(SPI1_MAPS_Pulse_commongenes.tail())

	#Number of rows/genes
	print (f'Number of common genes between MAPS and Pulse: {len(SPI1_MAPS_Pulse_commongenes)}')

	#Save the data
	#SPI1_MAPS_Pulse_commongenes.to_csv('../new_saves/SPI1_MAPS_Pulse_commongenes.csv', index=False)


	
	#capturing the sequences
	SPI1_MAPS_Pulse_commongenes_temp = SPI1_MAPS_Pulse_commongenes[['Gene_name','Genome', 'Start', 'End','Strand']]
	SPI1_MAPS_Pulse_commongenes_loc = SPI1_MAPS_Pulse_commongenes_temp.copy()

	# We want to get the 50nt before and 50nt after the start codon, keep in mind for genes on the - strand the 
	#  "End" location coordinate is the first nucleotide. 

	SPI1_MAPS_Pulse_commongenes_loc['from']=np.where(SPI1_MAPS_Pulse_commongenes_loc['Strand'] == '+', SPI1_MAPS_Pulse_commongenes_loc['Start'] - 51, SPI1_MAPS_Pulse_commongenes_loc['End'] - 51)

	SPI1_MAPS_Pulse_commongenes_loc['to']=np.where(SPI1_MAPS_Pulse_commongenes_loc['Strand'] == '+', SPI1_MAPS_Pulse_commongenes_loc['Start'] + 51, SPI1_MAPS_Pulse_commongenes_loc['End'] + 51)
   
	#print(SPI1_MAPS_Pulse_commongenes_loc.tail(10))

	# Save the location of the genes
	#SPI1_MAPS_Pulse_commongenes_loc.to_csv('new_saves/SPI1_MAPS_Pulse_commongenes_loc.csv', index=False)

#read and summarize IntaRNA output

def summarise():
	PinT_SPI1_MFE=called()
	print (len(PinT_SPI1_MFE))
	return PinT_SPI1_MFE

#shuffling

def shuffle():
	#shuffle_nts()
	pass #comment as necessary

#calling IntaRNA after shuffling
#figure out if IntaRNA can be packaged together
#probs with the ViennaRNA tarball
#w/ some instructions

def call_intarna():
	pass #comment as necessary
	print ("Calling intaRNA")
	#call("sh call_intarna.sh", shell=True)
	print("intaRNA run finished")
	

#start with p-value calculation

def call_pval(PinT_SPI1_MFE):
	calc_pval(PinT_SPI1_MFE)
	#return?

def exec(maps, pulse, rilseq, intarna):
	#file locations should come from landing script
	Pulse_SPI1_PinT_v1 = load_Pulse(pulse)
	MAPS_SPI1_PinT_v2 = load_MAPS(maps)
	merge_tables(Pulse_SPI1_PinT_v1, MAPS_SPI1_PinT_v2)
	PinT_SPI1_MFE = summarise()
	if intarna == True:
		call_intarna()
	call_pval(PinT_SPI1_MFE)


