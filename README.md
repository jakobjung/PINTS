# PINTS

- Project name: PINTS - **P**-value **IN**tegration for **T**argets of **s**RNAs

- Experiments: Hoda Kooshapour

- Supervision: Alexander Westermann

- Data analysis/algorithm development: Hoda Kooshapour, Jakob Jung

- Date: 07-2024

  

## Introduction
PINTS is a computational pipeline for the identification of sRNA targets. It integrates IntaRNA target predictions 
with experimental data, such as MAPS (MS2-affinity purification coupled with RNA sequencing) and other Sequencing data.
 

## Directory structure

The project is divided in 3 main directories:

-   [data](data) : contains all raw, intermediate and final data of the project.  
-   [analysis](analysis): contains analysis files, such as figures, plots, tables, etc. 
-   [scripts](scripts): contains scripts to process and analyze data from the data directory.

Some directories have their own README.md file with information on the respective files. 



## Workflow

Here I describe the workflow, which can be followed to fully understand the procedure.



### 1. Prerequisites

For running the whole analysis, one needs following packages/tools/software:

- IntaRNA (v3.4.0) 

- python (v3.12.3) along with packages from biopython, numpy, pandas, matplotlib, seaborn, scipy, statsmodels, etc.

- Linux shell (we used Ubuntu 20.04) for commands & bash scripts

- bedtools (v2.31.1)

- easel (v0.49)

  



### 2. ...

