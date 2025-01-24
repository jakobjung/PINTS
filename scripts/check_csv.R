# Script that reads in file from Rscript call, and renames all the columns
# as the filename_ + column name to avoid conflicts with other files.
# it also checks whether the columns "p_value" and "gene_name" are present. If not, it will
# cause an error as output of the script.

# load packages
library(tidyverse)


# import the csv file name from the bash Rscript call
args <- commandArgs(trailingOnly = TRUE)
csv_file <- args[1]
# csv_file <- "./data/MAPS_test.csv"

# extract the filename from the csv file without .csv
name_output <- gsub(".csv", "", basename(csv_file))

# read in the csv file
csv_file <- read_csv(csv_file)

# check if the columns p_value and gene_name are present
if (!("p_value" %in% colnames(csv_file)) | !("gene_name" %in% colnames(csv_file))) {
  stop("Columns 'p_value' and 'gene_name' are not present in the file.")
} else {
  print("Columns 'p_value' and 'gene_name' are present in the file.")
}


# rename all the columns to avoid conflicts with other files, except for gene_name
for (col in colnames(csv_file)) {
  if (col != "gene_name") {
    colnames(csv_file)[colnames(csv_file) == col] <- paste0(name_output, "_", col)
  }
}


# rename all "-" and " " to "_" in the column names
colnames(csv_file) <- gsub("-", "_", colnames(csv_file))
colnames(csv_file) <- gsub(" ", "_", colnames(csv_file))

# write the tibble to a file
write_tsv(csv_file, paste0(name_output, ".tsv"))

# print the head of the tibble
head(csv_file)
