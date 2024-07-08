#!/bin/bash
set -e
#set -o pipefail
################################################################################
#                              PINTS command line                              #
#                                                                              #
# PINTs (P-value Integration for Targets of sRNAs) is a computational pipeline #
# of bacterial sRNAs.                                                          #
#                                                                              #
################################################################################
################################################################################
################################################################################
#                                                                              #
#  Created by Hoda Kooshapour, Jakob Jung                                      #
#                                                                              #
################################################################################
################################################################################
################################################################################

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "PINTS can be used to ..."
   echo
   echo "Version:     v1.0.0"
   echo "About:       Developed by Hoda Kooshapour & Jakob Jung @westermannlab @Barquistlab, HIRI"
   echo "Docs & Code: https://github.com/BarquistLab/mason_commandline"
   echo "Mail:        jakobjung@tutanota.com"
   echo
   echo "Usage:         sh PINTS.sh -f <fasta_file.fasta> -a <file1.csv> -b <file2.csv> [...]"
   echo
   echo "Options:"
   echo "Required:      "
   echo "               -f  FASTA file of target organism"
   echo "                   - (PATH)"
   echo "Help:"
   echo "                -h print this help menu"
   echo
   echo "Version:"
   echo "                -V print version"
   echo
}

################################################################################
################################################################################
# Main program                                                                 #
################################################################################
################################################################################
################################################################################
# Process the input options. Add options as needed.                            #
################################################################################
# I start with assigning the flags (user inputs):

version="v1.0.0"

while getopts f:a:b:c:d:s:hV flag
do
    case "${flag}" in
      f) fasta=${OPTARG};;
      a) file1=${OPTARG};;
      b) file2=${OPTARG};;
      c) file3=${OPTARG};;
      d) file4=${OPTARG};;
      s) srna=${OPTARG};;
      h) Help
         exit;;
      V) echo "$version"
         exit;;
      *) # Default case for handling invalid options
               echo "Invalid option: -$OPTARG" >&2
               Help
               exit 1
               ;;
    esac
done



# make fasta +f1 required
if [ -z "$fasta" ]
then
    echo "Please provide a fasta file"
    exit 1
fi


# run py_fastavalidator -f $fasta  and exit if the output is not empty (0)
if ! ( py_fasta_validator -f $fasta )
then
    echo "ERROR: Please provide a valid fasta file"
    exit 1
fi

# check for all csv files whether they have colnames of genome, gene_name, p_value, start, end, strand, type in their first row.
# it does not matter if the order is different
for file in $file1 $file2 $file3 $file4
do
    # run python script that checks for colnames etc and outputs an error message if the file is not valid
    python3 ./scripts/check_csv.py "$file"
done


# run following code if srna is provided
if [ -f "$srna" ]
then
  # run python script to prepare bed file from first csv file ($file1)
  python3 ./scripts/csv_to_bed.py "$file1"
  # run bedtools getfasta to get the sequences of the genes. the bed file is the input file but with .bed instead of
  # .csv
  echo "runnning bedtools getfasta for $file1"
  bedtools getfasta -fi "$fasta" -bed "${file1%.csv}".bed -fo "${file1%.csv}".fasta -s -name+
  # run intarna prediction for the sRNA and the fasta file
  echo "running IntaRNA for $file1"
  IntaRNA -q "$srna" -t "${file1%.csv}".fasta --out "${file1%.csv}"_intarna.csv --outMode C
  # now shuffle the input genome file 1 times and run IntaRNA again
  for i in 1
  do
    # shuffle the genome file
    echo "shuffling $fasta"
    esl-shuffle -d -o "${fasta%.fasta}"_shuffled"${i}".fasta "$fasta"
    # remove "-shuffled" from the fasta header
    sed -i 's/-shuffled//g' "${fasta%.fasta}"_shuffled"${i}".fasta
    # run bedtools getfasta to get the sequences of the genes. the bed file is the input file but with .bed instead of
    # .csv
    bedtools getfasta -fi "${fasta%.fasta}"_shuffled"${i}".fasta -bed "${file1%.csv}".bed -fo "${file1%.csv}"_shuffled"${i}".fasta -s -name+
    # run intarna prediction for the sRNA and the fasta file
    echo "shuffle $i"
    IntaRNA -q "$srna" -t "${file1%.csv}"_shuffled"${i}".fasta --out "${file1%.csv}"_shuffled"${i}"_intarna.csv --outMode C
    echo "shuffle done"
  done
  # now run get_pvalues_from_intarna.py to get the p-values from the intarna output. input is the intarna output file
  # and the csv file
  python ./scripts/get_pvalues_from_intarna.py "${file1%.csv}"_intarna.csv "${file1%.csv}"_shuffled"${i}"_intarna.csv

# else print error message
else
  echo "No sRNA fastafile provided!"
  exit 1
fi

# ok, now combine all the P-values from the csv files using the combine_pvalues.py script
python ./scripts/combine_pvalues.py $file1 $file2 $file3 $file4












