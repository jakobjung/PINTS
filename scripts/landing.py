from argparse import ArgumentParser
import sys
from pathlib import Path
from data_prep_initial import exec
from combine import start
import os

def main(maps, pulse, rilseq, intarna):
	print ("Starting data preprocessing...")
	exec(maps, pulse, rilseq, intarna)#should contain files
	print ("End of data preprocessing")
	print ("p value combination...")
	start()
	
if __name__ == "__main__":
	parser = ArgumentParser(
	prog="catchy name",#name of program to be displayed in the help msg
	description="sRNA target hybridization. Please make sure that the MAPS and the PULSE datasets are in the following format:",#description
	epilog="",)#final msg to be displayed, if anything
	
	params = parser.add_argument_group("parameters")
	params.add_argument(
	    "-d", "--path",
	    default="../data_store/",
	    #metavar=("X", "Y"),#to add as help/usage what possible values can be selected
	    help="path to the data files (default: %(default)s)",
	)
	#following params accept 1 or more files
	params.add_argument("-m", "--maps", nargs="+", 
	    help="MAPS input data, accepts multiple datasets")
	params.add_argument("-p", "--pulse", nargs="+",
	    help="Pulse input data, accepts multiple datasets")
	params.add_argument("-ril", "--rilseq", nargs="+",
	    help="RILSeq input data, accepts multiple datasets")
	params.add_argument("-x",
	    help="Future: other datasets")
	params.add_argument("-i", "--intarna", action="store_true",
	    help="hybridization using intaRNA")

	params.add_argument(
	    "--version", action="version", version="%(prog)s 0.1")
    
	args = parser.parse_args()
	target_dir = Path(args.path)
	
	if not target_dir.exists():
		parser.exit(1, message="The target directory doesn't exist\n")
	#print(args)
	print(args.maps[0])#checking if storage works properly
	print(args.intarna)
	print(args.path)
	os.chdir(args.path)
	print (os.getcwd())
	#pass on the arguments, ie, filenames obtained
	    
	main(args.maps[0], args.pulse[0], args.rilseq[0], args.intarna)

	##if nothing is input, return error

