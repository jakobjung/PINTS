#the shuffle
#generate 256 shuffled sequences 32 times
#prerequisite easel: https://github.com/EddyRivasLab/easel
#install the whole package or just the esl-shuffle bin(or C script) will do?

from subprocess import call

def shuffle_nts():
	for x in range(1,33):
	    start = datetime.now()
	    call("/home/rituparno/Stockpile/tools/easel/miniapps/esl-shuffle -N 256 -d -o ../new_saves/SPI1_PinT_shuffled_256_"+str(x)+".fa ../data_files/SPI1_PinT_extracted_seqs.fa", shell=True)
	    print (datetime.now() - start)
	    
