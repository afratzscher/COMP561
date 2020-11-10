import sys
import itertools
import config

def getconfigs(conf):
	f = open(conf, 'r')
	line = f.readlines()
	config.__AVGINTLEN__ = float(line[0].strip())
	config.__AVGGENELEN__ = float(line[1].strip())
	config.__NTFREQ__ = eval(line[2].strip())
	config.__CODONFREQ__ = eval(line[3].strip())
	config.__STARTCODON__ = eval(line[4].strip())
	f.close()
	
	return


def getseq(seq):
	seqfile = open(seq, "r")
	lines = seqfile.readlines()

	num = 0
	seqdict = {}
	idx = 0
	for line in lines:
		if ">" in line:
			data = line.strip().split()
			num = int(data[0][12:].lstrip("0"))
			config.__NAMES__.append(data[0][1:])
		else:
			data = line.strip('\n')
			if num in seqdict:
				seqdict[num] = seqdict[num] + data
			else:
				seqdict[num] = data
	config.__SEQ__ = seqdict
	return


def getProbTables():
	config.__initprob__ = {'I': 1, 'START': 0, 'G': 0, 'STOP': 0} # given in instructions
	config.__transprob__ = {('I', 'I'): ((config.__AVGINTLEN__ - 1)/ config.__AVGINTLEN__), ('I', 'START'): (1/config.__AVGINTLEN__), ('I', 'G'): 0, ('I', 'STOP'): 0,
							('START', 'I'): 0, ('START', 'START'): 0, ('START', 'G'): 1, ('START', 'STOP'): 0,
							('G', 'I'): 0, ('G', 'START'): 0, ('G', 'G'): (((config.__AVGGENELEN__/3) - 1)/ (config.__AVGGENELEN__/3)), ('G', 'STOP'): (1/ (config.__AVGGENELEN__/3)),
							('STOP', 'I'): 1, ('STOP', 'START'): 0, ('STOP', 'G'): 0, ('STOP', 'STOP'): 0}
	#EMISSION probabilities
	config.__interEmitProb__ = {k: v/ sum(config.__NTFREQ__.values()) for k,v in config.__NTFREQ__.items()} 
	
	#find frequencies of start codons
	config.__startEmitProb__ = {k: 0 for k in config.__CODONFREQ__.keys()}
	startfreq = sum(config.__CODONFREQ__[i] for i in config.__STARTCODON__)
	for i in config.__STARTCODON__:
		config.__startEmitProb__[i] = config.__CODONFREQ__[i]/startfreq 

	config.__geneEmitProb__ = {k: v/ sum(config.__CODONFREQ__.values()) for k,v in config.__CODONFREQ__.items()}
	
	# find frequencies of stop codons
	STOPCODONS = ['TAG', 'TTA', 'TGA']
	config.__stopEmitProb__ = {k: 0 for k in config.__CODONFREQ__.keys()}
	stopfreq = sum(config.__CODONFREQ__[i] for i in STOPCODONS)
	for i in STOPCODONS:
		config.__stopEmitProb__[i] = config.__CODONFREQ__[i]/stopfreq 

	return

def viterbi(iter):
	print(config.__SEQ__[iter]) 



def main(argv):
	seq = argv[0]
	conf = argv[1]
	getconfigs(conf)
	getseq(seq)
	getProbTables()
	
	for i in range(1, len(config.__SEQ__)+1):
		name = config.__NAMES__[i]
		viterbi(i)
		break

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['Vibrio_cholerae.GFC_11.dna.toplevel.fa', 'configuration.txt'])
