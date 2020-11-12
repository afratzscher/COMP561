import sys
import config
import numpy as np
from math import inf
from math import log

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

	seqdict = {}
	idx = -1
	for line in lines:
		if ">" in line:
			data = line.strip().split()
			config.__NAMES__.append(data[0][1:])
			idx+=1
		else:
			data = line.strip('\n')
			if idx in seqdict:
				seqdict[idx] = seqdict[idx] + data
			else:
				seqdict[idx] = data
	config.__SEQ__ = seqdict
	return

def getProbTables():
	avgintprob = 1/ config.__AVGINTLEN__
	avggenprob = 3/ config.__AVGGENELEN__

	config.__INITPROB__ = {'I': 1, 'A': 0, 'G': 0, 'Z': 0} # given in instructions (A = start, Z = stop)
	config.__TRANSPROB__ = {('I', 'I'): 1-avgintprob, ('I', 'A'): avgintprob, ('I', 'G'): 0, ('I', 'Z'): 0,
							('A', 'I'): 0, ('A', 'A'): 0, ('A', 'G'): 1, ('A', 'Z'): 0,
							('G', 'I'): 0, ('G', 'A'): 0, ('G', 'G'): 1-avggenprob, ('G', 'Z'): avggenprob,
							('Z', 'I'): 1, ('Z', 'A'): 0, ('Z', 'G'): 0, ('Z', 'Z'): 0}
	#EMISSION probabilities
	interemit = {k: v/ sum(config.__NTFREQ__.values()) for k,v in config.__NTFREQ__.items()} 

	#find frequencies of start codons
	startemit = {k: 0 for k in config.__CODONFREQ__.keys()}
	startfreq = sum(config.__CODONFREQ__[i] for i in config.__STARTCODON__)
	for i in config.__STARTCODON__:
		startemit[i] = config.__CODONFREQ__[i]/startfreq 
	
	# find frequencies of stop codons
	stopemit = {k: 0 for k in config.__CODONFREQ__.keys()}
	stopfreq = sum(config.__CODONFREQ__[i] for i in config.__STOPCODONS__)
	for i in config.__STOPCODONS__:
		stopemit[i] = config.__CODONFREQ__[i]/stopfreq
		# stopemit[i] = 1/3

	config.__CODONFREQ__['TAA'] = 0
	config.__CODONFREQ__['TAG'] = 0
	config.__CODONFREQ__['TGA'] = 0
	geneemit = {k: v/ sum(config.__CODONFREQ__.values()) for k,v in config.__CODONFREQ__.items()}

	config.__EMITPROB__ = {'I': interemit, 'A': startemit, 'G': geneemit, 'Z': stopemit}
	return

def viterbi(it):
	if not (config.__SEQ__.get(it)):
		return
	seq = config.__SEQ__[it]
	initprob = config.__INITPROB__
	emitprob = config.__EMITPROB__
	transprob = config.__TRANSPROB__
	take = False
	codflag = 1
	currprob = {'I': emitprob['I'][seq[0]], 'A': -inf, 'G': -inf, 'Z': -inf} 
	path = [{'I': 'I','A': 'I','G': 'I','Z': 'I'}]
	states = ['I', 'A', 'G', 'Z']

	for i in range(1, len(seq)-2):
		last_prob = currprob.copy()
		currpath = {}

		if take:
			codflag += 1
			if codflag > 3:
				codflag = 1
		if (not take) or codflag==1:
			codon = seq[i:i+3]
			codflag = 1
			if emitprob['A'][codon]>0:
				take = True

		for curr in states:
			prob_list = {}
			if curr == 'I':
				for prev in states:
					if transprob[(prev,'I')]>0 and emitprob['I'][seq[i]]>0:
						prob_list[prev] = last_prob[prev] + log(transprob[(prev,'I')]) + log(emitprob['I'][seq[i]])
					else:
						prob_list[prev] = -inf
				currprob['I'] = max(prob_list.values())
				currpath['I'] = max(prob_list, key = prob_list.get)
		
		
			if curr=='A' or curr=='G' or curr=='Z':
				if take == True:
					if codflag==1:
						for prev in initprob:
							if transprob[(prev,curr)]>0 and emitprob[curr][codon]>0:
								prob_list[prev] = last_prob[prev] + log(transprob[(prev,curr)]) + log(emitprob[curr][codon])
							else:
								prob_list[prev] = -inf

						currprob[curr] = max(prob_list.values())
						currpath[curr] = max(prob_list, key = prob_list.get)
						
					elif codflag==2 or codflag==3:
						currpath[curr] = curr
				else:
					currprob[curr] = -inf
		if emitprob['Z'][codon]>0 and codflag==3:
			take = False
		path.append(currpath)

	traceback(path, currprob, it)

def traceback(path, currprob, it):
	states = ''
	curr = max(currprob, key = currprob.get)
	for i in range(0, len(path)):
		idx = len(path) - i - 1
		states = path[idx][curr] + states
		curr = path[idx][curr]
	states = states[1:]
	i = 0
	start = 0
	pts = [] # DONT subtract 1 because idx 0 is empty
	
	while i < len(states)-1:
		if states[i] != 'I':
			start = i
			while (states[i] != 'I') and (i < len(states)-1):
				i+=1
			pts.append((start, i-1))
		i+=1
	if states[len(states)-1] != 'I':
		point = pts[-1]
		pts.pop()
		pts.append((point[0], point[1]+1))

	with open('1c.gff3', 'a') as f:
		print('###', file = f)
		for i in pts:
			print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(config.__NAMES__[it], 'ena', 'CDS', i[0]+1, i[1]+1, '.', '+', '0', '.'), file = f)
	f.close()

def main(argv):
	seq = argv[0]
	conf = argv[1]
	getconfigs(conf)
	getseq(seq)
	getProbTables()
	for i in range(len(config.__SEQ__)):
		viterbi(i)

if __name__ == '__main__':
	main(sys.argv[1:])
	# main(['Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa', 'configuration.txt'])
