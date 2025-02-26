import sys
import config
import numpy as np
from math import log, inf

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
	idx = -1
	for line in lines:
		if ">" in line:
			data = line.strip().split()
			# num = int(data[0][12:].lstrip("0"))
			num = int(data[0][8:].lstrip("0"))
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
	# config.__STARTCODON__ = ['ATG']

	config.__INITPROB__ = {'I': 1, 'A': 0, 'G': 0, 'Z': 0} # given in instructions (A = start, Z = stop)
	config.__TRANSPROB__ = {('I', 'I'): 1-avgintprob, ('I', 'A'): avgintprob, ('I', 'G'): 0, ('I', 'Z'): 0,
							('A', 'I'): 0, ('A', 'A'): 0, ('A', 'G'): 1, ('A', 'Z'): 0,
							('G', 'I'): 0, ('G', 'A'): 0, ('G', 'G'): 1-avggenprob, ('G', 'Z'): avggenprob,
							('Z', 'I'): 1, ('Z', 'A'): 0, ('Z', 'G'): 0, ('Z', 'Z'): 0}#EMISSION probabilities
	config.__INTEREMIT__ = {k: v/ sum(config.__NTFREQ__.values()) for k,v in config.__NTFREQ__.items()} 
	
	#find frequencies of start codons
	config.__STARTEMIT__ = {k: 0 for k in config.__CODONFREQ__.keys()}
	startfreq = sum(config.__CODONFREQ__[i] for i in config.__STARTCODON__)
	for i in config.__STARTCODON__:
		config.__STARTEMIT__[i] = config.__CODONFREQ__[i]/startfreq 

	# find frequencies of stop codons
	config.__STOPEMIT__ = {k: 0 for k in config.__CODONFREQ__.keys()}
	stopfreq = sum(config.__CODONFREQ__[i] for i in config.__STOPCODONS__)
	for i in config.__STOPCODONS__:
		config.__STOPEMIT__[i] = config.__CODONFREQ__[i]/stopfreq 

	config.__CODONFREQ__['TAA'] = 0
	config.__CODONFREQ__['TAG'] = 0
	config.__CODONFREQ__['TGA'] = 0
	config.__GENEEMIT__ = {k: v/ sum(config.__CODONFREQ__.values()) for k,v in config.__CODONFREQ__.items()}
	
	return

def viterbi(it):
	if not (config.__SEQ__.get(it)):
		return
	seq = config.__SEQ__[it]
	# seq= 'GCGATGCGTCTCATTTATAAAATATGGAATTATTATAGATTGATTGCATAAGCTATTCTCCAGCTTTATTTGGCTAGAGTAACTTTATGAAAACTAAAATCATTTTATGTTCAGCGGTACTAGCAGTACTTTCTGGTTGTGCGTCAGTACCTATGGTTGATTCTGAACTCTCCGATCAAGCGAAACAGTTTGATGCGCCAACCGAAGGGAAAGCGGGCGTGTATGTATATCGTCCAGAATCTGGCATTGGTGGTGCACTGAAAAAAGATGTGCATATTGATGGTGAATGCATTGGTGAAACGGCACCGGGTGTTTTCTTCTACCACGAAGTGGATGGCGATAAAGAGCACATTGTCAGTACCGAATCTGAATTTTCTCCAAATGAAGTCACCTTGTTTACTGAGCAAGGACGCCTCTATTTTGTTCAGCAATACATCAAAATGGGCGCATTTGTTGGCGGTGCGGATTTAGTGGTTGTCGATGAGTCAACGGGTAAATCTGACGTTTACAAAA'
	# seq= 'GCGATGCGTCTCATTTATAAAATATGGAATTATTATAGATTGATTGCATAAGCTATTCTCCAGCTTTATTTGGCTAGAGTAACTTTATGAAAACTAAAATCATTTTATGTTCAGCGGTACTAGCAGTACTTTCTGGTTGTGCGTCAGTACCTATGGTTGATTCTGAACTCTCCGATCAAGCGAAACAGTTTGATGCGCCAACCGAAGGGAAAGCGGGCGTGTATGTATATCGTCCAGAATCTGGCATTGGTGGTGCACTGAAAAAAGATGTGCATATTGATGGTGAATGCATTGGTGAAACGGCACCGGGTGTTTTCTTCTACCACGAAGTGGATGGCGATAAAGAGCACATTGTCAGTACCGAATCTGAATTTTCTCCAAATGAAGTCACCTTGTTTACTGAGCAAGGACGCCTCTATTTTGTTCAGCAATACATCAAAATGGGCGCATTTGTTGGCGGTGCGGATTTAGTGGTTGTCGATGAGTCAACGGGTAAATCTGACGTTTACA'
	# seq = 'ACATGACGGGATAGGTTGAAATAGGGACGTGACATAGTTAA'
	# seq = 'ATGAAAGATA'
	# seq = 'GCTGAGGTGACGTGCAACAGTCGATTCGTGAATGCGAAGAGCTTGAGAAATCATCGCCTGACTCCAACCTTCAGACGCAAGTAATACCGCTTTGATGCGGTCACGCACTCGACCATCACGAGTGGAATCGTGCATCTCTTCGAGTTGTAGTTTCTGTTGGGAAGTCAGTATTATTTTCATGGTGAGTAGAATGATCCTGATTCCATGAAAAATCAAGCATCTTCAATGATCACGGGTATATGTGCTCGCGTTTTTGACCTTCGGTTTCAGTGTTCTCGGCACTTTTATTGTCCGCTCGGGGATTTTGACATCGGTCCATGCGTTTGCCGTGGATCCAACCAAAGGTATTGTGCTTTTGCTGGTCATGGCGTTCATTTTTTTACTCACTTTTGCGTTATTGATCCTCAAAAGCGATAGCATTCCCGCTAAAGCCATTACCCATTGGCTAAGTCGCCAATACCTTACGGTGGTGGCGATGGGACTGTTACTGATCGCAACCAGTACCGTGTTCCTTGGCACCTTCTACCCAATGATTTATGAAAAATGGAATAG'

	states = ['I', 'A', 'G', 'Z']

	codon = False
	stop = False
	end = False
	start = False
	firststart = False
	i = 0
	prev = {'I':'I', 'A':'I', 'G':'I', 'Z':'I'}
	path = [{'I':'I','A':'I','G':'I','Z':'I'}]
	currprob = {'I': config.__INTEREMIT__[seq[0]], 'A': 0, 'G': 0, 'Z': 0}
	
	while i < len(seq)-2:
		prevprob = currprob.copy()
		codonpath1 = {}
		codonpath2 = {}
		currpath = {}
		if seq[i:i+3] != '':
			for k in range(0, len(states)):
				curr = states[k]
				if seq[i:i+3] in config.__STARTCODON__:
					codon = True
					if not (firststart):
						start = True
						firststart = True
				if seq[i:i+3] in config.__STOPCODONS__ and (codon):
					stop = True
				for j in range(0, len(states)):
					prevstate = states[j]
					if curr == 'I':
						emitprob = config.__INTEREMIT__
					elif curr == 'G':
						emitprob = config.__GENEEMIT__
					elif curr == 'A':
						emitprob = config.__STARTEMIT__
					elif curr == 'Z':
						emitprob = config.__STOPEMIT__
					if codon:
						if curr == 'I':
							currprob[prevstate] = 0
						else:
							if stop:
								if curr == 'G':
									currprob[prevstate] = 0
								else:
									if config.__TRANSPROB__[(prevstate, curr)]>0 and emitprob[seq[i:i+3]]>0:
										currprob[prevstate] = prevprob[prevstate] + log(emitprob[seq[i:i+3]]) + log(config.__TRANSPROB__[(prevstate, curr)])
									else:
										currprob[prevstate] = 0
							else:
								if (len(seq[i:i+3]) == 3):
									if config.__TRANSPROB__[(prevstate, curr)]>0 and emitprob[seq[i:i+3]]>0:
										currprob[prevstate] = prevprob[prevstate] + log(emitprob[seq[i:i+3]]) + log(config.__TRANSPROB__[(prevstate, curr)])
									else:
										currprob[prevstate] = 0
								else:
									end = True
					else:
						if curr == 'I':
							if config.__TRANSPROB__[(prevstate, curr)]>0 and emitprob[seq[i]]>0:
								currprob[prevstate] = prevprob[prevstate] + log(emitprob[seq[i]]) + log(config.__TRANSPROB__[(prevstate, curr)])
							else:
								currprob[prevstate] = 0
						else:
							if (len(seq[i:i+3]) == 3):
								if config.__TRANSPROB__[(prevstate, curr)]>0 and emitprob[seq[i:i+3]]>0:
									currprob[prevstate] = prevprob[prevstate] + log(emitprob[seq[i:i+3]]) + log(config.__TRANSPROB__[(prevstate, curr)])
								else:
									currprob[prevstate] = 0
							else:
								end = True
				
				absval = {key : abs(val) for key, val in currprob.items()} 
				maxst = max(absval, key=absval.get)

				if max(absval.values()) == 0:
					maxst = ''
				if not maxst == '':
					currpath[curr] = maxst
				if codon:
					if maxst != '':
						prev[curr] = curr
						currpath[curr] = maxst
						codonpath1[curr] = curr
						codonpath2[curr] = curr
						
				else:
					if maxst != '':
						prev[curr] = maxst
						currpath[curr] = maxst
				
		path.append(currpath)
		if codon:
			path.append(codonpath1)
			path.append(codonpath2)
			if stop:
				stop = False
				codon = False
				firststart = False
			if start:
				start = False
			i+=3
		else:
			i+=1
	assignment = traceback(currprob, iter, path)

def traceback(currprob, iter, path):
	states = ''
	currstate = max(currprob, key = currprob.get)
	states += currstate
	
	for i in range(len(path)):
		idx = len(path) - i - 1
		states = path[idx][currstate] + states
		currstate = path[idx][currstate]
	states = states[1:]
	print(states)

	i = 0
	start = 0
	pts = [] # DONT subtract 1 because idx 0 is empty
	while i < len(path):
		if path[i] != 'I':
			start = i
			j = i
			while j < len(path) and path[j] != 'I':
				j+=1
			pts.append((start, j))
			i = j
		else:
			i+=1


	with open('test.gff3', 'a') as f:
		print('###', file = f)
		for i in pts:
			print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s" %(config.__NAMES__[iter-1], 'ena', 'CDS', i[0], i[1], '.', '+', '0'), file = f)
	f.close()

def main(argv):
	seq = argv[0]
	conf = argv[1]
	getconfigs(conf)
	getseq(seq)
	getProbTables()
	for i in range(len(config.__SEQ__)):
		viterbi(i)
		print('end iter ', i)
		break

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa', 'configuration.txt'])
