import sys
import config
import numpy as np
import math
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


	num = 0
	seqdict = {}
	idx = 0
	for line in lines:
		if ">" in line:
			data = line.strip().split()
			# num = int(data[0][12:].lstrip("0"))
			num = int(data[0][8:].lstrip("0"))
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
	config.__INITPROB__ = {'I': 1, 'A': 0, 'G': 0, 'Z': 0} # given in instructions (A = start, Z = stop)
	config.__TRANSPROB__ = {('I', 'I'): ((config.__AVGINTLEN__ - 1)/ config.__AVGINTLEN__), ('I', 'A'): (1/config.__AVGINTLEN__), ('I', 'G'): 0, ('I', 'Z'): 0,
							('A', 'I'): 0, ('A', 'A'): 0, ('A', 'G'): 1, ('A', 'Z'): 0,
							('G', 'I'): 0, ('G', 'A'): 0, ('G', 'G'): (((config.__AVGGENELEN__/3) - 1)/ (config.__AVGGENELEN__/3)), ('G', 'Z'): (1/ (config.__AVGGENELEN__/3)),
							('Z', 'I'): 1, ('Z', 'A'): 0, ('Z', 'G'): 0, ('Z', 'Z'): 0}
	#EMISSION probabilities
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

def viterbi(iter):
	seq = config.__SEQ__[iter]
	# seq= 'GCGATGCGTCTCATTTATAAAATATGGAATTATTATAGATTGATTGCATAAGCTATTCTCCAGCTTTATTTGGCTAGAGTAACTTTATGAAAACTAAAATCATTTTATGTTCAGCGGTACTAGCAGTACTTTCTGGTTGTGCGTCAGTACCTATGGTTGATTCTGAACTCTCCGATCAAGCGAAACAGTTTGATGCGCCAACCGAAGGGAAAGCGGGCGTGTATGTATATCGTCCAGAATCTGGCATTGGTGGTGCACTGAAAAAAGATGTGCATATTGATGGTGAATGCATTGGTGAAACGGCACCGGGTGTTTTCTTCTACCACGAAGTGGATGGCGATAAAGAGCACATTGTCAGTACCGAATCTGAATTTTCTCCAAATGAAGTCACCTTGTTTACTGAGCAAGGACGCCTCTATTTTGTTCAGCAATACATCAAAATGGGCGCATTTGTTGGCGGTGCGGATTTAGTGGTTGTCGATGAGTCAACGGGTAAATCTGACGTTTACAAAACAAAAATGGCGATCAAGGGCAATTGCTCTGCCAAGTAATTCGACGTTCTCTCGTCTCCCAGACCCAAGCTTTGTGCTTGGGTTTTTA'
	# seq = 'ACATGACGGGATAGGTTGAAATAGGGACGTGACATAGTTAA'
	# seq = 'ATGACAGAT'
	# seq= 'GCGATGCGTCTCATTTATAAAATATGGAATTATTATAGATTGATTGCATAAGCTATTCTCCAGCTTTATTTGGCTAGAGTAACTTTATGAAAACTAAAATCATTTTATGTTCAGCGGTACTAGCAGTACTTTCTGGTTGTGCGTCAGTACCTATGGTTGATTCTGAACTCTCCGATCAAGCGAAACAGTTTGATGCGCCAACCGAAGGGAAAGCGGGCGTGTATGTATATCGTCCAGAATCTGGCATTGGTGGTGCACTGAAAAAAGATGTGCATATTGATGGTGAATGCATTGGTGAAACGGCACCGGGTGTTTTCTTCTACCACGAAGTGGATGGCGATAAAGAGCACATTGTCAGTACCGAATCTGAATTTTCTCCAAATGAAGTCACCTTGTTTACTGAGCAAGGACGCCTCTATTTTGTTCAGCAATACATCAAAATGGGCGCATTTGTTGGCGGTGCGGATTTAGTGGTTGTCGATGAGTCAACGGGTAAATCTGACGTTTACA'
	# seq = 'GCGATGCGTCTCATTTATAAAATATGGAATTATTGA'
	states = ['I', 'A', 'G', 'Z']

	codonflag = 1
	start = False
	inf = math.inf
	path = [{'I':'I','A':'I','G':'I','Z':'I'}]
	currprob = {'I': config.__INTEREMIT__[seq[0]], 'A': -inf, 'G': -inf, 'Z': -inf}
	codon = ''
	for i in range(0, len(seq)-2):
		lastprob = currprob.copy()
		currpath = {}
		if start:
			codonflag +=1
			if codonflag > 3: 
				codonflag = 1
		if ((not start) or codonflag==1) and (len(seq[i:i+3]) == 3):
			codon = seq[i:i+3]
			codonflag = 1
			if config.__STARTEMIT__[codon]>0:
				start = True

		for k in range(0, len(states)):
			tmpprob = {}
			curr = states[k]
			if curr == 'I':
				emitprob = config.__INTEREMIT__
			elif curr == 'G':
				emitprob = config.__GENEEMIT__
			elif curr == 'A':
				emitprob = config.__STARTEMIT__
			elif curr == 'Z':
				emitprob = config.__STOPEMIT__
			if curr == 'I':
				if start:
					currprob[curr] = -inf
				else:
					for j in range(0, len(states)):
						prev = states[j]
						if config.__TRANSPROB__[(prev, curr)]>0 and emitprob[seq[i]]>0:
							tmpprob[prev] = lastprob[prev] + log(emitprob[seq[i]]) + log(config.__TRANSPROB__[(prev, curr)])
						else:
							tmpprob[prev] = -inf
					currprob[curr] = max(tmpprob.values())
					currpath[curr] = max(tmpprob, key = tmpprob.get)
			else:
				if start and len(seq[i:i+3]) == 3:
					if codonflag == 1:
						for j in range(0, len(states)):
							prev = states[j]
							if config.__TRANSPROB__[(prev, curr)]>0 and emitprob[codon]>0:
								tmpprob[prev] = lastprob[prev] + log(emitprob[seq[i:i+3]]) + log(config.__TRANSPROB__[(prev, curr)])
							else:
								tmpprob[prev] = -inf
						currprob[curr] = max(tmpprob.values())
						currpath[curr] = max(tmpprob, key = tmpprob.get)
					elif codonflag==2 or codonflag == 3:
						currpath[curr] = curr
				else:
					currprob[curr] = -inf
		if config.__STOPEMIT__[codon]>0 and codonflag==3:
			take = False
		path.append(currpath)
		endincodon = False
		if codonflag == 3:
			endincodon=True
	
	traceback(path, currprob)

def traceback(path_trellis, currprob):
	state_seq = ''
	curr_state = max(current_prob, key = current_prob.get)
	for i in range(len(path_trellis)):
		index = len(path_trellis) - i - 1
		state_seq = path_trellis[index][curr_state] + state_seq
		curr_state = path_trellis[index][curr_state]

	start = False
	start_index = 0
	end_index = 0
	result_file.write('###\n')
	for i in range(len(state_seq)):
		if start==False and state_seq[i]=='S':
			start_index = i
			start = True
		elif start == True and state_seq[i]=='T':
			end_index = i + 2
			start = False
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(start_index, end_index))
	if end_index==0:
		result_file.write(seq_name+'\t ena\t CDS \t .\t .\t . \t +\t 0\t .\n')

# # print(path)
# 	states = ''
# 	currstate = max(currprob, key = currprob.get)
# 	states += currstate
	
# 	for i in range(len(path)):
# 		idx = len(path) - i - 1
# 		states = path[idx][currstate] + states
# 		currstate = path[idx][currstate]

# 	print(states)

# 	i = 0
# 	start = 0
# 	pts = [] # DONT subtract 1 because idx 0 is empty
# 	while i < len(path):
# 		if path[i] != 'I':
# 			start = i
# 			j = i
# 			while j < len(path) and path[j] != 'I':
# 				j+=1
# 			pts.append((start, j))
# 			i = j
# 		else:
# 			i+=1

	# with open('1c.gff3', 'a') as f:
	# 	print('###', file = f)
	# 	for i in pts:
	# 		print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s" %(config.__NAMES__[iter-1], 'ena', 'CDS', i[0], i[1], '.', '+', '0'), file = f)
	# f.close()

def main(argv):
	seq = argv[0]
	conf = argv[1]
	getconfigs(conf)
	getseq(seq)
	getProbTables()

	for i in range(1, len(config.__SEQ__)+1):
		viterbi(i)
		print('end iter ', i)
		break

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa', 'configuration.txt'])
