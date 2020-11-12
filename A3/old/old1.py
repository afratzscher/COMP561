import sys
import config
import numpy as np

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
	if not (config.__SEQ__.get(iter)):
		return
	seq = config.__SEQ__[iter]
	# seq= 'GCGATGCGTCTCATTTATAAAATATGGAATTATTATAGATTGATTGCATAAGCTATTCTCCAGCTTTATTTGGCTAGAGTAACTTTATGAAAACTAAAATCATTTTATGTTCAGCGGTACTAGCAGTACTTTCTGGTTGTGCGTCAGTACCTATGGTTGATTCTGAACTCTCCGATCAAGCGAAACAGTTTGATGCGCCAACCGAAGGGAAAGCGGGCGTGTATGTATATCGTCCAGAATCTGGCATTGGTGGTGCACTGAAAAAAGATGTGCATATTGATGGTGAATGCATTGGTGAAACGGCACCGGGTGTTTTCTTCTACCACGAAGTGGATGGCGATAAAGAGCACATTGTCAGTACCGAATCTGAATTTTCTCCAAATGAAGTCACCTTGTTTACTGAGCAAGGACGCCTCTATTTTGTTCAGCAATACATCAAAATGGGCGCATTTGTTGGCGGTGCGGATTTAGTGGTTGTCGATGAGTCAACGGGTAAATCTGACGTTTACAAAA'
	# seq= 'GCGATGCGTCTCATTTATAAAATATGGAATTATTATAGATTGATTGCATAAGCTATTCTCCAGCTTTATTTGGCTAGAGTAACTTTATGAAAACTAAAATCATTTTATGTTCAGCGGTACTAGCAGTACTTTCTGGTTGTGCGTCAGTACCTATGGTTGATTCTGAACTCTCCGATCAAGCGAAACAGTTTGATGCGCCAACCGAAGGGAAAGCGGGCGTGTATGTATATCGTCCAGAATCTGGCATTGGTGGTGCACTGAAAAAAGATGTGCATATTGATGGTGAATGCATTGGTGAAACGGCACCGGGTGTTTTCTTCTACCACGAAGTGGATGGCGATAAAGAGCACATTGTCAGTACCGAATCTGAATTTTCTCCAAATGAAGTCACCTTGTTTACTGAGCAAGGACGCCTCTATTTTGTTCAGCAATACATCAAAATGGGCGCATTTGTTGGCGGTGCGGATTTAGTGGTTGTCGATGAGTCAACGGGTAAATCTGACGTTTACA'
	# seq = 'ACATGACGGGATAGGTTGAAATAGGGACGTGACATAGTTAA'
	# seq = 'AATGGAATGA'
	# seq = 'GCTGAGGTGACGTGCAACAGTCGATTCGTGAATGCGAAGAGCTTGAGAAATCATCGCCTGACTCCAACCTTCAGACGCAAGTAATACCGCTTTGATGCGGTCACGCACTCGACCATCACGAGTGGAATCGTGCATCTCTTCGAGTTGTAGTTTCTGTTGGGAAGTCAGTATTATTTTCATGGTGAGTAGAATGATCCTGATTCCATGAAAAATCAAGCATCTTCAATGATCACGGGTATATGTGCTCGCGTTTTTGACCTTCGGTTTCAGTGTTCTCGGCACTTTTATTGTCCGCTCGGGGATTTTGACATCGGTCCATGCGTTTGCCGTGGATCCAACCAAAGGTATTGTGCTTTTGCTGGTCATGGCGTTCATTTTTTTACTCACTTTTGCGTTATTGATCCTCAAAAGCGATAGCATTCCCGCTAAAGCCATTACCCATTGGCTAAGTCGCCAATACCTTACGGTGGTGGCGATGGGACTGTTACTGATCGCAACCAGTACCGTGTTCCTTGGCACCTTCTACCCAATGATTTATGAAAAATGGAATAG'
	seq = ' ' + seq
	states = ['I', 'A', 'G', 'Z']
	tb = np.zeros((4, len(seq)))
	pointer = np.zeros((4, len(seq)), dtype=np.dtype('U4'))
	
	#initialize
	tb[0,0] = config.__INITPROB__['I'] 
	pointer[0,0] = 'E'

	codon = False
	stop = False
	end = False
	start = False
	firststart = False
	i = 1
	prev = 'I'
	while i < len(seq):
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
				probtb = {'I': 0, 'A': 0, 'G': 0, 'Z': 0}
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
							prob = 0
						else:
							if stop:
								if curr == 'G':
									prob = 0
								else:
									prob = tb[j,i-1] * emitprob[seq[i:i+3]] * config.__TRANSPROB__[(prevstate, curr)]
							else:
								if (len(seq[i:i+3]) == 3):
									prob = tb[j,i-1] * emitprob[seq[i:i+3]] * config.__TRANSPROB__[(prevstate, curr)]
								else:
									end = True
					else:
						if curr == 'I':
							prob = tb[j,i-1] * emitprob[seq[i]] * config.__TRANSPROB__[(prevstate, curr)]
						else:
							if (len(seq[i:i+3]) == 3):
									prob = tb[j,i-1] * emitprob[seq[i:i+3]] * config.__TRANSPROB__[(prevstate, curr)]
							else:
								end = True

					if prob > probtb[curr]:
						probtb[curr] = prob

				maxst = max(probtb, key=probtb.get)
				if max(probtb.values()) == 0 and codon:
					maxst = ''
				if codon:
					if maxst == '':
						pointer[k,i:i+3] = prev
					else:
						pointer[k, i] = prev
						pointer[k, i+1:i+3] = maxst
						prev = maxst
					tb[k, i:i+3] = probtb[curr]

				else:
					# if end:
					# 	pointer[2, i:i+3] = 'G'
					tb[k, i] = probtb[curr]
					if maxst == '':
						pointer[k,i] = prev
					else:
						pointer[k, i] = prev
						prev = maxst
		if codon:
			if stop:
				stop = False
				codon = False
				firststart = False
				tb[:,i] = 0
				tb[3,i] = 1
			if start:
				start = False
			i+=3
		else:
			i+=1
	print(tb[:, 530:540])
	assignment = traceback(tb, iter, pointer)

def traceback(tb, iter, pointer):
	states = ['I', 'A', 'G', 'Z']
	indices = np.argwhere(tb[:,-1] == np.amax(tb[:,-1])).flatten().tolist()
	if indices[0] == 0 and pointer[indices[0], -2] == 'G':
		indices = [2]
	check = False
	for idx in indices:
		if pointer[idx,-1] == '':
			idx +=1
		else:
			path = ''
			i = tb.shape[1] - 1
			while i >= 0:
				path += states[idx]
				curr = pointer[idx, i]
				if curr == 'I':
					idx = 0
				elif curr == 'A':
					idx = 1
				elif curr == 'G':
					idx = 2
				elif curr == 'Z':
					idx = 3
				i-=1
			path = (path[::-1])[1:]
			break
	i = 0
	pts = []
	print(path)

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

	with open('1c.gff3', 'a') as f:
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

	for i in range(1, len(config.__SEQ__)+1):
		viterbi(i)
		print('done ', i)
		break

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa', 'configuration.txt'])
