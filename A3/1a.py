import sys
import itertools

def readanno(file):
	anno_file = open(file,'r')
	lines = anno_file.readlines()
	
	lengthsum = 0
	lengthdict = {}
	start = 0;
	for line in lines:
		if "sequence-region" in line:
			data = line.strip().split()
			num = int(data[1][11:].lstrip("0"))
			lengthdict[num] = int(data[3])
			lengthsum += int(data[3])
		start += 1
		if "#!" in line:
			break

	# assumes only keep if CDS AND on + strand
	annodict = {new_list: [] for new_list in range(1, len(lengthdict)+1)} 
	idx = 0
	genesum = 0
	intercount = 0
	for line in lines:
		if not "#!" in line and idx > start:
			if not "###" in line:
				data = line.strip().split()
				if data[3] != 'Archive':
					if data[2] == 'CDS' and data[6] == '+':
						num = int(data[0][11:].lstrip("0"))
						start_end = (int(data[3]), int(data[4]))
						seq = annodict[num]
						#if overlap, consider as SINGLE gene
						if annodict[num]:
							prev = annodict[num][-1]
							if (prev[1] > start_end[0]):
								annodict[num][-1] = (prev[0], start_end[1])
								# allseq = annodict[num]
								# allseq[-1][1] = start_end[1]
								# annodict[num] = allseq
							else:
								annodict[num].append(start_end)
						else:
							annodict[num].append(start_end)
		idx+=1

	for i in annodict:
		if annodict[i]:
			# check start -> add inter after if no gene at start
			if (annodict[i][0][0]) != 1:
				intercount += 1
			# check end -> add inter after if no gene at end
			if (annodict[i][-1][1]) != lengthdict[i]:
				intercount+=1
			intercount += len(annodict[i]) - 1

			for k in annodict[i]:
				genesum += (k[1] - k[0] + 1)
		else:
			intercount += 1

	intersum = lengthsum - genesum

	genecount = 0
	for anno in annodict:
		genecount += len(annodict[anno])

	avginter = intersum/intercount
	avggene = genesum/genecount

	with open('configuration.txt', 'w') as f:
		print(avginter, file = f)
		print(avggene, file = f)
	f.close()
	
	return lengthdict, annodict

def readseq(file):
	seqfile = open(file, "r")
	lines = seqfile.readlines()

	num = 0
	seqdict = {}
	for line in lines:
		if ">" in line:
			data = line.strip().split()
			num = int(data[0][12:].lstrip("0"))
		else:
			data = line.strip('\n')
			if num in seqdict:
				seqdict[num] = seqdict[num] + data
			else:
				seqdict[num] = data
	return seqdict

def freq(annodict, seqdict, lengthdict):
	genedict = {new_list: '' for new_list in range(1, len(annodict)+1)} 
	interdict = {new_list: '' for new_list in range(1, len(annodict)+1)} 

	#split into inter and genic regions
	for i in range(1, len(genedict)+1):
		for k in range(0, len(annodict[i])):
			start = annodict[i][k][0] - 1 #because first idx = 0
			end = annodict[i][k][1] - 1 #because first idx = 0
		if annodict[i]:
			genedict[i] = genedict[i] + (seqdict[i][start:end])
			if k == 0:
				interdict[i] = interdict[i] + (seqdict[i][0:start])
			else:
				interdict[i] = interdict[i] + (seqdict[i][annodict[i][k-1][1]-1:start])
			if k == len(annodict[i])-1:
				interdict[i] = interdict[i] + (seqdict[i][end:lengthdict[i]])

	#initializing dictionaries
	nt = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
	# keywords = [''.join(i) for i in itertools.product(('A', 'C', 'T', 'G'), repeat = 3)] #gets 64 codons
	keywords = ['ATG', 'AAT','GCG','TTC','AAA','TTA','GTA','CAT','TCT','AAG','ATT','AAC','TTT','TTG','TCG',
				'TCC','ATC','GGT','TGG','GAC','CAA','GCC','GTC','CCC','AGT','GGC','GCT','GAA','CGC','GAG','CTT',
				'GTG','CTA','ACC','CAG','CCA','CTG','GAT','GCA','ACG','CTC','GTT','CCT','TGT','CAC','CGT','AGC',
				'ACA','ACT','GGG','TAC','TAT','CCG','CGA', 'TCA', 'CGG', 'ATA', 'TAA', 'GGA', 'TGC', 'AGA', 'AGG', 'TGA', 'TAG']
	codons = {}
	startcodons = []
	for i in keywords:
		codons[i] = 0

	#get codon and nt freq
	for i in range(1, len(genedict)+1):
		k = 0
		while k < len(genedict[i]) and len(genedict[i][k:k+3]) == 3:
			if k == 0:
				triplet = genedict[i][k:k+3]
				if not triplet in startcodons:
					startcodons.append(triplet)
			codons[genedict[i][k:k+3]] += 1
			k+=3

	# get nt frequency
	for i in range(1, len(interdict)+1):
		nt['A'] += interdict[i].count('A')
		nt['C'] += interdict[i].count('C')
		nt['G'] += interdict[i].count('G')
		nt['T'] += interdict[i].count('T')

	with open('configuration.txt', 'a') as f:
		print(nt, file = f)
		print(codons, file = f)
		print(startcodons, file = f)
	f.close()
	return

def main(argv):
	anno = argv[0]
	seq = argv[1]
	lengthdict, annodict = readanno(anno)
	seqdict = readseq(seq)
	freq(annodict, seqdict, lengthdict)

if __name__ == '__main__':
	main(sys.argv[1:])
	# main(['Vibrio_cholerae.GFC_11.37.gff3', 'Vibrio_cholerae.GFC_11.dna.toplevel.fa'])