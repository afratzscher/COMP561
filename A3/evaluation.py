import sys
import config
import numpy as np

def readreal(file):
	anno_file = open(file,'r')
	lines = anno_file.readlines()
	
	names = []
	start = 0;
	for line in lines:
		if "sequence-region" in line:
			data = line.strip().split()
			names.append(data[1])
		start += 1
		if "#!" in line:
			break

	# assumes only keep if CDS AND on + strand
	annodict = {new_list: [] for new_list in names} 
	idx = 0
	for line in lines:
		if not "#!" in line and idx > start:
			if not "###" in line:
				data = line.strip().split()
				if data[3] != 'Archive':
					if data[2] == 'CDS' and data[6] == '+':
						name = data[0]
						start_end = (int(data[3]), int(data[4]))
						seq = annodict[name]
						#if overlap, consider as SINGLE gene
						if seq:
							prev = annodict[name][-1]
							if (prev[1] > start_end[0]):
								annodict[name][-1] = (prev[0], start_end[1])
							else:
								annodict[name].append(start_end)
						else:
							annodict[name].append(start_end)
		idx+=1
	
	return annodict, names

def readmyanno(file, names):
	anno_file = open(file,'r')
	lines = anno_file.readlines()

	# assumes only keep if CDS AND on + strand
	annodict = {new_list: [] for new_list in names} 
	idx = 0
	start = 0
	for line in lines:
		if not "#!" in line and idx > start:
			if not "###" in line:
				data = line.strip().split()
				if data[3] != 'Archive':
					if data[2] == 'CDS' and data[6] == '+':
						name = data[0]
						start_end = (int(data[3]), int(data[4]))
						seq = annodict[name]
						#if overlap, consider as SINGLE gene
						if seq:
							prev = annodict[name][-1]
							if (prev[1] > start_end[0]):
								annodict[name][-1] = (prev[0], start_end[1])
							else:
								annodict[name].append(start_end)
						else:
							annodict[name].append(start_end)
		idx+=1
	return annodict

def getseq(seq):
	seqfile = open(seq, "r")
	lines = seqfile.readlines()

	seqdict = {}
	for line in lines:
		if ">" in line:
			data = line.strip().split()
			name = data[0][1:]
		else:
			data = line.strip('\n')
			if name in seqdict:
				seqdict[name] = seqdict[name] + data
			else:
				seqdict[name] = data
	return seqdict

def fractionpred(annodict, mydict, names):
	perfect = {}
	start = {}
	end = {}
	neither = {}
	numgenes = 0
	numperfect = 0
	numstart = 0
	numend = 0
	numneither = 0

	for name in names:
		for i in mydict[name]:
			numgenes+=1
			if i in annodict[name]:
				if name in perfect:
					perfect[name].append(i)
				else:
					perfect[name] = [i]
				numperfect+=1
			elif i[0] in (i[0] for i in annodict[name]):
				if name in start:
					start[name].append(i)
				else:
					start[name] = [i]
				numstart+=1
			elif i[1] in (i[1] for i in annodict[name]):
				if name in end:
					end[name].append(i)
				else:
					end[name] = [i]
				numend+=1
			else:
				if name in neither:
					neither[name].append(i)
				else:
					neither[name] = [i]
				numneither+=1

	perfectmin = 1000000
	perfectmax = 0
	perfectlen = 0
	perfectcount = 0
	for name in names:
		if name in perfect:
			for i in perfect[name]:
				perfectcount += 1
				length = i[1] - i[0] + 1
				if length > perfectmax:
					perfectmax = length
				if length < perfectmin:
					perfectmin = length
				perfectlen += length
	perfectavg = perfectlen/perfectcount

	startmin = 1000000
	startmax = 0
	startlen = 0
	startcount = 0
	for name in names:
		if name in start:
			for i in start[name]:
				startcount += 1
				length = i[1] - i[0] + 1
				if length > startmax:
					startmax = length
				if length < startmin:
					startmin = length
				startlen += length
	startavg = startlen/startcount

	stopmin = 1000000
	stopmax = 0
	stoplen = 0
	stopcount = 0
	for name in names:
		if name in end:
			for i in end[name]:
				stopcount += 1
				length = i[1] - i[0] + 1
				if length > stopmax:
					stopmax = length
				if length < stopmin:
					stopmin = length
				stoplen += length
	stopavg = stoplen/stopcount

	nonemin = 1000000
	nonemax = 0
	nonelen = 0
	nonecount = 0
	for name in names:
		if name in neither:
			for i in neither[name]:
				nonecount += 1
				length = i[1] - i[0] + 1
				if length > nonemax:
					nonemax = length
				if length < nonemin:
					nonemin = length
				nonelen += length
	noneavg = nonelen/nonecount
	print(perfectmin, perfectmax, perfectavg)
	print(startmin, startmax, startavg)
	print(stopmin, stopmax, stopavg)
	print(nonemin, nonemax, noneavg)

	return perfect, start, end, neither

def startstopanalysis(cat, names, seq):
	startcnt = {'ATG': 0, 'TGG': 0, 'GTG': 0, 'other': 0}
	stopcnt = {'TAG': 0, 'TGA': 0, 'TAA': 0, 'other': 0}
	for name in names:
		if name in cat:
			for i in cat[name]:
				if (i[1]-i[0]) < 180:
					# print(name, i[0], i[1])
					start = seq[name][i[0]-1:i[0]+2] #-1 beacuse 1 is 0 in seq
					end = seq[name][i[1]:i[1]+3]
					if start in startcnt:
						startcnt[start] += 1
					else:
						startcnt['other'] += 1
					if end in stopcnt:
						stopcnt[end] +=1
					else:
						startcnt['other'] += 1
	print('START: ', startcnt)
	print('STOP: ', stopcnt)

def fractionanno(annodict, mydict, names):
	perfect = {}
	start = {}
	end = {}
	neither = {}
	numgenes = 0
	numperfect = 0
	numstart = 0
	numend = 0
	numneither = 0

	for name in names:
		for i in annodict[name]:
			numgenes+=1
			if i in mydict[name]:
				if name in perfect:
					perfect[name].append(i)
				else:
					perfect[name] = [i]
				numperfect+=1
			elif i[0] in (i[0] for i in mydict[name]):
				if name in start:
					start[name].append(i)
				else:
					start[name] = [i]
				numstart+=1
			elif i[1] in (i[1] for i in mydict[name]):
				if name in end:
					end[name].append(i)
				else:
					end[name] = [i]
				numend+=1
			else:
				if name in neither:
					neither[name].append(i)
				else:
					neither[name] = [i]
				numneither+=1

	perfectmin = 1000000
	perfectmax = 0
	perfectlen = 0
	perfectcount = 0
	for name in names:
		if name in perfect:
			for i in perfect[name]:
				perfectcount += 1
				length = i[1] - i[0] + 1
				if length > perfectmax:
					perfectmax = length
				if length < perfectmin:
					perfectmin = length
				perfectlen += length
	perfectavg = perfectlen/perfectcount

	startmin = 1000000
	startmax = 0
	startlen = 0
	startcount = 0
	for name in names:
		if name in start:
			for i in start[name]:
				startcount += 1
				length = i[1] - i[0] + 1
				if length > startmax:
					startmax = length
				if length < startmin:
					startmin = length
				startlen += length
	startavg = startlen/startcount

	stopmin = 1000000
	stopmax = 0
	stoplen = 0
	stopcount = 0
	for name in names:
		if name in end:
			for i in end[name]:
				stopcount += 1
				length = i[1] - i[0] + 1
				if length > stopmax:
					stopmax = length
				if length < stopmin:
					stopmin = length
				stoplen += length
	stopavg = stoplen/stopcount

	nonemin = 1000000
	nonemax = 0
	nonelen = 0
	nonecount = 0
	for name in names:
		if name in neither:
			for i in neither[name]:
				nonecount += 1
				length = i[1] - i[0] + 1
				if length > nonemax:
					nonemax = length
				if length < nonemin:
					nonemin = length
				nonelen += length
	noneavg = nonelen/nonecount
	print(perfectmin, perfectmax, perfectavg)
	print(startmin, startmax, startavg)
	print(stopmin, stopmax, stopavg)
	print(nonemin, nonemax, noneavg)

	return perfect, start, end, neither
	

def main(argv):
	real = argv[0]
	results = argv[1]
	fasta = argv[2]
	seq = getseq(fasta)
	annodict, names = readreal(real)
	mydict = readmyanno(results, names)
	print('For predicted genes')
	perfect, start, end, neither = fractionpred(annodict, mydict, names)
	print('perfect')
	startstopanalysis(perfect, names, seq)
	print('start')
	startstopanalysis(start, names, seq)
	print('end')
	startstopanalysis(end, names, seq)
	print('neither')
	startstopanalysis(neither, names, seq)
	print('For annotated genes')
	# perfect, start, end, neither = fractionanno(annodict, mydict, names)
	# print('perfect')
	# startstopanalysis(perfect, names, seq)
	# print('start')
	# startstopanalysis(start, names, seq)
	# print('end')
	# startstopanalysis(end, names, seq)
	# print('neither')
	# startstopanalysis(neither, names, seq)
	print(seq['contig_118'][19914-1:19973])

# 	(18717, 18812) contig_112
# (10328, 10423) contig_118
# (19914, 19973) contig_118
# (7864, 7962) contig_121
# (18068, 18154) contig_131
# (63192, 63287) contig_25
# (33500, 33598) contig_3
# (45334, 45420) contig_30
# (1388, 1477) contig_31
# (532, 596) contig_314
# (10917, 11009) contig_32
# (171, 232) contig_391
# (54613, 54708) contig_6
# (47203, 47301) contig_63
# (20029, 20112) contig_75
# (7688, 7783) contig_77
# (27117, 27188) contig_87
# (105755, 105832) contig_9
# (7773, 7850) contig_96


if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['Vibrio_vulnificus.ASM74310v1.37.gff3', '1c.gff3', 'Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa'])
