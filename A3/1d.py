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
	
	with open('fractionanno.gff3', 'a') as f:
		print('STATISTICS FOR ANNOTATED GENES: ', file = f)
		print('\ttotal genes: ', numgenes, file = f)
		print('\tperfect: {:.0%}'.format(numperfect/numgenes), file = f)
		print('\tstart only: {:.0%}'.format(numstart/numgenes), file = f)
		print('\tend only: {:.0%}'.format(numend/numgenes), file = f)
		print('\tneither: {:.0%}'.format(numneither/numgenes), file = f)
		print('The fraction of annotated genes on the positive strand that:', file = f)
		print('###Perfectly match predicted gene###', file = f)
		for name in names:
			if name in perfect:
				for i in perfect[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###Start only matches predicted gene###', file = f)
		for name in names:
			if name in start:
				for i in start[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###End only matches predicted gene###', file = f)
		for name in names:
			if name in end:
				for i in end[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###No match for predicted gene###', file = f)
		for name in names:
			if name in neither:
				for i in neither[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###', file = f)
	f.close()

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
	
	with open('fractionpred.gff3', 'a') as f:
		print('STATISTICS FOR PREDICTED GENES: ', file = f)
		print('\ttotal genes: ', numgenes, file = f)
		print('\tperfect: {:.0%}'.format(numperfect/numgenes), file = f)
		print('\tstart only: {:.0%}'.format(numstart/numgenes), file = f)
		print('\tend only: {:.0%}'.format(numend/numgenes), file = f)
		print('\tneither: {:.0%}'.format(numneither/numgenes), file = f)
		print('The fraction of predicted genes on the positive strand that:', file = f)
		print('###Perfectly match annotated gene###', file = f)
		for name in names:
			if name in perfect:
				for i in perfect[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###Start only matches annotated gene###', file = f)
		for name in names:
			if name in start:
				for i in start[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###End only matches annotated gene###', file = f)
		for name in names:
			if name in end:
				for i in end[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###No match for annotated gene###', file = f)
		for name in names:
			if name in neither:
				for i in neither[name]:
					print("%-20s %-5s %-5s %-5d %-5d %-5s %-5s %-5s %-5s" %(name, 'ena', 'CDS', i[0], i[1], '.', '+', '0', '.'), file = f)
		print('###', file = f)
	f.close()
	
def main(argv):
	real = argv[0]
	results = argv[1]
	annodict, names = readreal(real)
	mydict = readmyanno(results, names)
	fractionanno(annodict, mydict, names)
	fractionpred(annodict, mydict, names)

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['Vibrio_vulnificus.ASM74310v1.37.gff3', '1c.gff3'])
