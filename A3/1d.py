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

def compareannomy(anno_dict, my_anno_dict, seq_name_list):
	perfect_anno_dict = {}
	start_anno_dict = {}
	end_anno_dict = {}
	neither_anno_dict = {}
	for name in seq_name_list:
		perfect_anno_dict[name] = []
		start_anno_dict[name] = []
		end_anno_dict[name] = []
		neither_anno_dict[name] = []

	for seq_name in anno_dict:
		for gene in anno_dict[seq_name]:
			#1.perfectly match
			if gene in my_anno_dict[seq_name]:
				perfect_anno_dict[seq_name].append(gene)
			elif gene[0] in [i[0] for i in my_anno_dict[seq_name]]:
				start_anno_dict[seq_name].append(gene)
			elif gene[1] in [i[1] for i in my_anno_dict[seq_name]]:
				end_anno_dict[seq_name].append(gene)
			else:
				neither_anno_dict[seq_name].append(gene)

	result_file = open('fraction_anno.gff3', 'w')
	result_file.write('###  Perfectly match both ends of one of my predicted genes\n')
	for seq_name in anno_dict:
		for gene in perfect_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.write('###  Match the start but not the end of a predicted gene\n')
	for seq_name in anno_dict:
		for gene in start_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.write('###  Match the end but not the start of a predicted gene\n')
	for seq_name in anno_dict:
		for gene in end_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.write('###  Do not match neither the start not the end of a predicted gene\n')
	for seq_name in anno_dict:
		for gene in neither_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.close()

def report_my_anno(anno_dict, my_anno_dict, seq_name_list):
	perfect_anno_dict = {}
	start_anno_dict = {}
	end_anno_dict = {}
	neither_anno_dict = {}
	for name in seq_name_list:
		perfect_anno_dict[name] = []
		start_anno_dict[name] = []
		end_anno_dict[name] = []
		neither_anno_dict[name] = []

	for seq_name in my_anno_dict:
		for gene in my_anno_dict[seq_name]:
			#1.perfectly match
			if gene in anno_dict[seq_name]:
				perfect_anno_dict[seq_name].append(gene)
			elif gene[0] in [i[0] for i in anno_dict[seq_name]]:
				start_anno_dict[seq_name].append(gene)
			elif gene[1] in [i[1] for i in anno_dict[seq_name]]:
				end_anno_dict[seq_name].append(gene)
			else:
				neither_anno_dict[seq_name].append(gene)

	result_file = open('my_fraction_anno.gff3', 'w')
	result_file.write('###  Perfectly match both ends of one of my predicted genes\n')
	for seq_name in anno_dict:
		for gene in perfect_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.write('###  Match the start but not the end of a predicted gene\n')
	for seq_name in anno_dict:
		for gene in start_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.write('###  Match the end but not the start of a predicted gene\n')
	for seq_name in anno_dict:
		for gene in end_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.write('###  Do not match neither the start not the end of a predicted gene\n')
	for seq_name in anno_dict:
		for gene in neither_anno_dict[seq_name]:
			result_file.write(seq_name+'\t ena\t CDS \t {}\t {}\t . \t +\t 0\t .\n'.format(gene[0], gene[1]))
	result_file.close()

def main(argv):
	real = argv[0]
	results = argv[1]
	annodict, names = readreal(real)
	mydict = readmyanno(results, names)
	print(annodict)
	compareannomy(annodict, mydict, names)
	report_my_anno(annodict, mydict, names)

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['Vibrio_vulnificus.ASM74310v1.37.gff3', '1c.gff3'])
