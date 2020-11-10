
from math import inf
from math import log


def read_config_anno(anno):
	#open and read the annotation file
	anno_file = open(anno,'r')
	anno_data = anno_file.readlines()

	#calculate the overall length of all sequences
	seq_length_sum = 0
	seq_length_dict = {}
	for i in range(1,152):
		line = anno_data[i].strip().split()
		seq_length_sum += int(line[3])
		seq_length_dict[i] = int(line[3])

	#dict store all anno data, eg.{1:[(s,e),(s,e),(s,e),(s,e)],2:[(s,e),(s,e),(s,e)]}
	anno_dict = {}
	seq_index = 1
	line_index = 158

	while seq_index<=151 and line_index<len(anno_data):

		#handle the 151-th seq annotation
		if line_index+2>=len(anno_data):
			anno_dict[seq_index] = current_seq
			break

		current_seq = []
		#loop within the same seq
		while not (anno_data[line_index].startswith('###') and anno_data[line_index+2].startswith('###')):
			#if at the beginning of current segment, record the data
			current_line = anno_data[line_index].split()
			if anno_data[line_index-1].startswith('###') and current_line[6] == '+':
				segment = (int(current_line[3]), int(current_line[4]))
				current_seq.append(segment)
			line_index += 1
		anno_dict[seq_index] = current_seq
		seq_index += 1
		line_index += 2

	seg_count = 0 
	for anno in anno_dict:
		seg_count += len(anno_dict[anno])

	return anno_dict, seq_length_sum, seg_count

def read_config_seq(seq):
	#open and read the sequence file
	seq_file = open(seq,'r')
	seq_data = seq_file.readlines()

	seq_index = 1
	line_index = 1 #start from the 2nd line
	
	#dict store all anno data, eg.{1:'CATGGAAGTA,2:'ACATA'}
	seq_dict = {}

	while seq_index<=151:
		current_seq = ''
		while not (line_index>=len(seq_data) or seq_data[line_index].startswith('>')):
			current_seq += seq_data[line_index].strip()
			line_index +=1 
		seq_dict[seq_index] = current_seq
		seq_index += 1
		line_index += 1
	return seq_dict

def write_config(anno_dict, seq_length_sum, seg_count, seq_dict):
	#calculate average intergenic and gene region length
	inter_length_sum = 0
	inter_number_sum = 0
	gene_length_sum = 0
	gene_number_sum = 0
	inter_nuc_dict = {'A':0, 'C':0, 'G':0, 'T':0}
	gene_codon_dict = {}

	for seq_index in anno_dict:
		inter_start_index = 0
		for seg in anno_dict[seq_index]:
			#count intergenic area
			for i in seq_dict[seq_index][inter_start_index:seg[0]-1]:
				inter_nuc_dict[i] += 1
			inter_start_index = seg[1]

			#count gene area
			for i in range(seg[0]-1, seg[1], 3):
				if seq_dict[seq_index][i:i+3] in gene_codon_dict:
					gene_codon_dict[seq_dict[seq_index][i:i+3]] += 1
				else:
					gene_codon_dict[seq_dict[seq_index][i:i+3]] = 1

			#calculate length
			gene_length_sum += seg[1] - seg[0]
			gene_number_sum += 1
			inter_number_sum += 1

		for i in seq_dict[seq_index][inter_start_index:]:
			inter_nuc_dict[i] += 1
		inter_number_sum += 1

	inter_length_sum = seq_length_sum - gene_length_sum

	inter_length_avg = inter_length_sum / inter_number_sum
	gene_length_avg = gene_length_sum / gene_number_sum
	
	#write info: inter_length_avg, gene_length_avg, inter_nuc_dict, gene_codon_dict
	config_file_name = "text.txt"
	config_file = open(config_file_name, "w")
	config_file.write(str(inter_length_avg) + '\n')
	config_file.write(str(gene_length_avg) + '\n')
	config_file.write(str(inter_nuc_dict) + '\n')
	config_file.write(str(gene_codon_dict) + '\n')
	config_file.close()

def read_config():
	config_file_name = "text.txt"
	config_file = open(config_file_name, "r")
	config = config_file.readlines()
	inter_length_avg = float(config[0].strip())
	gene_length_avg = float(config[1].strip())
	inter_nuc_dict = yaml.load(config[2].strip())
	gene_codon_dict = yaml.load(config[3].strip())

	i_to_s_prob = 1 / inter_length_avg
	m_to_t_prob = 3 / gene_length_avg

	inter_nuc_prob = {}
	inter_length = 0
	for nuc in inter_nuc_dict:
		inter_length += inter_nuc_dict[nuc]
	for nuc in inter_nuc_dict:
		inter_nuc_prob[nuc] = inter_nuc_dict[nuc] / inter_length

	gene_codon_prob = {}
	gene_codon_dict['ATG'] -= seg_count #ATG as start codon is removed for prob dist
	gene_length = 0
	for gene in gene_codon_dict:
		gene_length += gene_codon_dict[gene]
	for gene in gene_codon_dict:
		gene_codon_prob[gene] = gene_codon_dict[gene] / gene_length
	return i_to_s_prob, m_to_t_prob, inter_nuc_prob, gene_codon_prob

def read_seq():
	seq_file = open(sequence_file_name, "r")
	lines = seq_file.readlines()
	seq = ''
	seq_list = []
	seq_name_list = []
	for i in range(0,len(lines)):
		if lines[i].strip().startswith('>'):
			line = lines[i].strip().split()
			seq_name_list.append(line[0][1:])
			if i!=0:
				seq_list.append(seq)
				seq = ''
		else:
			seq += lines[i].strip()
	seq_list.append(seq)
	return seq_list, seq_name_list


def main(argv):
	anno = argv[0]
	seq = argv[1]
	anno_dict, seq_length_sum, seg_count = read_config_anno(anno)
	seq_dict = read_config_seq(seq)
	write_config(anno_dict, seq_length_sum, seg_count, seq_dict)
	i_to_s_prob, m_to_t_prob, inter_nuc_prob, gene_codon_prob = read_config()
	seq_list, seq_name_list = read_seq()

if __name__ == '__main__':
	# main(sys.argv[1:])
	# main(['test.gff3', 'Vibrio_cholerae.GFC_11.dna.toplevel.fa'])
	main(['Vibrio_cholerae.GFC_11.37.gff3', 'Vibrio_cholerae.GFC_11.dna.toplevel.fa'])