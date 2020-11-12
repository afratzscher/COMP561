
import argparse
from math import inf
from math import log

def argument_parse():
	ap = argparse.ArgumentParser()
	ap.add_argument('sequence_file')
	ap.add_argument('config_file')
	args = vars(ap.parse_args())
	sequence_file_name = args['sequence_file']
	config_file_name = args['config_file']
	print('Sequence file name:', sequence_file_name)
	print('config file name:', config_file_name)
	return sequence_file_name, config_file_name


def read_config():
	config_file_name  = 'configuration.txt'
	config_file = open(config_file_name, "r")
	config = config_file.readlines()
	inter_length_avg = float(config[0].strip())
	gene_length_avg = float(config[1].strip())
	inter_nuc_dict = eval(config[2].strip())
	gene_codon_dict = eval(config[3].strip())

	i_to_s_prob = 1 / inter_length_avg
	m_to_t_prob = 3 / gene_length_avg

	inter_nuc_prob = {}
	inter_length = 0
	for nuc in inter_nuc_dict:
		inter_length += inter_nuc_dict[nuc]
	for nuc in inter_nuc_dict:
		inter_nuc_prob[nuc] = inter_nuc_dict[nuc] / inter_length

	gene_codon_prob = {}
	# gene_codon_dict['ATG'] -= seg_count #ATG as start codon is removed for prob dist
	gene_length = 0
	for gene in gene_codon_dict:
		gene_length += gene_codon_dict[gene]
	for gene in gene_codon_dict:
		gene_codon_prob[gene] = gene_codon_dict[gene] / gene_length
	return i_to_s_prob, m_to_t_prob, inter_nuc_prob, gene_codon_prob

def read_seq():
	sequence_file_name = 'Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa'
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
	print(seq)

	return seq_list, seq_name_list

def initialize():
	initial_prob = {'I':1, 'S':0, 'M':0, 'T':0}
	transition_prob = {('I','I'): 1-i_to_s_prob, ('I','S'):i_to_s_prob, ('I','M'):0, ('I','T'):0,
				('S','I'):0, ('S','S'):0, ('S','M'):1, ('S','T'):0,
				('M','I'): 0, ('M','S'):0, ('M','M'):1-m_to_t_prob, ('M','T'):m_to_t_prob,
				('T','I'): 1, ('T','S'):0, ('T','M'):0, ('T','T'):0}
	i_emission_prob = inter_nuc_prob.copy()
	s_emission_prob = gene_codon_prob.copy()
	for s in s_emission_prob:
		if s=='ATG':
			s_emission_prob[s] = 1
		else:
			s_emission_prob[s] = 0

	m_emission_prob = gene_codon_prob.copy()
	prob_sum = 0
	for m in m_emission_prob:
		if m!='TAA' and m!='TAG' and m!='TGA':
			prob_sum += m_emission_prob[m]
		else:
			m_emission_prob[m] = 0
	for m in m_emission_prob:
		m_emission_prob[m] /= prob_sum

	t_emission_prob = gene_codon_prob.copy()
	for t in t_emission_prob:
		if t=='TAA' or t=='TAG' or t=='TGA':
			t_emission_prob[t] = 1/3
		else:
			t_emission_prob[t] = 0
	emission_prob = {'I':i_emission_prob, 'S':s_emission_prob, 'M':m_emission_prob, 'T':t_emission_prob}
	return initial_prob, transition_prob, emission_prob

def viterbi(seq):
	take_codon = False
	codon_flag = 1
	current_prob = {'I':emission_prob['I'][seq[0]], 'S':-inf, 'M':-inf, 'T':-inf} 
	path_trellis = [{'I':'I','S':'I','M':'I','T':'I'}]

	for ch in range(1, len(seq)-2):
		current_path = {}
		last_prob = current_prob.copy()

		if take_codon == True:
			codon_flag += 1
			if codon_flag > 3:
				codon_flag = 1
		if take_codon == False or codon_flag==1:
			codon = seq[ch:ch+3]
			codon_flag = 1
			if emission_prob['S'][codon]>0:
				take_codon = True

		for curr_state in initial_prob:
			prob_list = {}
			if curr_state == 'I':
				for prev_state in initial_prob:
					if transition_prob[(prev_state,'I')]>0 and emission_prob['I'][seq[ch]]>0:
						prob_list[prev_state] = last_prob[prev_state] + log(transition_prob[(prev_state,'I')]) + log(emission_prob['I'][seq[ch]])
					else:
						prob_list[prev_state] = -inf
				current_prob['I'] = max(prob_list.values())
				current_path['I'] = max(prob_list, key = prob_list.get)
		
		
			if curr_state=='S' or curr_state=='M' or curr_state=='T':
				if take_codon == True:
					if codon_flag==1:
						for prev_state in initial_prob:
							if transition_prob[(prev_state,curr_state)]>0 and emission_prob[curr_state][codon]>0:
								prob_list[prev_state] = last_prob[prev_state] + log(transition_prob[(prev_state,curr_state)]) + log(emission_prob[curr_state][codon])
							else:
								prob_list[prev_state] = -inf

						current_prob[curr_state] = max(prob_list.values())
						current_path[curr_state] = max(prob_list, key = prob_list.get)
						
					elif codon_flag==2 or codon_flag==3:
						current_path[curr_state] = curr_state
				else:
					current_prob[curr_state] = -inf
		if emission_prob['T'][codon]>0 and codon_flag==3:
			take_codon = False
		path_trellis.append(current_path)

	return path_trellis, current_prob

def trace_back(seq_name):
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

if __name__=="__main__":
	#usage: python predict.py Vibrio_vulnificus.ASM74310v1.dna.nonchromosomal.fa config.txt
	i_to_s_prob, m_to_t_prob, inter_nuc_prob, gene_codon_prob = read_config()
	seq_list, seq_name_list = read_seq()
	initial_prob, transition_prob, emission_prob = initialize()
	result_file = open('result.gff3', 'w')
	print('{} sequences in total'.format(len(seq_list)))
	for i in range(len(seq_list)):
		print('searching {}-th sequence'.format(i))
		path_trellis, current_prob = viterbi(seq_list[i])
		trace_back(seq_name_list[i])
		break
	result_file.write('###\n')
	result_file.close()