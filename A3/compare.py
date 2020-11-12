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
	return seq_list

def myread_seq():
	seqfile = open('Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa', "r")
	lines = seqfile.readlines()


	num = 0
	seqdict = {}
	idx = 0
	for line in lines:
		if ">" in line:
			data = line.strip().split()
			# num = int(data[0][12:].lstrip("0"))
			num = int(data[0][8:].lstrip("0"))
			# config.__NAMES__.append(data[0][1:])
		else:
			data = line.strip('\n')
			if num in seqdict:
				seqdict[num] = seqdict[num] + data
			else:
				seqdict[num] = data
	# config.__SEQ__ = seqdict
	return seqdict

def main():
	my = myread_seq()
	other = read_seq()
	print(other[0])
	print(my[100] == other[100])
	shared_items = {k: my[k] for k in my if k in other and my[k] == other[k]}
	print(len(shared_items))

if __name__ == '__main__':
	main()

