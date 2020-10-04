'''
FILE: NW.py
PURPOSE: practice implementing N-W pairwise sequence alignment algo
INPUT: none
OUTPUT: none
USAGE: python3 NW.py hw1_long.fa 1 -1 -1
'''
import sys
import numpy as np

def getSeq(fasta):
	first = True
	firstSeq = True
	seq1 = ""
	seq2 = ""
	with open(fasta) as f:
		while True:
			line = f.readline()
			if len(line) > 1:
				if line[0]==">":
					if firstSeq:
						seq1 = line[1:-1]
						firstSeq = False
					else:
						seq2 = line[1:-1]
				else:
					if first:
						S = line[:-1]
						first = False
					else:
						T = line
						break
	return seq1, seq2, S, T

def align(seqs, match, mismatch, pen):
	seq1 = seqs[0]
	seq2 = seqs[1]
	S = seqs[2]
	T = seqs[3]
	m = len(S)
	n = len(T)
	X = np.zeros((m+1,n+1))
	pointer = np.zeros((m+1, n+1), dtype=np.dtype('U3'))

	# initialize 0 column and row
	for j in range(0, n+1):
		X[0][j] = pen * j
		pointer[0,j] = 'L' #go to left
	for i in range(0, m+1):
		X[i][0] = pen * i
		pointer[i,0] = 'U' # go up
	pointer[0,0] = '0'

	# go through matrix
	for i in range(1, m+1):
		for j in range(1, n+1):
			diag = X[i-1, j-1]
			if S[i-1] == T[j-1]:
				diag += match
			else:
				diag += mismatch
			maxval = max(diag, X[i][j-1] + pen, X[i-1][j] + pen)
			if maxval == diag:
				pointer[i,j] = (pointer[i,j]+("D"))
			if maxval == (X[i][j-1] + pen):
				pointer[i,j] = (pointer[i,j]+("L"))
			if maxval == (X[i-1][j] + pen):
				pointer[i,j] = (pointer[i,j]+("U"))
			X[i,j] = maxval
	
	# start traceback
	S_prime = ['']
	T_prime = ['']
	pointercopy = pointer.copy()

	iters = 1
	num = 0 
	prev = []
	branch = []
	newFlag = False
	while iters > 0:
		if num != 0: # add new list for new branch
			S_prime.append('')
			T_prime.append('')
			newFlag = True
		i = m
		j = n
		lastS = ''
		lastT = ''
		# do for each branch
		while (i > 0 or j > 0):
			if (len(pointer[i,j]) > 1):
				iters += len(pointer[i,j]) - 1
				branch.append((i,j))
				move = pointer[i,j][-1]
			else:
				if (len(pointer[i,j]) == 0):
					pointer[i,j] = pointercopy[i,j][0]
				move = pointer[i,j]
			if move == 'D':
				if newFlag:
					S_prime[num] = S[i-1] + S_prime[num]
					T_prime[num]  = T[j-1] + T_prime[num]
					newFlag = False
				else:
					S_prime[num] = S[i-1] + S_prime[num]
					T_prime[num] = T[j-1] + T_prime[num]
				i -= 1
				j -= 1
				lastS = S[i-1]
				lastT = T[j-1]
			elif move == 'L':
				if (lastS != '-'):
					if newFlag:
						S_prime[num] = "-" + S_prime[num]
						T_prime[num] = T[j-1] + T_prime[num]
						newFlag = False
					else:
						S_prime[num] = "-" + S_prime[num]
						T_prime[num] = T[j-1] + T_prime[num]
					lastS = "-"
					lastT = T[j-1]
					j -= 1
				else:
					S_prime[num] = ""
					T_prime[num] = ""
					i = 0
					j = 0
			elif move == 'U':
				if (lastT != '-'):
					if newFlag:
						S_prime[num] = S[i-1] + S_prime[num]
						T_prime[num] = "-" + T_prime[num]
						newFlag = False
					else:
						S_prime[num]  = S[i-1] + S_prime[num]
						T_prime[num] = "-" + T_prime[num]
					lastS = S[i-1]
					lastT = "-"
					i -= 1
				else:
					S_prime[num] = ""
					T_prime[num] = ""
					i = 0
					j = 0
		
		# print only first found optimal (for speed for hw1_long.fa)
		if (S_prime[num] != ''):
			print(seq1 + ": " + S_prime[num])
			print(seq2 + ": ", T_prime[num])
			print("Score: ", X[m,n])
			break;
		#comment above out if want ALL backtrack solutions

		iters -=1
		num+=1

		if branch:
			if prev:
				if prev[-1] < branch[-1]:
					pointer[prev[-1][0], prev[-1][1]] = pointercopy[prev[-1][0], prev[-1][1]]
					prev = prev[:-1]
			pointer[branch[-1][0], branch[-1][1]] = pointer[branch[-1][0], branch[-1][1]][:-1]
			prev.append(branch[-1])
			branch = branch[:-1]
		
		if pointer[m,n] == '':
			break

	# if want to print all
	# for k in range(0, len(S_prime)):
	# 	if not (S_prime[k] == ''):
	# 		print("S': ", S_prime[k])
	# 		print("T': ", T_prime[k])
	# 		print("SCORE: ",  X[m,n])
	# 		print("\n")

def main(argv):
	fasta = argv[0]
	match = int(argv[1])
	mismatch = int(argv[2])
	penalty = int(argv[3])
	info = getSeq(fasta)
	align(info, match, mismatch, penalty)
	print("DONE")

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['long2.fa', '1', '-1', '-1'])
