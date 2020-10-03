'''
FILE: NW.py
PURPOSE: practice implementing N-W pairwise sequence alignment algo
INPUT: none
OUTPUT: none
USAGE: python3 NW.py hw1_long.fa 1 -1 2
'''
import sys
import numpy as np

def getSeq(fasta):
	first = True
	with open(fasta) as f:
		while True:
			line = f.readline()
			if len(line) > 0 and line[0]==">":
				line = f.readline() 
				if first:
					S = line[:-1]
					first = False
				else:
					T = line
					break
	return S, T

def align(seqs, match, mismatch, pen):
	S = seqs[0]
	T = seqs[1]
	m = len(S)
	n = len(T)
	X = np.zeros((m+1,n+1))
	pointer = np.zeros((m+1, n+1), dtype=np.dtype('U3'))

	# initialize 0 column and row
	for j in range(0, n+1):
		X[0][j] = -pen * j
		pointer[0,j] = 'L' #go to left
	for i in range(0, m+1):
		X[i][0] = -pen * i
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
			maxval = max(diag, X[i][j-1] - pen, X[i-1][j] - pen)
			if maxval == diag:
				pointer[i,j] = (pointer[i,j]+("D"))
			if maxval == (X[i][j-1] - pen):
				pointer[i,j] = (pointer[i,j]+("L"))
			if maxval == (X[i-1][j] - pen):
				pointer[i,j] = (pointer[i,j]+("U"))
			X[i,j] = maxval
	

	# start traceback
	S_prime = [[]]
	T_prime = [[]]
	pointercopy = pointer.copy()

	iters = 1
	num = 0 
	prev = []
	branch = []
	newFlag = False
	while iters > 0:
		if num != 0: # add new list for new branch
			S_prime.append([])
			T_prime.append([])
			newFlag = True
		i = m
		j = n
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
					S_prime[num].append(S[i-1])
					T_prime[num].append(T[j-1])
					newFlag = False
				else:
					S_prime[num].insert(0, S[i-1])
					T_prime[num].insert(0, T[j-1])	
				i -= 1
				j -= 1
			elif move == 'L':
				if newFlag:
					S_prime[num].append("-")
					T_prime[num].append(T[j-1])
					newFlag = False
				else:
					S_prime[num].insert(0, "-")
					T_prime[num].insert(0, T[j-1])
				j -= 1
			elif move == 'U':
				if newFlag:
					S_prime[num].append(S[i-1])
					T_prime[num].append("-")
					newFlag = False
				else:
					S_prime[num].insert(0, S[i-1])
					T_prime[num].insert(0, "-")
				i -= 1
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

	for k in range(0, len(S_prime)):
		print(S_prime[k])
		print(T_prime[k])
		print("SCORE: ",  X[m,n])
		print("\n")

def main(argv):
	fasta = argv[0]
	match = int(argv[1])
	mismatch = int(argv[2])
	penalty = int(argv[3])
	seqs = getSeq(fasta)
	align(seqs, match, mismatch, penalty)
	print("DONE")

if __name__ == '__main__':
	# main(sys.argv[1:])
	main(['hw1_long.fa', '1', '-1', '1'])
