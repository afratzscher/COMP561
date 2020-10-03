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
	# S_prime = np.zeros((m+1, n+1), dtype=object)
	# T_prime = np.zeros((m+1, n+1), dtype=object)
	# # S_prime[0][1] = 'AT'
	# # print(S_prime)
	# # S_prime[0][1] = 'T' + S_prime[0][1]
	# # print(S_prime)

	# stopped = []
	pointercopy = pointer.copy()
	# print(pointer)
	# num = 0
	# stopped.append((m,n))
	# loop = 1
	# flag = 0
	# while num > -1:
	# 	i = stopped[num][0]
	# 	j = stopped[num][1]
	# 	print(i,j)
	# 	while (i > 0 or j > 0):
	# 		repeats = len(pointer[i,j])
	# 		if (repeats > 1):
	# 			S_prime.append(S_prime[num].copy())
	# 			T_prime.append(T_prime[num].copy())
	# 			num += 1
	# 			print(num)
	# 			stopped.append((i, j))
	# 			print(stopped)
	# 			move = pointer[i,j][-1]
	# 			pointer[i,j] = pointer[i,j][:-1]
	# 		else:
	# 			move = pointer[i,j]
	# 		if move == 'D':
	# 			S_prime[num].insert(0, S[i-1])
	# 			T_prime[num].insert(0, T[j-1])
	# 			i -= 1
	# 			j -= 1
	# 		elif move == 'L':
	# 			S_prime[num].insert(0, "-")
	# 			T_prime[num].insert(0, T[j-1])
	# 			j -= 1
	# 		elif move == 'U':
	# 			S_prime[num].insert(0, S[i-1])
	# 			T_prime[num].insert(0, "-")
	# 			i -= 1
	# 	num -= 1

	iters = 0
	newFlag = False
	num = 0
	while iters >= 0:
		if num != 0:
			S_prime.append([])
			T_prime.append([])
			newFlag = True
		i = m
		j = n
		while (i > 0 or j > 0):
			if (len(pointer[i,j]) > 1):
				iters += len(pointer[i,j]) - 1
			if (len(pointer[i,j]) == 0):
				pointer[i,j] = pointercopy[i,j]
			move = pointer[i,j][-1]
			if move == 'D':
				if newFlag:
					S_prime[num].append(S[i-1])
					T_prime[num].append(T[j-1])
					newFlag = False
				else:
					S_prime[num].insert(0, S[i-1])
					T_prime[num].insert(0, T[j-1])	
				pointer[i,j] = pointer[i,j][:-1]
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
				pointer[i,j] = pointer[i,j][:-1]
				j -= 1
			elif move == 'U':
				if newFlag:
					S_prime[num].append(S[i-1])
					T_prime[num].append("-")
					newFlag = False
				else:
					S_prime[num].insert(0, S[i-1])
					T_prime[num].insert(0, "-")
				pointer[i,j] = pointer[i,j][:-1]
				i -= 1
		iters -=1
		num+=1
		
	
	
	for k in range(0, len(S_prime)):
		print ("S': ", S_prime[k])
		print(len(S_prime[k]))
		print("T': ", T_prime[k])
		print(len(T_prime[k]))
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
	main(['test3.fa', '1', '-1', '1'])
