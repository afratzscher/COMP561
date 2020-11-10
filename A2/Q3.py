import numpy as np
class leaf:
	def __init__(self, val, M, k):
		self.left = None
		self.right = None
		data = np.full((k, M+1), np.Inf)
		for i in range(0, len(val)):
			data[i][val[i]] = 0
		self.data = data

class Node:
	def __init__(self, u, v, M, k):
		self.left = u
		self.right = v
		data = np.full((k, M+1), np.Inf)
		self.data = data
		self.values = [np.Inf] * k

	def score(i):
		score(self, i)
	def setTraits(self):
		numk = len(self.data)
		for i in range(0, len(self.data)):
			scores = self.data[i].copy()
			idx = np.argmin(scores)
			self.values[i] = idx
		print(self.values)

def score(u, k, M):
	for i in range (0, M+1):
		v = u.left.data.copy()[k]
		w = u.right.data.copy()[k]
		vi = v[i]
		v[i] = np.Inf
		vmin = min(vi, min(v) + 1)
		wi = w[i]
		w[i] = np.Inf
		wmin = min(wi, min(w) + 1)
		u.data[k][i] = vmin + wmin
	print(u.data)
def main():
	M = 4
	n = 5
	k = 6
	s1 = (2, 1, 3, 3, 4, 1)
	s2 = (2, 1, 3, 1, 2, 1)
	s3 = (2, 3, 3, 3, 4, 1)
	s4 = (4, 3, 2, 1, 2, 3)
	s5 = (4, 3, 2, 3, 4, 1)
	l1 = leaf(s1, M, k)
	l2 = leaf(s2, M, k)
	l3 = leaf(s3, M, k)
	l4 = leaf(s4, M, k)
	l5 = leaf(s5, M, k)

	n1 = Node(l1,l2, M, k)
	n2 = Node(n1, l3, M, k)
	n3 = Node(l4, l5, M, k)
	n4 = Node(n2, n3, M, k)
	for i in range(0, k):
		score(n1, i, M)
		score(n2, i, M)
		score(n3, i, M)
		score(n4, i, M)
	n1.setTraits()
	n2.setTraits()
	n3.setTraits()
	n4.setTraits()
	
	

if __name__ == '__main__':
	main()