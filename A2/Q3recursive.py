import numpy as np
class leaf:
	def __init__(self, val, name, k, M):
		self.left = None
		self.right = None
		self.data = np.full((k, M+1), np.Inf)
		self.val = val
		self.name = name

class Node:
	def __init__(self, u, v, name, k, M):
		self.left = u
		self.right = v
		self.data = np.full((k, M+1), np.Inf)
		self.name = name

def score(u, i, M, k):
	if (u.left is not None): #not leaf
		score(u.left, i, M, k)
	if (u.right is not None):
		score(u.right, i, M, k)
	if (u.left is None and u.right is None): # leaf
		val = u.val
		u.data[i][val[i]] = 0
		return
	
	for m in range(0, M+1):
		minv = u.left.data[i][0] + abs(m - 0)
		minw = u.right.data[i][0] + abs(m-0)
		for n in range(0, k+1):
			minv = min(minv, (u.left.data[i][n] + abs(m-n)))
			minw = min(minw, (u.right.data[i][n] + abs(m-n)))
		u.data[i][m] = minv + minw
	return

def main():
	M = 6
	N = 5
	k = 6
	s1 = (3, 1, 3, 3, 4, 1)
	s2 = (2, 1, 3, 3, 4, 1)
	s3 = (2, 1, 3, 3, 4, 1)
	s4 = (3, 1, 3, 3, 4, 1)
	s5 = (3, 1, 3, 3, 4, 1)
	# s2 = (2, 1, 3, 1, 2, 1)
	# s3 = (2, 3, 3, 3, 4, 1)
	# s4 = (4, 3, 2, 1, 2, 3)
	# s5 = (4, 3, 2, 3, 4, 1)
	l1 = leaf(s1, 's1', k, M)
	l2 = leaf(s2, 's2', k, M)
	l3 = leaf(s3, 's3', k, M)
	l4 = leaf(s4, 's4', k, M)
	l5 = leaf(s5, 's5', k, M)

	n1 = Node(l1,l2, 'n1', k, M)
	n2 = Node(n1, l3, 'n2', k, M)
	n3 = Node(l4, l5, 'n3', k, M)
	n4 = Node(n2, n3, 'n4', k, M)

	result = [0, 0, 0, 0, 0, 0]
	d = [100, 100, 100, 100, 100, 100]
	for i in range(0, k):
		score(n4, i, M, k)
		result[i] = min(n4.data[i]) # pars score

		d[i] = n4.data[i].argmin()
	# print(result)
	print(d)



	
	

if __name__ == '__main__':
	main()