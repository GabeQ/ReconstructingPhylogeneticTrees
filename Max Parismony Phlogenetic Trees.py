#Max Parismony Phlogenetic Trees

import random

def isLeaf(tree):
	"""is the tree a leaf, returns true or false"""
	if tree[1] == () and tree[2] == ():
		return True
	else:
		return False

def different(char, char1):
	"""determines if we need to add a 1 or a 0 to the score of a letter """
	if char == char1:
		return 0
	else:
		return 1

def bestSinglePosition(tree, character, position, tipMapping, memo):
	"""Returns the minimum parsimony score for a given position at a given nucleotide position"""
	if tree in memo:
		return memo[tree]
	else: 
		if isLeaf(tree) == True:
			char1 = tipMapping[tree[0]][position] # this would find the leaf in the dictionary, and the return the letter of the sequence of that leaf at the given position.
			num1 = different(character, char1)
			if num1 == 1:
				num1 = float("inf")
			return num1
		else:
			root = tree[0]
			leftTree = tree[1]
			rightTree = tree[2]
			miniL = float("inf")
			for leftChar in ["A", "T", "C", "G"]:
				ansL = bestSinglePosition(leftTree, leftChar, position, tipMapping, memo) + different(leftChar, character)
				if ansL < miniL:
					miniL = ansL
			memo[leftTree[0]] = miniL
			miniR = float("inf")
			for rightChar in ["A", "T", "C", "G"]:
				ansR = bestSinglePosition(rightTree, rightChar, position, tipMapping, memo) + different(rightChar, character)
				if ansR < miniR:
					miniR = ansR 
			memo[rightTree[0]] = miniR
			return miniR + miniL

def maxParsimony(tree, tipMapping):
	"""computes the best score for all positions in a sequence and for all possible characters"""
	list1 = tipMapping.keys()
	length = len(tipMapping[list1[0]])
	ans = 0
	for position in range(length):
		minNum = float("inf")
		for char in ["A", "T", "C", "G"]:
			num1 = bestSinglePosition(tree, char, position, tipMapping, {})
			if num1 < minNum:
				minNum = num1
		ans += minNum
	return ans

def RLRtoNewick(tree):
	"""converts from Root,Left,Right format trees to Newick format trees """
	if isLeaf(tree) == True:
		return tree[0]
	else:
		return (RLRtoNewick(tree[1]), RLRtoNewick(tree[2]))

def NewicktoRLR(tree):
	"""converts from Newick format to RLR"""
	if type(tree[0]) != tuple and type(tree[1]) != tuple:
		return ('Anc', (tree[0], (), ()), (tree[1], (), ()))
	elif type(tree[0]) != tuple:
		return ('Anc', (tree[0], (), ()), NewicktoRLR(tree[1]))
	elif type(tree[1]) != tuple:
		return ('Anc', NewicktoRLR(tree[0]), (tree[1], (), ()))
	else:
		return ('Anc', NewicktoRLR(tree[0]), NewicktoRLR(tree[1]))

def goLeft(tree):
	"""goes to the left side of the given tree to find all possible rerooted trees"""
	LL = tree[0][0]
	LR = tree[0][1]
	R = tree[1]
	if type(LL) != tuple and type(LR) != tuple:
		return []
	elif type(LL) != tuple:
		return [(LR, (LL, R))] + goLeft((LR, (LL, R)))
	elif type(LR) != tuple:
		return [(LL, (LR, R))] + goLeft((LL, (LR, R)))
	else:
		return [(LL, (LR, R))] + goLeft((LL, (LR, R))) + [(LR, (LL, R))] + goLeft((LR, (LL, R)))

def goRight(tree):
	"""goes to the right side of the given tree to find all possible rerooted trees"""
	RL = tree[1][0]
	RR = tree[1][1]
	L = tree[0]
	if type(RL) != tuple and type(RR) != tuple:
		return []
	elif type(RL) != tuple:
		return [((L, RL), RR)] + goRight(((L, RL), RR))
	elif type(RR) != tuple:
		return [((L, RR), RL)] + goRight(((L, RR), RL))
	else:
		return [((L, RL), RR)] + goRight(((L, RL), RR)) + [((L, RR), RL)] + goRight(((L, RR), RL))

def go(tree):
	"""outputs a list of all the rerooted trees using goLeft and goRight"""
	return goLeft(tree) + goRight(tree)

def allNNIs(tree):
	""" "takes a root-left-right tree as input, finds all of the re-rootings of this tree, 
	and then returns a list of all of the NNI trees for those re-rootings in root-left-right format" """
	tree1 = RLRtoNewick(tree) #convert tree to Newick format
	reRooted = go(tree1)
	NNIs = []
	for tree in reRooted:
		LL = tree[0][0]
		LR = tree[0][1]
		RL = tree[1][0]
		RR = tree[1][1]
		NNIs.append(NewicktoRLR(((LL, RL), (LR, RR))))
		NNIs.append(NewicktoRLR(((LL, RR), (RL, LR))))
	return NNIs

def makeTree(tips):
	if len(tips) < 2:
		return "invaid dictionary"
	elif len(tips)==2:
		return ("Anc", (tips[0], (), ()), (tips[1], (), ()))
	else:
		return ("Anc", (tips[0], (), ()), makeTree(tips[1:]))

def NNIheuristic(tipMapping, sampleSize):
	""""Find the maximum parsimony score for that tree" """
	tips = tipMapping.keys()
	tree = makeTree(tips)
	score = maxParsimony(tree, tipMapping)
	iteration = 0
	while True:
		NNIs = allNNIs(tree)
		if len(NNIs)-1 < sampleSize:
			sampleSize = len(NNIs)-1
		toScore = random.sample(NNIs, sampleSize)
		count = 0
		for sample in toScore:
			newScore = maxParsimony(sample, tipMapping)
			if newScore < score:
				score = newScore
				tree = sample
			else:
				count += 1
			if count == sampleSize:
				iteration += 1
				break
		if iteration > 2:
			break
	outputTree = RLRtoNewick(tree)
	print score
	return outputTree
























