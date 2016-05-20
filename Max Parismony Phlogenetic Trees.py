#Max Parismony Phlogenetic Trees


def isLeaf(tree):
	if tree[1] == () and tree[2] == ():
		return True
	else:
		return False


def different(char, char1):
	""" determines if we need to add a 1 or a 0 to the score of a letter """

	if char == char1:
		return 0
	else:
		return 1


def bestSinglePosition(tree, character, position, tipMapping, memo):
	"""memo is a dictionary that will be used for memoization"""

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
