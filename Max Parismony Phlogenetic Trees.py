Max Parismony Phlogenetic Trees


def isLeaf(tree):
	if tree[1] and tree[2] == ():
		return true
	else:
		return false


def different(char, char1):
	""" determines if we need to add a 1 or a 0 to the score of a letter """

	if char == char1:
		return 0
	else:
		return 1


def bestSinglePosition(tree, character, position, tipMapping, memo):
	"""
	memo is a dictionary that will be used to support memoization

	"""
	if isLeaf(tree) == true:
		char1 = tipMapping[tree[0]][position] # this would find the leaf in the dictionary, and the return the letter of the sequence of that leaf at the given position.
		num1 = different(character, char1)
		if num1 == 1:
			num = float(inf)
		return num

	else:
		root = tree[0]
		leftTree = tree[1]
		rightTree = tree[2]

		miniL = float(inf)
		for leftChar in ["A", "T", "C", "G"]:
			ansL = different(leftChar, character) + bestSinglePosition(leftTree, leftChar, position, tipMapping, memo)
			if ansL < miniL:
				miniL = ansL

		miniR = float(inf)
		for rightChar in ["A", "T", "C", "G"]:
			ansR = different(rightChar, character) + bestSinglePosition(rightTree, rightChar, position, tipMapping, memo)
			if ansR < miniR:
				miniR = ansR

		return miniR + miniL










