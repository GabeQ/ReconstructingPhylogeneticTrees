from threading import Thread
import Queue

def isLeaf(tree):
    """is the tree a leaf, returns true or false"""
    if tree[1] == () and tree[2] == ():
        return True
    else:
        return False

def isDifferent(char1, char2):
    """determines if we need to add a 1 or a 0 to the score of a letter """
    if char1 == char2:
        return 0
    else:
        return 1

def bestSinglePosition(tree, character, position, tipMapping, memo):
    """Returns the minimum parsimony score for a given position at a given nucleotide position"""
    if (tree, character) in memo:
        return memo[(tree, character)]
    else:
        if isLeaf(tree) == True:
            currentChar = tipMapping[tree[0]][position] # this would find the leaf in the dictionary, and the return the letter of the sequence of that leaf at the given position.
            if isDifferent(character, currentChar) == 1:
                return float("inf")
            else:
                return 0
        else:
            leftTree = tree[1]
            rightTree = tree[2]
            minL = float("inf")
            for leftChar in ["A", "T", "C", "G"]:
                ansL = bestSinglePosition(leftTree, leftChar, position, tipMapping, memo) + isDifferent(leftChar, character)
                if ansL < minL:
                    minL = ansL
            minR = float("inf")
            for rightChar in ["A", "T", "C", "G"]:
                ansR = bestSinglePosition(rightTree, rightChar, position, tipMapping, memo) + isDifferent(rightChar, character)
                if ansR < minR:
                    minR = ansR 
            ans = minR + minL
            memo[(tree, character)] = ans
            return ans

def maxParsimony(tree, tipMapping):
    """computes the best score for all positions in a sequence and for all possible characters"""
    keys = tipMapping.keys()
    length = len(tipMapping[keys[0]])
    ans = 0
    for position in range(length):
        threads = []
        minScore = float("inf")
        que = Queue.Queue()
        for char in ["A", "T", "C", "G"]:
            args = (tree, char, position, tipMapping, {})
            thread = Thread(target=lambda q, arg: q.put(bestSinglePosition(args[0],args[1],args[2],args[3],args[4])), args=(que, args))
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join()
        while not que.empty():
            score = que.get()
            if score < minScore:
                minScore = score
        ans += minScore
    return ans