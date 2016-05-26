#Max Parismony Phlogenetic Trees
from Bio import Phylo
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from ast import literal_eval
import re
import random
import cStringIO
import networkx as nx
import pylab
import os
import tempfile

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
            currentChar = tipMapping[tree[0]][position]
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
        minScore = float("inf")
        for char in ["A", "T", "C", "G"]:
            score = bestSinglePosition(tree, char, position, tipMapping, {})
            if score < minScore:
                minScore = score
        ans += minScore
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
    tree1 = RLRtoNewick(tree)
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
    """constructs a basic tree given a list of the leafs"""
    if len(tips) < 2:
        return "invalid dictionary"
    elif len(tips) == 2:
        return ("Anc", (tips[0], (), ()), (tips[1], (), ()))
    else:
        return ("Anc", (tips[0], (), ()), makeTree(tips[1:]))
        
        
def RLRtoPhyloxmlConverter(tree):
    """ converts from RLR format to Phyloxml"""
    tree = RLRtoNewick(tree)
    fp = tempfile.TemporaryFile()
    fp.write(str(tree))
    fp.write(b";")
    buf2 = cStringIO.StringIO()
    fp.seek(0)
    Phylo.convert(fp, 'newick', buf2, 'phyloxml')
    tree = buf2.getvalue()
  #  if type(tree) == str:
  #      print True
  #  print "HERE is tree:", tree
    return tree

def getGeneGroups(leaves):
    leaves = sorted(leaves)
    groups = []
    while leaves != []:
        firstLeaf = leaves.pop(0)
        group = [firstLeaf]
        toRemove = []
        for leaf in leaves:
            if firstLeaf in leaf:
                group.append(leaf)
                toRemove.append(leaf)
        for leaf in toRemove:
            leaves.remove(leaf)
        groups.append(group)
    return groups
        
                                       
def isFeasible(graph, leaves):
    geneGroups = getGeneGroups(leaves)
    groupPaths = []
    notFeasable = []
    feasable = True
    for group in geneGroups:
        groupName = group[0]
        pathNodes = set()
        for leaf in group[1:]:
            pathNodes = pathNodes.union(set(nx.shortest_path(graph, group[0], leaf, None)))
        groupPaths.append((groupName, pathNodes))
    for i in range(len(groupPaths)):
        geneName = groupPaths[i][0][:groupPaths[i][0].find('_')]
        for j in range(i+1, len(groupPaths)):
            if groupPaths[j][0].startswith(geneName + '_'):
                if groupPaths[i][1].intersection(groupPaths[j][1]) != set([]):
                    feasable = False
                    notFeasable.append((groupPaths[i][0], groupPaths[j][0]))
    return feasable
                
        
        
    

def NNIheuristic(FASTAFile, sampleSize):
    """"Find the maximum parsimony score for that tree"""
    # Import fasta alignment file
    myAlignment = AlignIO.read(FASTAFile, "fasta")
    
    # Create a tip mapping from the fasta file
    tipMapping = {}
    for record in myAlignment:
        tipMapping[record.id] = str(record.seq)
        
    # Compute a distance matrix and construct tree
    calculator = DistanceCalculator("identity") 
    myMatrix = calculator.get_distance(myAlignment)
    constructor = DistanceTreeConstructor()
    upgmaTree = constructor.upgma(myMatrix)
    Phylo.draw(upgmaTree)
    
    # Convert phyloxml tree to newick
    # biopython does not provide a function to do this so it was necessary
    # to write to a buffer in newick to convert then get rid of unneeded info
    for clade in upgmaTree.get_terminals():
        clade.name = "\"" + clade.name + "\""
    buf = cStringIO.StringIO()
    Phylo.write(upgmaTree, buf, 'newick', plain = True)
    tree = buf.getvalue()
    tree = re.sub(r'Inner\d*', '', tree)
    tree = tree.replace(";", "")
    tree = literal_eval(tree)    

    # RLR tree required for maxParsimony function
    tree = NewicktoRLR(tree)
    
    score = maxParsimony(tree, tipMapping)
    print score
    
    # Perform NNI heuristic
    while True:
        NNIs = allNNIs(tree)
        if len(NNIs)-1 < sampleSize:
            sampleSize = len(NNIs)-1
        toScore = random.sample(NNIs, sampleSize)
        # add feasibility test
        feasible = []
        for x in toScore:
            testTree = RLRtoPhyloxmlConverter(x)
            if isFeasible(testTree) == true:
                feasible.append(testTree)
        
        scoredList = map(lambda x: (maxParsimony(x, tipMapping), x), feasible)
        sortedlist = sorted(scoredList)
        if sortedlist[0][0] < score:
            score = sortedlist[0][0]
            tree = sortedlist[0][1]
	sortedlist = sorted(scoredList)
	if sortedlist[0][0] < score:
	       score = sortedlist[0][0]
	       tree = sortedlist[0][1]
        else:
            break            
    outputTree = RLRtoNewick(tree)
    print score
    return outputTree
    

	
testList1 = ["a_1", "a_1_1", "a_2", "a_2_1", "b_1", "b_1_1"]
testList2 = ["a_1", "a_1_1", "b_2", "b_2_1", "b_1", "b_1_1"]




