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
import uuid
import time
import datetime

#from threading import Thread
#import Queue

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
'''           
def maxParsimony(tree, tipMapping):
    """computes the best score for all positions in a sequence and for all possible characters - threaded version"""
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
'''    
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
    """converts from Root,Left,Right (RLR) format trees to Newick format trees """
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
    if type(LL) != tuple and type(LR) != tuple: #tuples indicate that we havent hit a leaf yet in Newick Format
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
    """ input: root-left-right tree
    finds all of the re-rootings of this tree, 
    returns: list of of the NNI re-rootings in RLR format"""
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

def getGeneGroups(leaves):
    '''takes a list of the leaves (consisting of genes with locus numbers and individual numbers) and 
    creates a list of lists grouping same genes together with different loci'''
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
    '''checks if the tree is possible by taking a graph (of the tree using makeGraph) and constructing a graph connecting 
    genes that have overlapping paths within the tree; if there is a gene that has multiple loci connected within this 
    graph (meaning at least two loci from the same gene can be reached within the graph), the tree is considered infeasible'''
    geneGroups = getGeneGroups(leaves)
    groupPaths = []
    intersectGraph = nx.Graph()
    for group in geneGroups:
        groupName = group[0]
        pathNodes = set()
        for leaf in group[1:]:
            pathNodes = pathNodes.union(set(nx.shortest_path(graph, group[0], leaf, None)))
        groupPaths.append((groupName, pathNodes))
    for i in range(len(groupPaths)):
        geneLocus1 = groupPaths[i][0]
        for j in range(i+1, len(groupPaths)):
            geneLocus2 = groupPaths[j][0]
            if groupPaths[i][1].intersection(groupPaths[j][1]) != set([]):
                intersectGraph.add_edge(geneLocus1, geneLocus2)
    for i in range(len(geneGroups)):
        for j in range(i+1, len(geneGroups)):
            geneLocus1 = geneGroups[i][0]
            geneLocus2 = geneGroups[j][0]
            if geneLocus1[:groupPaths[i][0].find('_')] == geneLocus2[:groupPaths[i][0].find('_')]:
                if intersectGraph.has_node(geneLocus1) and intersectGraph.has_node(geneLocus2) \
                   and nx.has_path(intersectGraph, geneLocus1, geneLocus2):
                    return False
    return True
                  
def makeGraph(graph, tree): 
    '''takes RLR tree and recursively constructs a networkx graph'''
    if isLeaf(tree):
        return tree[0]
    else:
        uniqueId = str(uuid.uuid1())
        graph.add_edge(uniqueId, makeGraph(graph, tree[1]))
        graph.add_edge(uniqueId, makeGraph(graph, tree[2]))
        return uniqueId
        
def getLeaves(tree): 
    '''takes RLR tree and gets a list of leaves'''
    if isLeaf(tree):
        return [tree[0]]
    else:
        return getLeaves(tree[1]) + getLeaves(tree[2])
        
def getEdges(tree):
    '''takes a RLR tree and gets the edges of the tree'''
    if isLeaf(tree):
        return []
    elif isLeaf(tree[1]) and isLeaf(tree[2]):
        return [(tree[0], tree[1][0]), (tree[0], tree[2][0])]
    else: return [(tree[0], tree[1][0]), (tree[0], tree[2][0])] + getEdges(tree[1]) + getEdges(tree[2])

def NNIHeuristic(FASTAFile, sampleSize, threshold, outputDir):
    """"takes a FASTAFile, constructs a UPGMA Tree from the file data, converts this tree to RLR format,
    tests if the tree is feasible and scores it using maxParsimony, and then goes into a while loop constructing
    NNIs, tests to see if the trees are feasible or infeasible; tries to find the tree with the lowest parsimony
    score that is also feasible"""
    random.seed(0)
    outputFile = FASTAFile.replace(".align", ".out")
    if "/" in outputFile:
        outputFile = outputFile[outputFile.rfind("/"):]
    output = open(outputDir + "/" + outputFile, 'w')
    output.write("*****************RUN STARTS HERE!*****************")
    #start time
    startTime = time.clock()
    output.write("\n" + "Filename: " + FASTAFile + "\n")
    output.write("Program Start: {:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now()) + "\n")
    output.write("Sample Size: " + str(sampleSize) + "\nThreshold: " + str(threshold) + "\n\n")
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
    tree = literal_eval(tree)    #newick format

    # RLR tree required for maxParsimony function
    tree = NewicktoRLR(tree)
    score = maxParsimony(tree, tipMapping)
    
    #Graph required for feasibility check
    graph = nx.Graph()
    makeGraph(graph, tree)
    leaves = getLeaves(tree)
    currentFeasible = isFeasible(graph,leaves)
        
    # Perform NNI heuristic
    counter = 0
    loopCounter = 0
    while True:
        loopCounter += 1
        output.write("Loop Iteration: " + str(loopCounter) + "\n")
        output.write("Loop Start Time: {:%H:%M:%S}".format(datetime.datetime.now()) + "\n")
        output.write("Current Tree\nFeasibility: " + str(currentFeasible) + "\nScore: " + str(score) + "\nTree:\n" + str(tree) + "\n\n")
        NNIs = allNNIs(tree)
        if len(NNIs)-1 < sampleSize:
            sampleSize = len(NNIs)-1
        toScore = random.sample(NNIs, sampleSize)
        
        # add feasibility test
        feasible = []
        infeasible = []
        for tree in toScore:
            graph = nx.Graph()
            makeGraph(graph, tree)
            leaves = getLeaves(tree)
            if isFeasible(graph, leaves): #if this tree is possible
                feasible.append(tree)
            else:
                infeasible.append(tree) #if this tree is not possible
                
        output.write("Number of Feasible Neighbor Trees: " + str(len(feasible)) + "\n")
        output.write("Number of Infeasible Neighbor Trees: " + str(len(infeasible)) + "\n")
        if len(feasible) != 0: #if feasible NNIs were found
            scoredList = map(lambda x: (maxParsimony(x, tipMapping), x), feasible)
            sortedList = sorted(scoredList)
            counter = 0
            if not currentFeasible or sortedList[0][0] < score:
                score = sortedList[0][0]
                tree = sortedList[0][1]
                currentFeasible = True
                output.write("Found a New Feasible Tree!\n\n")
            else:
                output.write("Best Possible Feasible Tree Found\n" + str(tree) + "\n" + "Score: " + str(score) + "\n\n")
                break
        else: #if no possible trees we're found
            if currentFeasible: #checks if the original tree was feasible
                output.write("No Feasible Neighbors, Best Possible Feasible Tree\n" + str(tree) + "\n\n")
                break
            counter += 1
            output.write("Threshold counter: " + str(counter) + "\n\n")
            if counter >= threshold:
                output.write("Threshold Met: No Feasible Tree Found\n")
                stopTime = (time.clock() - startTime)
                output.write("Program Stop: " + str(stopTime) + " seconds\n\n")
                return
            output.write("Searching Infeasible Space\n")
            scoredList = map(lambda x: (maxParsimony(x, tipMapping), x), infeasible)
            sortedList = sorted(scoredList)
            choseNeighbor = False    
            for neighbor in sortedList: #if the original tree was infeasible and no feasible neighbors were found, take the next best infeasible tree and run again
                if neighbor[0] > score:
                    score = neighbor[0]
                    tree = neighbor[1]
                    choseNeighbor = True
                    break
            if not choseNeighbor: 
                score = sortedList[-1][0]
                tree = sortedList[-1][1]
            currentFeasible = False
            output.write("Next Best Infeasible Tree\n\n")
    endTime = (time.clock() - startTime)
    output.write("Program End: " + str(endTime) + " seconds\n\n")
                
    #outputTree = RLRtoNewick(tree)
    #print "Final score", score
    return
    
	
testList1 = ["a_1", "a_1_1", "a_2", "a_2_1", "b_1", "b_1_1"]
testList2 = ["a_1", "a_1_1", "b_2", "b_2_1", "b_1", "b_1_1"]
testTree = (1, ("a", (), ()), (2, ("b_1", (), ()), ("b_2", (), ())))