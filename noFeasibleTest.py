import NNIHeuristic as NNI
import random
from Bio import Phylo
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import cStringIO
import re
import time
import datetime
from ast import literal_eval

def noFeasibleTest(FASTAFile, sampleSize, outputDir):
    """"takes a FASTAFile, constructs a UPGMA Tree from the file data, converts this tree to RLR format,
    tries to find the tree with the lowest parsimony score (ignores feasibility check)"""
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
    output.write("Sample Size: " + str(sampleSize) + "\n\n")
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
    tree = NNI.NewicktoRLR(tree)
    score = NNI.maxParsimony(tree, tipMapping)
        
    # Perform NNI heuristic
    loopCounter = 0
    while True:
        loopCounter += 1
        output.write("Loop Iteration: " + str(loopCounter) + "\n")
        output.write("Loop Start Time: {:%H:%M:%S}".format(datetime.datetime.now()) + "\n")
        output.write("Current Tree\nScore: " + str(score) + "\nTree:\n" + str(tree) + "\n\n")
        NNIs = NNI.allNNIs(tree)
        if len(NNIs)-1 < sampleSize:
            sampleSize = len(NNIs)-1
        toScore = random.sample(NNIs, sampleSize)
        
        scoredList = map(lambda x: (NNI.maxParsimony(x, tipMapping), x), toScore)
        sortedlist = sorted(scoredList)
        if sortedlist[0][0] < score:
            score = sortedlist[0][0]
            tree = sortedlist[0][1]
            output.write("Found A More Parsimonious Tree!\n\n")
            
        else:
            break
            output.write("No Neighbors With Better Scores Found\n\n")
    output.write("Final Tree:\n" + str(tree) + "\nScore: " + str(score) + "\n\n")
    endTime = (time.clock() - startTime)
    output.write("Program End: " + str(endTime) + " seconds\n\n")
    return
    

	
testList1 = ["a_1", "a_1_1", "a_2", "a_2_1", "b_1", "b_1_1"]
testList2 = ["a_1", "a_1_1", "b_2", "b_2_1", "b_1", "b_1_1"]
testTree = (1, ("a", (), ()), (2, ("b_1", (), ()), ("b_2", (), ())))