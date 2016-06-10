#compareTreeScript.py

#Created by Gabriel Quiroz
#June 2016

#compareTreeScript is a script that takes in two directories, one containing the
#output files from noFeasibleTest (where the files are ones that are known to
#begin with infeasible trees). The second directory contains output files from
#NNIHeuristic with the same files in the same order. It writes into one file
#the final tree output without the feasiblity test, whether or not that tree is
#feasible, and the score of that tree. It does the same for trees that ran with 
#the feasibility test on. It then computes the percent difference of the scores.

import NNIHeuristic as NNI
from os import listdir
from ast import literal_eval
import networkx as nx

#gets the files from two directories (directories must contain only files)
def getFiles(DirNoFeas, DirFeas):
    Dir1Files = [f for f in listdir(DirNoFeas)]
    Dir2Files = [i for i in listdir(DirFeas)]
    return Dir1Files, Dir2Files


def getTreeAndScore(treeNoFeas, treeFeas):
    f = open(treeNoFeas, 'r')
    content1 = f.readlines()
    for i in range(len(content1)):
        if "Final Tree" in content1[i]:
            tree1 = literal_eval(content1[i+1])
            scoreNoFeas = content1[i+2]
            score1 = scoreNoFeas.split()
            tup1 = (tree1, int(score1[1]))
    i = open(treeFeas, 'r')
    content2 = i.readlines()
    for j in range(len(content2)):
        if "Best Possible Feasible Tree Found" in content2[j]:
            tree2 = literal_eval(content2[j+1])
            scoreFeas = content2[j+2]
            score2 = scoreFeas.split()
            tup2 = (tree2, int(score2[1]))
    return tup1, tup2
    
def getFeasibility(tree):
    G = nx.Graph()
    leaves = NNI.getLeaves(tree)
    NNI.makeGraph(G, tree)
    return NNI.isFeasible(G, leaves)
    

def compareTreesScript(Dir1NoFeas, Dir2Feas, outputDir):
    output = open(outputDir + '/compareTrees.out', 'w')
    output.write('Compare Trees\n\n')
    files1, files2 = getFiles(Dir1NoFeas, Dir2Feas)
    counter = 1
    for i in range(len(files1)):
        output.write('********** Comparison ' + str(counter) + ': ' + str(files1[i]) + ' **********\n\n')
        tup1, tup2 = getTreeAndScore(Dir1NoFeas +'/' + files1[i], Dir2Feas + '/' + files2[i])
        output.write('Tree Without Feasibility Test\n' + str(tup1[0]) + '\nIs Tree Feasible: ' + str(getFeasibility(tup1[0])) + '\nScore: ' + str(tup1[1]) + '\n\n')
        output.write('Tree With Feasibility Test\n' + str(tup2[0]) + '\nIs Tree Feasible: ' + str(getFeasibility(tup2[0])) + '\nScore ' + str(tup2[1]) + '\n\n')
        difference = tup1[1] - tup2[1]
        percentage = difference / float(tup2[1]) * 100
        output.write('Difference Percentage: ' + str(percentage) + '%\n\n')
        counter += 1
    output.write('Finished Running Data')
    return
        
        
        
        
        
        
        
        