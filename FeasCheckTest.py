# Script to test Prof Wu's feasibility checker against ours
# Matt Dohlen, June 2016

import re
import cStringIO
import os
from ast import literal_eval
import plct.rasmus.treelib as treelib
import plct.plctlib as plctlib
from NNIHeuristic import *
from Bio import Phylo
import networkx as nx

def subDirPaths(d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])


outFile = open("FeasCheckTest.out", 'w')
outSum = open("FeasCheckTestSummary.out", 'w')
buf = cStringIO.StringIO()
directory = "/Users/cssummer16/Source/sim-flies/100e6-1x"
wuInfeasCount, pwnInfeasCount, difCount, totCount = 0, 0, 0, 0
for sub in subDirPaths(directory):
    for f in os.listdir(sub):
        if f.endswith(".tree"):
            # Read in Prof Wu's tree
            wuTree = treelib.Tree()
            wuTree.read_newick(sub + "/" + f)

            # Read in North Pawn's tree
            pwnTree = Phylo.read(sub + "/" + f, 'newick')
            i = 0
            for clade in pwnTree.get_nonterminals():
                clade.name = "Inner" + str(i)
                i += 1
            for clade in pwnTree.get_terminals():
                clade.name = "\"" + clade.name + "\""
            buf = cStringIO.StringIO()
            Phylo.write(pwnTree, buf, 'newick', plain=True)
            pwnTree = buf.getvalue()
            pwnTree = re.sub(r'Inner\d*', '', pwnTree)
            pwnTree = pwnTree.replace(";", "")
            pwnTree = literal_eval(pwnTree)
            pwnTree = NewicktoRLR(pwnTree)

            # Check Feasibility of trees
            wuFeas = plctlib.is_reconcilable(wuTree, mapping='sli_')
            pwnLeaves = getLeaves(pwnTree)
            pwnGraph = nx.Graph()
            makeGraph(pwnGraph, pwnTree)
            pwnFeas = isFeasible(pwnGraph, pwnLeaves)
            if not wuFeas:
                wuInfeasCount += 1
            if not pwnFeas:
                pwnInfeasCount += 1


            # Output Findings
            outFile.write("File Name: " + f + "\n")
            if wuFeas != pwnFeas:
                difCount += 1
                outFile.write("=" * 54 + "\n")
                outFile.write("PROF WU'S FINDINGS WERE DIFFERENT FROM THE NORTH PAWNS\n")
                outFile.write("=" * 54 + "\n")
            outFile.write("Prof Wu's findings: " + str(wuFeas) + "\n")
            outFile.write("North Pawn's findings: " + str(pwnFeas) + "\n\n")
            outFile.write(str(pwnTree) + "\n\n")
            totCount += 1
            buf.flush()

# Output Summary
outSum.write("Summary of feasibility check test\n\n")
outSum.write("# Files: " + str(totCount) + "\n")
outSum.write("Prof Wu # infeasible: " + str(wuInfeasCount) + "\n")
outSum.write("Prof Wu % infeasible: {0:.2f}\n".format((float(wuInfeasCount)/totCount) * 100))
outSum.write("North Pawn # infeasible: " + str(pwnInfeasCount) + "\n")
outSum.write("North Pawn % infeasible: {0:.2f}\n".format((float(pwnInfeasCount)/totCount) * 100))
outSum.write("# different: " + str(difCount))
