import os
from Bio import AlignIO

def getKey(item):
    return item.id
    
def getGeneAndLocus(fileName):
    if fileName.endswith('.align'):
        myAlignment = AlignIO.read(fileName, "fasta")
        geneLocusList = []
        sorted(myAlignment, key=getKey)
        for align in myAlignment:
            geneLocusSpec = align.id.split('_')
            geneLocus = (geneLocusSpec[0], geneLocusSpec[1])
            if geneLocus not in geneLocusList:
                geneLocusList.append(geneLocus)
        return geneLocusList
    return []
    

def possiblyInfeasible(fileName):
    geneLocusList = getGeneAndLocus(fileName)
    infeasibleGenes = set([])
    for i in range(len(geneLocusList)):
        for j in range(i+1, len(geneLocusList)):
            if geneLocusList[i][0] == geneLocusList[j][0] and geneLocusList[i][1] != geneLocusList[j][1]:
                infeasibleGenes.union(geneLocusList[i][0])
                break
    return infeasibleGenes == set([])
        

def getPossiblyInfeasible(dataDirectory):
    potentialFiles = []
    for dirName, subDirList, fileList in os.walk(dataDirectory):
        if len(fileList) != 0:
            for subDir in subDirList:
                potentialFiles + getPossiblyInfeasible(subDir)
        else:
            for f in fileList:
                if possiblyInfeasible(f):
                    potentialFiles.append(f)
    return potentialFiles