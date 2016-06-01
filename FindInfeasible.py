import os
from Bio import AlignIO

def getKey(item):
    return item.id
    
def getGeneAndLocus(fileName):
    myAlignment = AlignIO.read(fileName, "fasta")
    geneLocusList = []
    myAlignment = sorted(myAlignment, key=getKey)
    for align in myAlignment:
        geneLocusSpec = align.id.split('_')
        print geneLocusSpec
        geneLocus = (geneLocusSpec[0], geneLocusSpec[1])
        if geneLocus not in geneLocusList:
            geneLocusList.append(geneLocus)
    return geneLocusList
    

def possiblyInfeasible(fileName):
    geneLocusList = getGeneAndLocus(fileName)
    infeasibleGenes = set([])
    for i in range(len(geneLocusList)):
        for j in range(i+1, len(geneLocusList)):
            if geneLocusList[i][0] == geneLocusList[j][0] and geneLocusList[i][1] != geneLocusList[j][1]:
                infeasibleGenes.union(geneLocusList[i][0])
                break
    if infeasibleGenes != set([]):
        return fileName
        
        

def getPossiblyInfeasible(dataDirectory):
    potentialFiles = []
    filePath = dataDirectory
    for dirName, subDirList, fileList in os.walk(dataDirectory):
        for subDir in subDirList:
            subFilePath = filePath + "/" + subDir
            for dirName, subDirList, fileList in os.walk(dataDirectory):
                for subDir in subDirList:
                    subSubFilePath = subFilePath + "/" + subDir
                    for f in fileList:
                        if f.endswith('.align') and possiblyInfeasible(subSubFilePath + "/" + f):
                            potentialFiles.append(f)
    return potentialFiles