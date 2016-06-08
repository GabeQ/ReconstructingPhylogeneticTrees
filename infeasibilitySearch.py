# Searches outfiles from NNIHeuristic to look for cases where infeasible trees were found.
# Chen Pekker, June 2016

import shutil
import os


for folderName, subfolders, filenames in os.walk('/Users/cssummer16/Documents/100e6-1x/'): #need to change this for every directory
    for filename in filenames:
        stringName = str(filename)
        while True:
            print('File inside: ' + folderName + ' <- of that folder: '+ filename)
            os.chdir(folderName)
            print "HERE"
            if 'Feasibility: False' in open(stringName).read():
                print "File contains infeasibility: " + stringName
                shutil.copy(filename, "/Users/cssummer16/GitHub/ReconstructingPhylogeneticTrees/treesWithInfeasibilities") #need to change the output file
            print "after while loop before break" + folderName
            break
                
    print('')

