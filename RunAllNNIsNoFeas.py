import os
import sys
from multiprocessing import Pool
import errno
import noFeasibleTest as NNI

#def getNNIArgs(dataDir, sampleSize, threshold, outputDir):
#    args = []
#    midDirs = next(os.walk(dataDir))[1] # sim-flies
#    for midDir in midDirs:
#        bottomDirs = next(os.walk(dataDir + "/" + midDir))[1] # ___-1x folders
#        for bottomDir in bottomDirs:
#            files = next(os.walk(dataDir + "/" + midDir + "/" + bottomDir))[2] # numbered dir
#            for f in files:
#                if f.endswith('.align'):
#                    dataPath = dataDir + "/" + midDir + "/" + bottomDir + "/" + f
#                    outputPath = outputDir + "/" + midDir + "/" + bottomDir + "/"
#                    if not os.path.exists(os.path.dirname(outputPath)):
#                        try:
#                            os.makedirs(os.path.dirname(outputPath))
#                        except OSError as exc: # Guard against race condition
#                            if exc.errno != errno.EEXIST:
#                                raise
#                    args.append((dataPath, sampleSize, threshold, outputPath))
#    return args
    
"""searching through the output files of the feasibility test
        collecting the names of the files that had infeasible trees
        finding the files with the same name in the dataDir,
        running NNIHeuristic on those
        returning the outputs in outputDir"""
def getNNIArgs(feasOutputDir, dataDir, sampleSize, outputDir):
    args = []
    for folderName, subfolders, filenames in os.walk(feasOutputDir): 
        for filename in filenames:
            stringName = str(filename)
            while True:
                print('File inside: ' + folderName + ' <- of that folder: '+ filename)
                os.chdir(folderName)
                if 'Feasibility: False' in open(stringName).read():
                    print "This file contains infeasibility: " + stringName
                    dataFile = stringName.replace(".out", ".align")
                    os.chdir(dataDir)
                    subDirs = next(os.walk(dataDir))[1] # ___-1x folders
                    for subDir in subDirs:
                        files = next(os.walk(dataDir + "/" + subDir))[2] # numbered dir
                        for f in files:
                            if f == dataFile:
                                dataPath = dataDir + "/" + subDir + "/" + f
                                outputPath = outputDir + "/" + subDir + "/"
                                if not os.path.exists(os.path.dirname(outputPath)):
                                    try:
                                        os.makedirs(os.path.dirname(outputPath))
                                    except OSError as exc: # Guard against race condition
                                        if exc.errno != errno.EEXIST:
                                            raise
                                args.append((dataPath, sampleSize, outputPath))
                break
    return args
    
def callNNI(args):
    NNI.noFeasibleTest(args[0], args[1], args[2])

if __name__ == '__main__':
    feasOutputDir = sys.argv[1]
    dataDir = sys.argv[2]
    sampleSize = int(sys.argv[3])
    outputDir = sys.argv[4]
    num_processors = int(sys.argv[5])
    p = Pool(num_processors)
    p.map(callNNI, getNNIArgs(feasOutputDir, dataDir, sampleSize, outputDir))