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
    
def getNNIArgs(dataDir, sampleSize, outputDir):
    args = []
    subDirs = next(os.walk(dataDir))[1] # ___-1x folders
    for subDir in subDirs:
        files = next(os.walk(dataDir + "/" + subDir))[2] # numbered dir
        for f in files:
            if f.endswith('.align'):
                dataPath = dataDir + "/" + subDir + "/" + f
                outputPath = outputDir + "/" + subDir + "/"
                if not os.path.exists(os.path.dirname(outputPath)):
                    try:
                        os.makedirs(os.path.dirname(outputPath))
                    except OSError as exc: # Guard against race condition
                        if exc.errno != errno.EEXIST:
                            raise
                args.append((dataPath, sampleSize, outputPath))
    return args
    
def callNNI(args):
    NNI.noFeasibleTest(args[0], args[1], args[2])

if __name__ == '__main__':
    dataDir = sys.argv[1]
    sampleSize = int(sys.argv[2])
    outputDir = sys.argv[3]
    num_processors = int(sys.argv[4])
    p = Pool(num_processors)
    p.map(callNNI, getNNIArgs(dataDir, sampleSize, outputDir))