import os
import NNIHeuristic as NNI

def runAll(dataDir, outputDir):
    filePath = dataDir
    outputPath = outputDir
    midDirs = next(os.walk(dataDir))[1] # sim-flies
    for midDir in midDirs:
        bottomDirs = next(os.walk(dataDir + "/" + midDir))[1] # ___-1x
        for bottomDir in bottomDirs:
            files = next(os.walk(dataDir + "/" + midDir + "/" + bottomDir))[2] # numbered dir
            for f in files:
                if f.endswith('.align'):
                    NNI.NNIHeuristic(dataDir + "/" + midDir + "/" + bottomDir + "/" + f, "")
                    