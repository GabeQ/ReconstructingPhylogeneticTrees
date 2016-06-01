import os
from multiprocessing import Pool
import NNIHeuristic as NNI

def runAll(dataDir, sampleSize, threshold, outputDir):
    p = Pool(5)
    args = []
    midDirs = next(os.walk(dataDir))[1] # sim-flies
    for midDir in midDirs:
        bottomDirs = next(os.walk(dataDir + "/" + midDir))[1] # ___-1x
        for bottomDir in bottomDirs:
            files = next(os.walk(dataDir + "/" + midDir + "/" + bottomDir))[2] # numbered dir
            for f in files:
                if f.endswith('.align'):
                    dataPath = dataDir + "/" + midDir + "/" + bottomDir + "/" + f
                    outputPath = outputDir + "/" + midDir + "/" + bottomDir + "/"
                    os.makedirs(outputPath, exist_ok=True)
                    args.append((dataPath, sampleSize, threshold, outputPath))
    p.map(NNI.NNIHeuristic, args)
                    