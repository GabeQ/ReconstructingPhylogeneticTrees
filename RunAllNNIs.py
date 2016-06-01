import os

def runAll(dataDir, outputDir):
    filePath = dataDir
    outputPath = outputDir
    dirList = next(os.walk(dataDir))[1] # sim-flies
    print dirList
    #for subDir in topDir[1]:
    #    midDir = os.walk(subDir) # ___-1x
    #    for midSubDir in midDir[1]:
    #        bottomDir = os.walk(midSubDir) # numbered dir
    #        for f in bottomDir[2]:
    #            if f.endswith('.align'):
                    