import sys
import os
import time
from trainingMethods import AZorngRF
from AZutilities import dataUtilities
from AZutilities import paramOptUtilities



if __name__ == "__main__":
    """
    Algorithm to optimize the model parameters and use the optimized parameters to train a model that is saved to disk. 
    Usage;
    python buildOptParamModel.py trainDataPath optDirPath modelPath
    """

    # Full path to train data file (with descriptors) in Orange format
    #trainDataFile = "/home/jonna/projects/M-Lab/scfbmPaper/data/trainData.tab"
    trainDataFile = sys.argv[1]

    # Define which AZOrange learner to use by instantiating the learner object
    AZOrangeLearner = AZorngRF.RFLearner()

    # Without access to a distributed computational environment set to 'NoSGE'
    queueType = 'NoSGE'

    # Directory (will be created) in which to write the results files from the parameter optimization. The file optimazationLog.txt summarizes these results. 
    #optDirRoot = "/home/jonna/projects/M-Lab/scfbmPaper/data/paramOpt"
    optDirRoot = sys.argv[2]
    runPath = optDirRoot+str(time.time())
    os.system("mkdir "+runPath)

    # Get a learner object with optimized parameters (default settings)
    print "Optimizing model hyper-parameters"
    optLearner, isOptimized = paramOptUtilities.getOptParam(AZOrangeLearner, trainDataFile, verbose = 0, queueType = queueType, runPath = runPath)

    print "Parameters successfully optimized?"
    print isOptimized

    # Load the data on which to train the model
    trainData = dataUtilities.DataTable(trainDataFile)

    # Build the model with optimized parameters
    print "Building model with optimized parameters"
    model = optLearner(trainData)

    # Save the model
    #modelPath = "/home/jonna/projects/M-Lab/scfbmPaper/data/optRF.model"
    modelPath = sys.argv[3]
    model.write(modelPath)
    

