import os
import sys
import time
import orngTest
import orngStat
from trainingMethods import AZorngRF
from trainingMethods import AZBaseClasses
from AZutilities import dataUtilities
from AZutilities import paramOptUtilities



if __name__ == "__main__":
    """
    Script to calculate the accuracy on a temporal test set with a saved model. 
    Usage;
    python getTempAcc.py testDataPath modelPath
    """

    # Full path to temporal test data file (with descriptors) in Orange format
    #testDataFile = "/home/jonna/projects/M-Lab/scfbmPaper/data/trainData.tab"
    testDataFile = sys.argv[1]

    # Read the test data
    testData = dataUtilities.DataTable(testDataFile)

    # Full path to the model 
    #modelFile = "/home/jonna/projects/M-Lab/scfbmPaper/data/optRF.model" 
    modelFile = sys.argv[2]

    # Read the model
    model = AZBaseClasses.modelRead(modelFile)

    # Use Orange methods to get the accuracy (Please see Orange doc)
    results = orngTest.testOnData([model], testData) 

    print "Classification accuracy"
    print orngStat.CA(results)
    
    

