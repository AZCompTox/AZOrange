import sys,orange
from AZutilities import  getUnbiasedAccuracy
from trainingMethods import AZorngRF
from AZutilities import dataUtilities
from AZutilities import getCinfonyDesc


def getDesc(trainDataFile):

    # Read the SAR for which to calculate descriptors
    sarData = dataUtilities.DataTable(trainDataFile) 

    # Get names of RDK descriptors
    rdkDescs = getCinfonyDesc.getAvailableDescs("rdk")
   
    # Calculate the descriptors 
    trainData = getCinfonyDesc.getCinfonyDescResults(sarData, rdkDescs)

    # Deselect the SMILES attribute
    attrList = [attr.name for attr in trainData.domain.attributes if attr.varType == orange.Variable.String]
    
    trainData = dataUtilities.attributeDeselectionData(trainData, attrList)
     
    # Save the trainData set
    trainData.save("trainData.tab") 

    return trainData
     


def getAccParamOpt(trainData, AZOrangeLearner, paramList):

    # Parameter optimization and accuracy assessment
    evaluator = getUnbiasedAccuracy.UnbiasedAccuracyGetter(data = trainData, learner = AZOrangeLearner, paramList = paramList, nExtFolds = 5, nInnerFolds = 5)

    # Calculate a results object
    result = evaluator.getAcc()

    # Print the accuracy
    print "Classification Accuracy"
    print result["CA"]
    print "Confusion matrix"
    print result["CM"]


if __name__ == "__main__":
    """
    Script to assess the expected generalization accuracy for a model with optimized parameters. 
    Usage;
    python getAccParamOpt.py sarFilePath
    """

    # Full path to data file with a structure activity relationship (SAR) on which to train, optimize and validate the model
    # Format; SMILES\t ResponseVariable(Any name. The response is recognized by being in the last column)
    # There may be any number of variables inbetween the SMILES and the ResponseVariable columns, which will be used in training togeather with any new descriptors
    #trainDataFile = "/home/jonna/projects/M-Lab/scfbmPaper/data/639sarSmall.txt"
    trainDataFile = sys.argv[1]

    # Calculate descriptors
    trainData = getDesc(trainDataFile)
    
    # Define which AZOrange learner to use by instantiating the learner object
    AZOrangeLearner = AZorngRF.RFLearner()

    # Define the corresponding model parameters to optimize (all possible parameters defined in $AZORANGEHOME/azorange/AZLearnersParamsConfig.py)
    paramList = ["nActVars"]

    # Assess the generalization accuracy of the learner with optimized model parameters in a double loop cross validation
    getAccParamOpt(trainData, AZOrangeLearner, paramList)
