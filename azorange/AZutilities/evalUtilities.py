import types, os
import time
import math,struct
import statc
import string
import orngStat
import copy
import orngTest
import orange
import numpy
import commands
from statlib import stats
from AZutilities import miscUtilities
from rdkit import Chem
from rdkit.Chem import Draw
version = 2
verbose = 0

CLASSIFICATION = 1  # 0b00000001
REGRESSION = 2      # 0b00000010

def tanimotoSimilarity(A, B):
    """Calculates the tanimnoto similarity defined by:
             ts = c/(a+b-c)
        where a is the number of bits in A
              b is the number of bits in B
              c is the number of bits in common
             
        A and B are the strings containing the bits
        This returns a value between 0 and 1:
                0 - Nothing in common
                1 - Identical
        len(A) and len(A) and len(B) are expected to be the same
    """
    #Not validating input since this function has to called very often.
    nAB = 0.0 + miscUtilities.countOnBits(int(A,2) & int(B,2))
    nAnB = A.count("1") + B.count("1")
    if not nAnB: return 0.0
    else: return  nAB / (nAnB - nAB)

def fastTanimotoSimilarity(A,nA,B):
    """Same as tanimotoSimilarity but expecting in A and B a long integer representing the binary fingetprint
       and nA the number of ONbits in A"""
    #Not validating input since this function has to called very often.
    AandB = A & B
    nAB = miscUtilities.countOnBits(AandB)      # Number of match ON Bits
    nAnB = 0.0 + nA + miscUtilities.countOnBits(B)            # nA+nB
    if not nAnB: return 0.0
    else: return  nAB / (nAnB - nAB)



def getNearestNeighbors(query, n, NNDataPath, FPPath = None, resPath = None, idx = 0):
    """ get the n nearest neighbors
        query: bin string with query fingerprint
        returns an ordered list with the n top neighbors (each one in a dict):
            [ {
                "id"          : ID, 
                "expVal"      : ExpValues, 
                "similarity"  : TanimotoSimilarity, 
                "smi"         : smiles, 
                "imgPath"     : imgPath},  ... ]        

        It will saves the images in resPath:
             NN_1.png    #1 neighbor
             NN_2.png    #2 neighbor
             ...
             NN_n.png    #n neighbor
    """
    if not query or not n or not  NNDataPath or not  FPPath:
        return []
    #if resPath and not os.path.isdir(resPath):
    #    os.makedirs(resPath)

    # get the correct header
    file = open(NNDataPath,"r")
    header = file.readline().strip().split('\t')
    file.close()

    if "Molecule SMILES" not in header or "Compound Name" not in header:
        print "NN dataset ",NNDataPath, " have not the correct header. It must contain 'Molecule SMILES' and 'Compound Name' attributes."
        return [] 
    # Index will have to be sum 1 because the TS will be prepended
    idxID = header.index("Compound Name") + 1
    idxExpVal = len(header) 
    idxSMILES = header.index("Molecule SMILES") + 1
    idxSimilarity = 0


    Nbits = 2048
    status,output = commands.getstatusoutput('echo "' + query + '" | fpin ' + FPPath + " "  +NNDataPath + ' 0.0 '+str(n))
    if status:
        print status
        print output
        raise Exception(str(output))
    #             TS              SMILES                    AZID         DATE       expRes
    # output = "0.7117   CCCC(C)C1(C(=O)NC(=O)NC1=O)CC   AZ10046012   2009-12-02   3.480007"
    TS=[]
    for ts in output.split("\n"):
        TS.append(ts.strip().split('\t'))
    # in TS:
    #    TS[n][0] - tanimoto similarity
    #    TS[n][1] - SMILES
    #    TS[n][2] - AZID
    #    TS[n][-1]- expRes
    res = []
    timeStamp=str(time.time()).replace(".",'')
    for fidx,nn in enumerate(TS):
        ID= nn[idxID]
        expVal = nn[idxExpVal]
        SMILES = nn[idxSMILES]
        if resPath and os.path.isdir(resPath):
            imgPath = os.path.join(resPath,"NN"+str(idx)+"_"+str(fidx+1)+"_"+timeStamp+".png")
            mol = Chem.MolFromSmiles(SMILES)
            # save the respective imgPath...  
            Draw.MolToImageFile(mol,imgPath,size=(300, 300), kekulize=True, wedgeBonds=True)
        else:
            imgPath = ""
        res.append( {
                "id": ID, 
                "expVal": expVal, 
                "similarity": nn[idxSimilarity], 
                "smi": SMILES, 
                "imgPath": imgPath} )
    return res
    


def ConfMat(res = None):
        """ Returns a confusion matrix in the form of a vector:
                For Binary classifiers:  
                        [[TP, FN],
                         [FP, TN]]
                For classifiers with class having N values:
                                             Predicted class
                                        |   A       B       C
                                     ---------------------------
                           known     A  |  tpA     eAB     eAC
                           class     B  |  eBA     tpB     eBC
                                     C  |  eCA     eCB     tpC

                    [[tpA, eAB, ..., eAN],
                     [eBA, tpB, ..., eBN],
                      ...,
                     [eNA, eNB, ..., tpN ]]

                 where A, B, C are the class values in the same order as testData.domain.classVar.values
        """
        if res == None:
            return {"type":CLASSIFICATION}

        confMat = orngStat.confusionMatrices(res)[0]
        if len(res.classValues) == 2:
            cm = [[confMat.TP, confMat.FN],[confMat.FP, confMat.TN]]
        else:
            cm = confMat
        return cm

def getConfMat(testData, model):
        return ConfMat(orngTest.testOnData([model], testData))

def calcKappa(_CM):
    """Returns the Kappa statistical coefficient for the agreement between measured and predicted classes"""
    def calcClassAccuracy(CMx):
        """Returns the proportion correctly classified"""
        correctClass = 0
        for i in range(len(CMx)) : correctClass = correctClass + CMx[i][i]
        total = sum(sum(CMx))
        return correctClass / float(total)
    CM = numpy.array(_CM)
    p_correctlyClassed = calcClassAccuracy(CM)
    measured = sum(CM)
    p_measured = numpy.array(map(lambda x : x/sum(measured), measured))
    predicted = sum(CM.T)
    p_predicted = numpy.array(map(lambda x : x/sum(predicted), predicted))
    prob_chance = sum(p_measured * p_predicted)
    return (p_correctlyClassed - prob_chance ) / (1 - prob_chance)

def Kappa(res=None):
    if res == None:
        return {"type":CLASSIFICATION}
    return [calcKappa(ConfMat(res))]

def generalCVconfMat(data, learners, nFolds = 5):
    """
    General method for printing the X fold CV confusion matrix of an Orange data set (data)
    with any number of classes. learners is a list of AZorange learners.
    """

    res = orngTest.crossValidation(learners, data, strat=orange.MakeRandomIndices.StratifiedIfPossible, folds = nFolds)
    classes = data.domain.classVar.values

    for idx in range(len(learners)):
        cm = orngStat.computeConfusionMatrices(res)[idx]
        print "Results for "+learners[idx].name
        print "\t"+"\t".join(classes)
        for className, classConfusions in zip(classes, cm):
            print ("%s" + ("\t%i" * len(classes))) % ((className, ) + tuple(classConfusions))

def getClassificationAccuracy(testData, classifier):
    correct = 0.0
    for ex in testData:
        #print str(classifier(ex)) + "->" + str(ex.getclass())
        if str(classifier(ex)) == str(ex.getclass()):
            correct = correct + 1.0
    ClassificationAccuracy = correct/len(testData)
    return ClassificationAccuracy


def getRMSE(testData, predictor):

    accSum = 0.0
    nPredEx = 0
    for ex in testData:
        #print "Class "+str(ex.getclass())+" "+str(string.atof(str(predictor(ex))))
        try:
            #print "(",ex.getclass(),"-",string.atof(str(predictor(ex))),")^2= ",(ex.getclass()-string.atof(str(predictor(ex))))**2,"(",math.pow(ex.getclass()-string.atof(str(predictor(ex))),2),")" 
            accSum = accSum + math.pow(ex.getclass()-string.atof(str(predictor(ex))),2)
            nPredEx = nPredEx + 1 
        except:
            if verbose > 0: print "Warning!!!!"
            if verbose > 0: print "No prediction could be made for the following example:"
            if verbose > 0: print ex
    
    accuracy = math.sqrt(accSum/nPredEx)
    return accuracy


def getRsqrt(testData, predictor):
    """Calculate the coefficient of determination (R-squared) for the orange model predictor on the data set testData. 
        This uses the Test Set Activity Mean
        R^2 = 1 - sum((pred - actual)^2)/(sum((testMean - actual)^2))"""

    # Calc average of the response variable
    actualValuesList = []
    for ex in testData:
        actualValuesList.append(ex.getclass().value)
    testMean = statc.mean(actualValuesList)

    errSum = 0.0
    meanSum = 0.0
    for ex in testData:
        errSum = errSum + math.pow(ex.getclass() - string.atof(str(predictor(ex))),2)
        meanSum = meanSum + math.pow(testMean - ex.getclass(),2)

    Rsqrt = 1 - errSum/meanSum
    return Rsqrt


def getQ2(testData, predictor):
    """Calculate the predictive squared correlation coefficient Q2, which uses the Training Set Activity Mean.. 
        Q^2 = 1 - sum((pred - actual)^2)/(sum((trainMean - actual)^2))
        This is according:   
                Comments on the Definition of the Q2 Parameter for QSAR Validation
                Viviana Consonni, Davide Ballabio, Roberto Todeschini
                Journal of Chemical Information and Modeling 2009 49 (7), 1669-1678"""

    # Test if the predictor is compatible with getQ2
    if not testData.domain.classVar or testData.domain.classVar.varType != orange.VarTypes.Continuous:
        print "The dataset is not suitable for calculatinmg Q2. Data must have Continuous Class" 
        return None
    elif not hasattr(predictor,"basicStat") or not predictor.basicStat:
        print "Q2 Error: The predictor is not compatible with the use fo getQ2. It has no basicStat defined."
        return None
    # Calc average of the training class variable
    trainMean = predictor.basicStat[testData.domain.classVar.name]["avg"]

    errSum = 0.0
    meanSum = 0.0
    for ex in testData:
        errSum = errSum + math.pow(ex.getclass() - string.atof(str(predictor(ex))),2)
        meanSum = meanSum + math.pow(trainMean - ex.getclass(),2)

    Q2 = 1 - errSum/meanSum
    return Q2


def Sensitivity(confMatrixList, classes):
##scPA  Added last line of next comment ##ecPA
    """
    Takes a orngStat.confusionMatrices output object with N classes and computes the sensitivity 
    for each class returned as a dictionary indexed by the class name.
    The return object is a list with the lenght of the number of confusion matrices in confMatrixList. 
    The dictionary holding the sensitivities constitutes the objects of the returned list. 
    Sensitivity is obtained by dividing the diagonal element by the row sum.
    If the row sum is 0, it sets the Sensitivity to "N/A"
    """
   
    sensitivityList = []
    for confMatrix in confMatrixList:
        sensitivityDict = {}
        # Loop over the rows of the confusion matrix
        for idx in range(len(classes)):
##scPA
            if sum(confMatrix[idx])==0:
                sensitivityDict[classes[idx]] = "N/A"
            else:        
##ecPA
                sensitivityDict[classes[idx]] = confMatrix[idx][idx]/sum(confMatrix[idx])
        sensitivityList.append(sensitivityDict)

    #print "End sensitivity "+str(sensitivityList)
    return sensitivityList
        

def Predictivity(confMatrixList, classes):
##scPA  Added last line of next comment ##ecPA
    """
    Takes a orngStat.confusionMatrices output object with N classes and computes the sensitivity 
    for each class returned as a dictionary indexed by the class name.
    The return object is a list with the lenght of the number of confusion matrices in confMatrixList. 
    The dictionary holding the sensitivities constitutes the objects of the returned list. 
    Predictivity is obtained by dividing the diagonal element by the column sum of the confusion matrix.
    If one class was never predicted, the predictivity is set to "N/A"
    """

    PredictivityList = []
    for confMatrix in confMatrixList:
        PredictivityDict = {}
        # Loop over the rows of the confusion matrix
        for idx in range(len(classes)):
            colSum = 0
            # Loop over the columns of confMatrix 
            for innerIdx in range(len(classes)):
                colSum = colSum + confMatrix[innerIdx][idx]
##scPA
            if colSum==0:
                PredictivityDict[classes[idx]] = "N/A"
            else:
##ecPA
               PredictivityDict[classes[idx]] = confMatrix[idx][idx]/colSum
        PredictivityList.append(PredictivityDict)

    #print "End Predictivity "+str(PredictivityList)
    return PredictivityList            


def AUConeVSall(results, values, AUCList):
    """Retrun a list of lists where each element is a dictionary with one element per class with the value AUC_single. 
       Uses the orngTest.ExperimentResults object and data.domain.classVar.values"""
    #AUCList = []
    # Loop over the results for each learner in results
    # This is a work around as AUC_single does not seem to work with multiple learners.
    for keepIdx in range(results.numberOfLearners):
        result = copy.deepcopy(results)
        # Remove all learners except keepIdx
        for idx in range(results.numberOfLearners):
            if idx != keepIdx: result.remove(idx)
        AUCdict = {}
        for value in values:
            # Uses own modified version of this function because of orange bug.
            if result.numberOfIterations > 1:
                AUCdict[value] = AUC_single(result, classIndex = values.index(value))[0]
            else:
                AUCdict[value] = AUC_single(result, classIndex = values.index(value))[0][0]
        AUCList.append([AUCdict])

    return AUCList


# Computes AUC; in multivalued class problem, AUC is computed as one against all
# Results over folds are averages; if some folds examples from one class only, the folds are merged
def AUC_single(res, classIndex = -1, useWeights = True):
    if classIndex<0:
        if res.baseClass>=0:
            classIndex = res.baseClass
        else:
            classIndex = 1

    if res.numberOfIterations > 1:
        return orngStat.AUC_iterations(orngStat.AUC_i, orngStat.splitByIterations(res), (classIndex, useWeights, res, res.numberOfIterations))
    else:
        return orngStat.AUC_i(res, classIndex, useWeights)
        #return AUC_i([res], classIndex, useWeights)

def Q2(res = None, trainMeans = None):
    """Calculate the predictive squared correlation coefficient Q2, which uses the Training Set Activity Mean.. 
        Q^2 = 1 - sum((pred - actual)^2)/(sum((trainMean - actual)^2))
        res will have n Learner predictions and trainMeans the mean of the class on the trainSet resectively to the n Learners
        This is according:   
                Comments on the Definition of the Q2 Parameter for QSAR Validation
                Viviana Consonni, Davide Ballabio, Roberto Todeschini
                Journal of Chemical Information and Modeling 2009 49 (7), 1669-1678"""
    if res == None:
        return {"type":REGRESSION}
    nLearners = len(res.results[0].classes)
    if trainMeans == None or len(trainMeans) != nLearners:
        print "Q2 ERROR: The train mean must be provided in order to calculate Q2"
        return None

    if res.numberOfIterations > 1:
        print "ERROR: Q2 does not support more than one iteration!"
        return None
    else:
        errSums = [0.0]*res.numberOfLearners
        meanSums = [0.0]*res.numberOfLearners
        for tex in res.results:
            errSums = map(lambda res, cls, ac = float(tex.actualClass):
                       res + (float(cls) - ac)**2, errSums, tex.classes)
            meanSums = map(lambda res, mean, ac = float(tex.actualClass):
                       mean and res + (ac - mean)**2 or 0.0, meanSums, trainMeans)

        Q2s =  [ x[1] and 1-(x[0]/x[1]) or None for x in zip(errSums, meanSums)]

    return [x!=None and round(x,3) or None for x in Q2s]


def R2(res = None):
    """
    Truncate the orange method to 3 decimals. Allow for no input arguments. Used by the optimizer.
    """
    if res == None:
        return {"type":REGRESSION}
    else:
        scores = orngStat.R2(res)    
        return [round(x,3) for x in scores]

def RMSE(res = None):
    """
    Truncate the orange method to 3 decimals. Allow for no input arguments. Used by the optimizer.
    """
    if res == None:
        return {"type":REGRESSION}
    else:
        scores = orngStat.RMSE(res)    
        return [round(x,3) for x in scores]

def CA(res = None):
    """
    Truncate the orange method to 3 decimals. Allow for no input arguments. Used by the optimizer.
    """
    if res == None:
        return {"type":CLASSIFICATION}
    else:
        scores = orngStat.CA(res)    
        return [round(x,3) for x in scores]



##scPA
def Rsqrt_obsolete(res = None):
    """
    Calculates the R-squared (Coefficient of determination) of orngTest.ExperimentResults in res
    The results res must be from a learner
    """
    # If Called without arguments, return the type of problems this method can be used for: 
    # 1 - Classification problems (Discrete Class)
    # 2 - Regression problems (Continuous Class)
    # 3 - Both Regression and Classification problems (Continuous or Discrete Class)
    if res == None:
        return {"type":REGRESSION}

    if res.numberOfIterations > 1:
        Rs = [[0.0] * res.numberOfIterations for i in range(res.numberOfLearners)]
        errSum = [[0.0] * res.numberOfIterations for i in range(res.numberOfLearners)]
        meanSum = [[0.0] * res.numberOfIterations for i in range(res.numberOfLearners)]
        means = [[0.0] * res.numberOfIterations for i in range(res.numberOfLearners)]
        nIter = [0]*res.numberOfIterations
        for tex in res.results:
            ac = float(tex.actualClass)
            nIter[tex.iterationNumber] += 1
            for i, cls in enumerate(tex.classes):
                means[i][tex.iterationNumber] += ac
        for nit, it in enumerate(nIter):
            for i, cls in enumerate(tex.classes):
                means[i][nit] /=it

        for tex in res.results:
            ac = float(tex.actualClass)
            for i, cls in enumerate(tex.classes):
                errSum[i][tex.iterationNumber] += (float(cls) - ac)**2
                meanSum[i][tex.iterationNumber] += (means[i][tex.iterationNumber] - ac)**2
        for learner in range(res.numberOfLearners):
            for it in range(len(nIter)):
                if meanSum[learner][it]==0:
                    return "N/A"
                Rs[learner][it] = 1-(errSum[learner][it] / meanSum[learner][it])
        return [statc.mean(x) for x in Rs]

    else:
        RsqrtList=[]
        for nLearner in range(len(res.results[0].classes)):
            # Calc average of the prediction variable
            testMean = 0
            for ex in res.results:
                testMean = testMean + ex.actualClass
            testMean = testMean/len(res.results)
            errSum = 0.0
            meanSum = 0.0
            for ex in res.results:
                errSum = errSum + math.pow(ex.actualClass - ex.classes[nLearner],2)
                meanSum = meanSum + math.pow(testMean - ex.actualClass,2)
            if meanSum==0:
                return "N/A"
            RsqrtList.append(1 - errSum/meanSum)
        return RsqrtList

def RMSE_obsolete(res = None):
    """
    Calculates the Root Mean Squared Error of orngTest.ExperimentResults in res
    The results res must be from a regressor
    """
    # If Called without arguments, return the type of problems this method can be used for: 
    # 1 - Classification problems (Discrete Class)
    # 2 - Regression problems (Continuous Class)
    # 3 - Both Regression and Classification problems (Continuous or Discrete Class)
    if res == None:
        return {"type":REGRESSION}

    if res.numberOfIterations > 1:
        MSEs = [[0.0] * res.numberOfIterations for i in range(res.numberOfLearners)]
        nIter = [0]*res.numberOfIterations
        for tex in res.results:
            ac = float(tex.actualClass)
            nIter[tex.iterationNumber] += 1
            for i, cls in enumerate(tex.classes):
                MSEs[i][tex.iterationNumber] += (float(cls) - ac)**2
        MSEs = [[x/ni for x, ni in zip(y, nIter)] for y in MSEs]
        MSEs = [[math.sqrt(x) for x in y] for y in MSEs]

        # Print output from each fold to tem file
        RMSEfoldList = MSEs
        RMSE = [statc.mean(x) for x in RMSEfoldList]
        RMSEstd = stats.stdev(RMSEfoldList[0])
        #print str(RMSE[0])+"\t"+str(RMSEstd)+"\t"+string.join( [str(x) for x in RMSEfoldList[0]] , "\t")

        return [round(statc.mean(x),2) for x in MSEs]

    else:
        MSEs = [0.0]*res.numberOfLearners
        for tex in res.results:
            MSEs = map(lambda res, cls, ac = float(tex.actualClass):
                       res + (float(cls) - ac)**2, MSEs, tex.classes)

        MSEs = [x/(len(res.results)) for x in MSEs]
        return [round(math.sqrt(x),2)  for x in MSEs]


def CA_obsolete(res = None, returnFoldStat = False):
    """
    Calculates the classification Accuracy of orngTest.ExperimentResults in res
    The results res must be from a classifier
    """
    # If Called without arguments, return the type of problems this method can be used for: 
    # 1 - Classification problems (Discrete Class)
    # 2 - Regression problems (Continuous Class)
    # 3 - Both Regression and Classification problems (Continuous or Discrete Class)
    if res == None:
        return {"type":CLASSIFICATION}

    if res.numberOfIterations > 1:
        CAs = [[0.0] * res.numberOfIterations for i in range(res.numberOfLearners)]
        nIter = [0]*res.numberOfIterations
        for tex in res.results:
            ac = tex.actualClass
            nIter[tex.iterationNumber] += 1
            for i, cls in enumerate(tex.classes):
                if cls == ac:
                    CAs[i][tex.iterationNumber] += 1
        CAs = [[x/ni for x, ni in zip(y, nIter)] for y in CAs]

        CAfoldList = CAs
        CA = [statc.mean(x) for x in CAs]
        CAstd = stats.stdev(CAfoldList[0])

        if returnFoldStat:
            return [round(statc.mean(x),3) for x in CAs], CAfoldList
        else:
            return [round(statc.mean(x),3) for x in CAs]

    else:
        CAs = [0.0]*res.numberOfLearners
        for tex in res.results:
            CAs = map(lambda res, cls, ac = tex.actualClass:
                       res + types.IntType(cls == ac), CAs, tex.classes)
        return [round(x/(len(res.results)),3) for x in CAs]


##ecPA

def getRMSEstd(res, nFolds):
    """
    Method for calculating the std of RMSE of nFolds in a crossvalidation (returned).
    res is the object containing the results from orngTest methods such as crossValidation.
    """

    # Initialize a list to contain lists of errors for each fold.
    errorList = []
    for idx in range(nFolds):
        errorList.append([])

    # ex contains info on the fold number, prediction and actural responses for exah example used in the CV
    # Append ex error to correct fold list
    for ex in res.results:
         error = (ex.classes[0]- ex.actualClass)**2
         errorList[ex.iterationNumber].append(error)

    # RMSE of the different folds
    RMSElist = []
    for idx in range(nFolds):
        average =  sum(errorList[idx])/len(errorList[idx])
        RMSElist.append(math.sqrt(average))
    RMSEstd = stats.stdev(RMSElist)
    RMSEmean = statc.mean(RMSElist)
    if verbose > 0: print str(RMSEmean)+"\t"+str(RMSEstd)+"\t"+string.join( [str(x) for x in RMSElist], "\t")
    return RMSEstd, RMSElist


def WilcoxonRankTest(accLearner1, accLearner2):
    """
    The input is two list with the value pairs to be compared!
    Single sided Wilcoxon rank sum test. 
    See critical values: http://www.euronet.nl/users/warnar/demostatistiek/tables/WILCOXONTABEL.htm
    http://web.anglia.ac.uk/numbers/biostatistics/wilcoxon/local_folder/critical_values.html
    """
    # Learner 1 is the most accurate
    diffPlus = []
    # Learner 2 is the most accurate
    diffMinus = []
    for idx in range(len(accLearner2)):
        diff = accLearner1[idx]-accLearner2[idx]
        if diff > 0:
            diffPlus.append(abs(diff))
        elif diff < 0:
            diffMinus.append(abs(diff))
        else:
            diffPlus.append(abs(diff))
            diffMinus.append(abs(diff))
    diffPlus.sort()
    diffMinus.sort()

    # Rank the differences according to absolute values
    # R is a dictionary indexed by the rank number and with the values +, - or +/-
    # indicating which learner the rank number will be assigned to
    # The greater the diff the greater the rank idx
    R = {}
    for idx in range(len(accLearner1)):
        # Get the smallest value in each diff list (small diff -> small idx)
        try: diffPlusMin = diffPlus[0]
        except: diffPlusMin = 10000  # No more diffPlus elements, always take diffMinus
        try: diffMinusMin = diffMinus[0]
        except: diffMinusMin = 10000
        if diffPlusMin < diffMinusMin:
            if len(diffPlus) > 0: min = diffPlus.pop(0)
            R[str(idx)] = "+"
        elif diffPlusMin == diffMinusMin:
            if len(diffPlus) > 0: min = diffPlus.pop(0)
            if len(diffMinus) > 0: min = diffMinus.pop(0)
            R[str(idx)] = "+/-"
        else:
            if len(diffMinus) > 0: min = diffMinus.pop(0)
            R[str(idx)] = "-"

    # Get rank sums for the two learners - The greater the sum, the more accurate the learner
    Rplus = 0
    Rminus = 0
    for key, value in R.iteritems():
        if value == "+":
            Rplus = Rplus + int(key)
        elif value == "-":
            Rminus = Rminus + int(key)
        elif value == "+/-":
            Rplus = Rplus + (1.0/2)*int(key)
            Rminus = Rminus + (1.0/2)*int(key)

    Rlist = [Rplus, Rminus]
    # Does not work!!
    #print min(Rlist)
    Rlist.sort()
    # ***** Already in Orange - don't use the above *************
    #T = Rlist.pop(0)
    T = statc.wilcoxont(accLearner1, accLearner2)[0]
    N = len(R)

    print "Rank sum of learner 1"
    print Rplus
    print "Rank sum of learner 2"
    print Rminus
    if Rplus < Rminus:
        print "The hypothesis is that learner 2 is the most accurate"
    else:
        print "The hypothesis is that learner 1 is the most accurate"

    info = "If the number of data sets (N) is equal to 16 (our regression suite):\n"
    info += "N " + str(N) +"\n"
    info += "If T < 35 there is a 10% chance that the hypothesis is not true\n"
    info += "If T < 29 there is a 5% chance that the hypothesis is not true\n"
    info += "T " + str(T) + "\n"

    info += "If the number of data sets (N) is equal to 17 (our classification suite):\n"
    info += "N " + str(N) +"\n"
    info += "If T < 41 there is a 10% chance that the hypothesis is not true\n"
    info += "If T < 34 there is a 5% chance that the hypothesis is not true\n"
    info += "T " + str(T) + "\n"
    # If N > 20
    #z = (T - (1.0/4)*N*(N+1))/math.sqrt((1.0/24)*N*(N+1)*(2*N+1))
    #print z
    print info

    #Return the index of the best LEarner
    if Rplus < Rminus:
        return (1, info)
    else:
        return (0, info)


