from AZutilities import dataUtilities
from trainingMethods import AZorngRF
import Orange
import orange
import math
import copy

"""
Module for calculation of non conformity scores and the corresponding p-values and
conformal predictions for binary classifiers. 
getPvalue
	|
	|
	getScore
		|
		|
		{Methods to calculate the non-conf score}
"""


def meanStd(data):
    """ Calculate mean and standard deviation of data data[]: """
    length, mean, std = len(data), 0, 0
    for elem in data:
        mean = mean + elem
    mean = mean / float(length)
    for elem in data:
        std = std + (elem - mean) ** 2
    std = math.sqrt(std / float(length))
    mean = round(mean, 3)
    std = round(std, 3)
    return mean, std


def getScore(idx, extTrain, method = "minNN", maxDistRatio = None, measure = None):
    """
    Calculates non-conformity score for the example with index idx in the data set extTrain
    method:
    1) minNN - Get relative (all ex with diff labels) min distance in feature space from ex with idx in extTrain to the rest of extTrain with the same label as idx
    2) avgNN - average distance to 10 NN of the two diff classes
    """

    if method == "minNN":
        alpha = minNN(idx, extTrain, measure)
    elif method == "avgNN":
        alpha = avgNN(idx, extTrain, measure)
    elif method == "scaledMinNN":
        print "There is some problem with the scaling"
        alpha = minNN(idx, extTrain, maxDistRatio, measure) 
    elif method == "kNNratio":
        alpha = kNNratio(idx, extTrain, measure)
    elif method == "kNNratioStruct":
        alpha = kNNratioStruct(idx, extTrain, measure)
    elif method == "probPred":
        alpha = probPred(idx, extTrain)
    elif method == "LLOO":
        alpha = LLOO(idx, extTrain, measure)
    elif method == "LLOOprob":
        alpha = LLOOprob(idx, extTrain, measure)
    elif method == "LLOOprob_b":
        alpha = LLOOprob_b(idx, extTrain, measure)
    else:
        alpha = None
        print "Method not implemented"

    return alpha




def descRange(idx, extTrain):
    """
    Use the fraction of descriptors in the train set range.
    Not possible to use. Alpha must reflect the non-conformity with the rest of the train set with a given lable.
    Inside our outside the range is not predictive for which class the example belongs to. 
    """

    # Deselect example idx in extTrain
    idxList = range(0,idx)
    idxList.extend(range(idx+1,len(extTrain)))
    train = extTrain.get_items(idxList)

    # Get the idx example
    idxEx = extTrain.get_items([idx])

    # Loop over att attributes to see if the idxEx values are within the range of train
    outRangeCount = 0
    stat = Orange.statistics.basic.Domain(train)
    #print "%20s %5s %5s %5s" % ("feature", "min", "max", "avg")
    for a in stat:
        if a:
            #print "%20s %5.3f %5.3f %5.3f" % (a.variable.name, a.min, a.max, a.avg)
            #print idxEx[0][a.variable.name]
            idxValue = idxEx[0][a.variable.name]
            trainMin = a.min
            trainMax = a.max
            try:
                if idxValue < trainMin:
                    outRangeCount = outRangeCount + 1
                elif idxValue > trainMax:
                    outRangeCount = outRangeCount + 1
            except: pass

    alpha = float(outRangeCount)/len(extTrain.domain.attributes)

    return alpha


def probPred(idx, extTrain):
    """
    Use the RF prediction probability to set the non-conf score
    """
    attrList = ["SMILES_1"]
    extTrain = dataUtilities.attributeDeselectionData(extTrain, attrList)

    # Deselect example idx in extTrain
    idxList = range(0,idx)
    idxList.extend(range(idx+1,len(extTrain)))
    train = extTrain.get_items(idxList)

    # Train a model
    model = AZorngRF.RFLearner(train)

    # Predict example idx
    predList = model(extTrain[idx], returnDFV = True)
    pred = predList[0].value
    prob = predList[1]
    actual = extTrain[idx].get_class().value
    #print pred, actual, prob

    # More non conforming if prediction is different from actual label
    if pred != actual:
        alpha = 1.0 + abs(prob)
    else:
        alpha = 1.0 - abs(prob)

    #print alpha
    return alpha


def minNN(idx, extTrain, maxDistRatio = None, measure = None):
    """
    Use the ratio between the distance to the nearest neighbor of the same and of the other class
    Two versions exist, with and without scaling with the max distance ratio within the train set. 
    """

    attrList = ["SMILES_1"]
    extTrain = dataUtilities.attributeDeselectionData(extTrain, attrList)

    distListSame = [] 
    distListDiff = [] 
    #measure = Orange.distance.Euclidean(extTrain)
    if not measure:
        measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            dist = measure(extTrain[idx], extTrain[runIdx]) 
            if extTrain[idx].get_class().value == extTrain[runIdx].get_class().value:
                distListSame.append(dist)
            else:
                distListDiff.append(dist)
    minDistSame = min(distListSame)
    minDistDiff = min(distListDiff)
    if minDistDiff == 0:
        if maxDistRatio:
            alpha = 1.0
        else:
            alpha = max(distListDiff)
    else:
        if maxDistRatio:
            alpha = minDistSame/(float(minDistDiff)*maxDistRatio)
        else:
            alpha = minDistSame/float(minDistDiff)

    #fid = open("tempFile.txt", "a") 
    #fid.write(str(minDistSame)+"\t"+str(minDistDiff)+"\t"+str(maxDistRatio)+"\t"+str(alpha)+"\n")
    #fid.close()

    return alpha
        

def avgNN(idx, extTrain, measure = None):
    """
    Use the ratio between the distance to the kNN of the same and of the other class
    """
    attrList = ["SMILES_1"]
    extTrain = dataUtilities.attributeDeselectionData(extTrain, attrList)

    distListSame = [] 
    distListDiff = [] 
    #measure = Orange.distance.Euclidean(extTrain)
    if not measure:
        measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            dist = measure(extTrain[idx], extTrain[runIdx]) 
            if extTrain[idx].get_class().value == extTrain[runIdx].get_class().value:
                distListSame.append(dist)
            else:
                distListDiff.append(dist)
    distListSame.sort()
    avgSame = sum(distListSame[0:10])/10.0
    distListDiff.sort()
    avgDiff = sum(distListDiff[0:10])/10.0
    if avgDiff == 0:
        alpha = max(distListDiff)
    else:
        alpha = avgSame/float(avgDiff)

    return alpha


def kNNratio(idx, extTrain, measure = None):
    """
    Use the fraction of kNN with the same response.
    """
    attrList = ["SMILES_1"]
    extTrain = dataUtilities.attributeDeselectionData(extTrain, attrList)

    distList = [] 
    if not measure:
        #measure = instances.MahalanobisConstructor(extTrain)
        measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            dist = measure(extTrain[idx], extTrain[runIdx]) 
            distList.append(dist)

    # Get the distance of the 10th NN
    distList.sort()
    thresDist = distList[9]

    # Find the labels of the 10 NN
    sameCount = 0
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            dist = measure(extTrain[idx], extTrain[runIdx]) 
            if dist <= thresDist:
                if extTrain[idx].get_class().value == extTrain[runIdx].get_class().value:
                    sameCount = sameCount + 1
    alpha = 1.00 - float(sameCount)/10.0

    return alpha


def kNNratioInd(train, calSet, measure = None):
    """
    Use the fraction of kNN with the same response.
    """
    if not measure:
        #measure = instances.MahalanobisConstructor(extTrain)
        measure = orange.ExamplesDistanceConstructor_Euclidean(train)

    alphaList = []
    for predEx in calSet:
        distList = [] 
        for runIdx in range(len(train)):
            dist = measure(predEx, train[runIdx]) 
            distList.append(dist)

        # Get the distance of the 10th NN
        distList.sort()
        thresDist = distList[9]

        # Find the labels of the 10 NN
        sameCount = 0
        for runIdx in range(len(train)):
            dist = measure(predEx, train[runIdx]) 
            if dist <= thresDist:
                if predEx.get_class().value == train[runIdx].get_class().value:
                    sameCount = sameCount + 1
        alpha = 1.00 - float(sameCount)/10.0
        alphaList.append(alpha)

    return alphaList, train


def kNNratioStruct(idx, extTrain, measure = None):
    """
    Use the fraction of kNN with the same response.
    """

    from rdkit import Chem
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit import DataStructs

    # Daylight like fp
    smiles = extTrain[idx]["SMILES_1"].value
    mol = Chem.MolFromSmiles(smiles)
    fp = FingerprintMols.FingerprintMol(mol)

    simList = [] 
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            smiles0 = extTrain[runIdx]["SMILES_1"].value
            mol0 = Chem.MolFromSmiles(smiles0)
            fp0 = FingerprintMols.FingerprintMol(mol0)
            tanSim = DataStructs.FingerprintSimilarity(fp,fp0)
            simList.append(tanSim)

    # Get the distance of the 10th NN
    simList.sort(reverse = True)
    thresDist = simList[9]

    # Find the labels of the 10 NN
    sameCount = 0
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            smiles0 = extTrain[runIdx]["SMILES_1"].value
            mol0 = Chem.MolFromSmiles(smiles0)
            fp0 = FingerprintMols.FingerprintMol(mol0)
            tanSim = DataStructs.FingerprintSimilarity(fp,fp0)
            if tanSim >= thresDist:
                if extTrain[idx].get_class().value == extTrain[runIdx].get_class().value:
                    sameCount = sameCount + 1
    alpha = 1.00 - float(sameCount)/10.0

    return alpha


def LLOO(idx, extTrain, measure = None):
    """
    Use the fraction of kNN correctly predicted by a local model
    Hard coded to 20 NN.
    Modeling method. RF of Tree?
    """
    attrList = ["SMILES_1"]
    extTrain = dataUtilities.attributeDeselectionData(extTrain, attrList)

    distList = [] 
    if not measure:
        measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            dist = measure(extTrain[idx], extTrain[runIdx]) 
            distList.append(dist)

    # Get the distance of the 20th NN
    distList.sort()
    thresDist = distList[19]

    # Find the labels of the 20 NN
    kNN = []
    for runIdx in range(len(extTrain)):
        dist = measure(extTrain[idx], extTrain[runIdx]) 
        if dist <= thresDist:
            kNN.append(extTrain[runIdx])
    kNNtrain = dataUtilities.DataTable(kNN)
    
    # Find the fraction of correctly predicted ex in a LOO over kNN
    corrPred = 0
    for idx in range(len(kNNtrain)):

        # Deselect example idx in extTrain
        idxList = range(0,idx)
        idxList.extend(range(idx+1,len(kNNtrain)))
        train = kNNtrain.get_items(idxList)

        # Train a model
        model = AZorngRF.RFLearner(train)
        #model = Orange.classification.tree.TreeLearner(train)

        pred = model(kNNtrain[idx]).value
        actual = kNNtrain[idx].get_class().value
        if pred == actual:
            corrPred = corrPred + 1
    alpha = 1.0 - float(corrPred)/len(kNNtrain)

    return alpha


def LLOOprob(idx, extTrain, measure = None):
    """
    Use the fraction of kNN correctly predicted by a local model
    Hard coded to 20 NN.
    Modeling method. RF of Tree?
    """

    distList = [] 
    if not measure:
        measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            dist = measure(extTrain[idx], extTrain[runIdx]) 
            distList.append(dist)

    # Get the distance of the 20th NN
    distList.sort()
    thresDist = distList[50]    # Smaller number of NN does not work with returnDFV

    # Find the predEx and the 20 NN
    kNN = []
    for runIdx in range(len(extTrain)):
        dist = measure(extTrain[idx], extTrain[runIdx]) 
        if dist <= thresDist:
            kNN.append(extTrain[runIdx])
    kNNtrain = dataUtilities.DataTable(kNN)
    
    # Find the fraction of correctly predicted ex in a LOO over kNN
    alphaList = []
    for iidx in range(len(kNNtrain)):

        # Deselect example idx in extTrain
        idxList = range(0,iidx)
        idxList.extend(range(iidx+1,len(kNNtrain)))
        train = kNNtrain.get_items(idxList)

        # Get prediction and pred probability 
        model = AZorngRF.RFLearner(train)
        predList = model(kNNtrain[iidx], returnDFV = True)
        pred = predList[0].value
        prob = predList[1]
        actual = kNNtrain[iidx].get_class().value
        # alpha should be greater the less certain the model
        try:
            if pred != actual:
                alpha = 1.0 + abs(prob)
            else:
                alpha = 1.0 - abs(prob)
            alphaList.append(alpha)
        except: pass

    alpha = sum(alphaList)/float(len(alphaList))

    return alpha


def LLOOprob_b(idx, extTrain, measure = None):
    """
    Use the fraction of kNN correctly predicted by a local model
    Hard coded to 50 NN.
    Modeling method. RF of Tree?
    """

    distList = [] 
    if not measure:
        measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for runIdx in range(len(extTrain)):
        if runIdx != idx:
            dist = measure(extTrain[idx], extTrain[runIdx]) 
            distList.append(dist)

    # Get the distance of the 50th NN
    distList.sort()
    thresDist = distList[50]    # Smaller number of NN does not work with returnDFV

    # Find the predEx and the 20 NN
    kNN = []
    for runIdx in range(len(extTrain)):
        dist = measure(extTrain[idx], extTrain[runIdx]) 
        if dist <= thresDist:
            kNN.append(extTrain[runIdx])
    kNNtrain = dataUtilities.DataTable(kNN)
    
    # Find the fraction of correctly predicted ex in a LOO over kNN
    alphaList = []
    alphaEx = 0
    for iidx in range(len(kNNtrain)):

        # Deselect example idx in extTrain
        idxList = range(0,iidx)
        idxList.extend(range(iidx+1,len(kNNtrain)))
        train = kNNtrain.get_items(idxList)

        # Get prediction and pred probability 
        model = AZorngRF.RFLearner(train)
        predList = model(kNNtrain[iidx], returnDFV = True)
        pred = predList[0].value
        prob = predList[1]
        actual = kNNtrain[iidx].get_class().value

        # The prob of the predEx is more important
        dist = measure(extTrain[idx], kNNtrain[iidx]) 

        # alpha should be greater the less certain the model
        try:
            if pred != actual:
                alpha = 1.0 + abs(prob)
                if dist < 0.001:
                    alphaEx = alpha
            else:
                alpha = 1.0 - abs(prob)
                if dist < 0.001:
                    alphaEx = alpha
            alphaList.append(alpha)
        except: pass

    alpha = alphaEx + sum(alphaList)/float(len(alphaList))

    return alpha


def getMeanStd(extTrain):

    # Get the min dist for all ex in the data set
    minSame = []
    minDiff = []
    measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for idx in range(len(extTrain)):
        distListSame = []
        distListDiff = []
        for iidx in range(len(extTrain)):
            if idx != iidx:
                dist = measure(extTrain[idx], extTrain[iidx]) 
                if extTrain[idx].get_class().value == extTrain[iidx].get_class().value:
                    distListSame.append(dist)
                else:
                    distListDiff.append(dist)
        minSame.append(min(distListSame))
        minDiff.append(min(distListDiff))

    # Calculate mean and std of all the min distances
    meanSame, stdSame = meanStd(minSame)
    meanDiff, stdDiff = meanStd(minDiff)

    return meanSame, stdSame, meanDiff, stdDiff


def getMinDistRatio(train):
    """
    Calculate the minDistSame and minDistDiff ratio for all ex in the data set and select the greatest quotient.
    Used to scale the minDist ratios in the non-conf score.
    """

    # Get the min dist for all ex in the data set
    minSame = []
    minDiff = []
    minRatio = []
    measure = orange.ExamplesDistanceConstructor_Euclidean(extTrain)
    for idx in range(len(train)):
        distListSame = []
        distListDiff = []
        for iidx in range(len(train)):
            if idx != iidx:
                dist = measure(train[idx], train[iidx]) 
                if train[idx].get_class().value == train[iidx].get_class().value:
                    distListSame.append(dist)
                else:
                    distListDiff.append(dist)
        minSame.append(min(distListSame))
        minDiff.append(min(distListDiff))
        if min(distListDiff) == 0:
            alpha = max(distListDiff)
        else:
            minRatio.append(min(distListSame)/float(min(distListDiff)))

    # Calculate min, mean and std of all the min distances
    meanSame, stdSame = meanStd(minSame)
    meanDiff, stdDiff = meanStd(minDiff)
    maxDistRatio = max(minRatio)

    return maxDistRatio 
        

def getPvalue(train, predEx, label, method = "avgNN", measure = None):
    """
    method; avgNN, scaledMinNN, minNN, kNNratio
    """

    # Set label to class of predEx
    newPredEx = Orange.data.Table(predEx.domain, [predEx])
    newPredEx[0][newPredEx.domain.classVar] = label

    # Add predEx to train
    extTrain = dataUtilities.concatenate([train, newPredEx])
    extTrain = extTrain[0]

    # Calculate a non-conf score for each ex in train + predEx with given label
    if method == "scaledMinNN":
        # Calculate average and std of min distanses in train set
        maxDistRatio = getMinDistRatio(train)
    nonConfList = []
    for idx in range(len(extTrain)):
        if method == "scaledMinNN":
            alpha = getScore(idx, extTrain, method, maxDistRatio)
        else:
            alpha = getScore(idx, extTrain, method, None, measure)
        nonConfList.append(alpha)
        #if idx == 1: 
        #    print "Breaking after one ex!!"
        #    break

    nonConfListSorted = copy.deepcopy(nonConfList)
    nonConfListSorted.sort()
    fid = open("NonConf.txt", "w")
    for ex in nonConfListSorted:
        fid.write(str(ex)+"\n")
    fid.close()

    # The last non-conf score is that of predEx
    # The p-value is the fraction of ex with alpha gt that of predEx
    trainList = nonConfList[0:len(nonConfList)-1]
    alphaPredEx = nonConfList[len(nonConfList)-1] 
    moreNonConfList = []
    for score in trainList:
        if score > alphaPredEx:
            moreNonConfList.append(score)
    pvalue = len(moreNonConfList)/float(len(trainList))

    return pvalue


def printResults(pvalues, labels, actualLabel, method, resultsFile):

    confLevel = 0.95
    
    #print "OBS! Assuming two labels!!!"
    #print "Confidence level in predicting label ", labels[0]
    conf1 = round(1-pvalues[1], 3)

    #print "Confidence level in predicting label ", labels[1]
    conf2 = round(1-pvalues[0], 3)

    #print "Requiering 95% confidence gives "
    if conf1 > confLevel and conf2 < confLevel:
    #    print "Label ", labels[0], "is predicted with at least 95% conf"
        prediction = labels[0]
    elif conf1 < confLevel and conf2 > confLevel:
    #    print "Label ", labels[1], "is predicted with at least 95% conf"
        prediction = labels[1]
    elif conf1 < confLevel and conf2 < confLevel:
    #    print "No prediction can be given at 95% confidence. Both?"
        prediction = "Both"
    else:           # if conf1 > confLevel and conf2 > confLevel:
    #    print "Predicting both labels. Empty?"
        prediction = "Empty"

    fid = open(resultsFile, "a")
    fid.write(actualLabel+"\t"+labels[0]+"\t"+labels[1]+"\t"+str(pvalues[0])+"\t"+str(pvalues[1])+"\t"+str(conf1)+"\t"+str(conf2)+"\t"+prediction+"\n")
    fid.close()

    return prediction


def printStat(resDict, labels):

    # Print statistics
    T0 = 0
    T1 = 0
    F0 = 0
    F1 = 0
    Both = 0
    Empty = 0
    for key, values in resDict.iteritems():
        if values["actualLabel"] == labels[0]:
            if values["actualLabel"] == values["prediction"]:
                T0 = T0 + 1
            elif  values["prediction"] == "Both":
                Both = Both + 1
            elif  values["prediction"] == "Empty":
                Empty = Empty + 1
            elif values["prediction"] == labels[1]:
                F1 = F1 + 1
        if values["actualLabel"] == labels[1]:
            if values["actualLabel"] == values["prediction"]:
                T1 = T1 + 1
            elif  values["prediction"] == "Both":
                Both = Both + 1
            elif  values["prediction"] == "Empty":
                Empty = Empty + 1
            elif values["prediction"] == labels[0]:
                F0 = F0 + 1
    print "True ", labels[0], ": ", T0
    print "True ", labels[1], ": ", T1
    print "False ", labels[0], ": ", F0
    print "False ", labels[1], ": ", F1
    print "Both: ", Both
    print "Empty: ", Empty


def getRFAcc(train, work):

    model = AZorngRF.RFLearner(train)
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    for ex in work:
        pred = model(ex).value
        actual = ex.get_class().value
        if actual == "A":
            if pred == "A":
                TP = TP + 1
            else:
                FN = FN + 1
        elif actual == "N":
            if pred == "N":
                TN = TN + 1
            else:
                FP = FP + 1
    print "TP\tTN\tFP\tFN\n"    
    print str(TP)+"\t"+str(TN)+"\t"+str(FP)+"\t"+str(FN)+"\n"    

    fid = open("RFresults.txt", "a")
    fid.write(str(TP)+"\t"+str(TN)+"\t"+str(FP)+"\t"+str(FN)+"\n")
    fid.close()


def getRFprobAcc(train, work, probThres):

    model = AZorngRF.RFLearner(train)
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    noPred = 0
    for ex in work:
        actual = ex.get_class().value
        predList = model(ex, returnDFV = True)
        pred = predList[0].value
        prob = predList[1]

        if abs(prob) > probThres:
            if actual == "A":
                if pred == "A":
                    TP = TP + 1
                else:
                    FN = FN + 1
            elif actual == "N":
                if pred == "N":
                    TN = TN + 1
                else:
                    FP = FP + 1
        else: 
            noPred = noPred + 1
    print "TP\tTN\tFP\tFN\tnoPred\n"    
    print str(TP)+"\t"+str(TN)+"\t"+str(FP)+"\t"+str(FN)+"\t"+str(noPred)+"\n"    

    fid = open("RFprob"+str(probThres)+"Results.txt", "a")
    fid.write(str(TP)+"\t"+str(TN)+"\t"+str(FP)+"\t"+str(FN)+"\t"+str(noPred)+"\n")
    fid.close()


def getProbPredAlpha(model, ex):
        predList = model(ex, returnDFV = True)
        pred = predList[0].value
        prob = predList[1]
        actual = ex.get_class().value

        # More non conforming if prediction is different from actual label
        if pred != actual:
            alpha = 1.0 + abs(prob)
        else:
            alpha = 1.0 - abs(prob)

        return alpha


def probPredInd(trainSet, calSet):
    """
    Use the RF prediction probability to set the non-conf score
    """
    attrList = ["SMILES_1"]
    trainSet = dataUtilities.attributeDeselectionData(trainSet, attrList)

    # Train a model
    model = AZorngRF.RFLearner(trainSet)

    # Get the list of NC for all ex in calSet
    alphaList = []
    for ex in calSet:
        alpha = getProbPredAlpha(model, ex)
        alphaList.append(alpha)

    return alphaList, model


def getScores(trainSet, calSet, method):

    #if method == "minNN":
    #    alpha = minNN(idx, extTrain, measure)
    #elif method == "avgNN":
    #    alpha = avgNN(idx, extTrain, measure)
    #elif method == "scaledMinNN":
    #    print "There is some problem with the scaling"
    #    alpha = minNN(idx, extTrain, maxDistRatio, measure) 
    if method == "kNNratio":
        alphaList, model = kNNratioInd(trainSet, calSet)
    #elif method == "kNNratioStruct":
    #    alpha = kNNratioStruct(idx, extTrain, measure)
    elif method == "probPred":
        alphaList, model = probPredInd(trainSet, calSet)
    #elif method == "LLOO":
    #    alpha = LLOO(idx, extTrain, measure)
    #elif method == "LLOOprob":
    #    alpha = LLOOprob(idx, extTrain, measure)
    #elif method == "LLOOprob_b":
    #    alpha = LLOOprob_b(idx, extTrain, measure)
    #else:
    #    alpha = None
    #    print "Method not implemented"

    return alphaList, model


def getIndPvalue(model, NClist, predEx, label, method = "avgNN", measure = None):
    """
    method; avgNN, scaledMinNN, minNN, kNNratio
    """

    # Set label to class of predEx
    newPredEx = Orange.data.Table(predEx.domain, [predEx])
    newPredEx[0][newPredEx.domain.classVar] = label
 
    # Calculate the NC score of predEx
    if method == "probPred":
        alpha = getProbPredAlpha(model, newPredEx[0])
    elif method == "kNNratio":
        alphaList, model = kNNratioInd(model, newPredEx)  # model is the train set for NN methods
        alpha = alphaList[0]

    # The p-value is the fraction of ex with alpha gt that of predEx
    moreNonConfList = []
    for score in NClist:
        if score > alpha:
            moreNonConfList.append(score)
    pvalue = len(moreNonConfList)/float(len(NClist))

    return pvalue


def getConfPred(train, work, method, measure = None, resultsFile = "CPresults.txt", verbose = False):
    """
    method - non-conformity score method
    """

    # Get conformal predictions
    resDict = {}
    idx = 0
    for predEx in work:
        labels = predEx.domain.classVar.values
        pvalues = []
        for label in labels:
            if method == "combo":
                pvalue1 = getPvalue(train, predEx, label, "kNNratio", measure)
                pvalue2 = getPvalue(train, predEx, label, "probPred")
                pvalue = (pvalue1 + pvalue2)/2.0
            else:
                pvalue = getPvalue(train, predEx, label, method, measure)
            pvalues.append(pvalue)
        actualLabel = predEx.get_class().value
        prediction = printResults(pvalues, labels, actualLabel, method, resultsFile)
        idx = idx + 1
        resDict[idx] = {"actualLabel": actualLabel, "prediction": prediction}

        #print "Break after the first example"
        #if idx == 1: break

    if verbose:
        printStat(resDict, labels)


def getIndConfPred(train, work, method, measure = None, resultsFile = "CPresults.txt", verbose = False):
    """
    Partition train into a training and a calibration set (10%). Use the non-conf scores of the cal set to predict 
    all examples in work.
    method - non-conformity score method
    """

    # Randomily select 10% of train as a cal set
    indices2 = Orange.data.sample.SubsetIndices2(p0=0.10)
    ind = indices2(train)
    calSet = train.select(ind, 0)
    trainSet = train.select(ind, 1)

    # Calculate NC for the calibration set
    NClist, model = getScores(trainSet, calSet, method)

    # Calculate p-values for all ex in work
    resDict = {}
    idx = 0
    for predEx in work:
        labels = predEx.domain.classVar.values
        pvalues = []
        for label in labels:
            if method == "combo":
                pvalue1 = getIndPvalue(model, NClist, predEx, label, "kNNratio", measure)
                pvalue2 = getIndPvalue(model, NClist, predEx, label, "probPred")
                pvalue = (pvalue1 + pvalue2)/2.0
            else:
                pvalue = getIndPvalue(model, NClist, predEx, label, method, measure)
            pvalues.append(pvalue)
        actualLabel = predEx.get_class().value
        prediction = printResults(pvalues, labels, actualLabel, method, resultsFile)
        idx = idx + 1
        resDict[idx] = {"actualLabel": actualLabel, "prediction": prediction}

        #print "Break after the first example"
        #if idx == 1: break


if __name__ == "__main__":
    """
    Assumptions;
    Binary classification 
    Class labels not generalized, assumed to be 'A' and 'N'

    This main will test the implemented CP methods in a 10 fold CV
    """

    data = dataUtilities.DataTable("trainData.tab")
    descList = ["SMILES", "SMILES_1"]
    data = dataUtilities.attributeDeselectionData(data, descList)

    print "Please note that the class labels are not generalized and need to be checked for a new data set"
    print "Assumed to be A and N"
    methods = ["kNNratio", "minNN", "avgNN", "probPred", "combo", "LLOO", "LLOOprob"]   # Non-conformity score method
    #methods = ["kNNratio"]
    cpMethod = "transductive"   # inductive or transductive

    #print "Temp position to save comp time!!"
    # Append to python path /home/kgvf414/dev/AZOrange0.5.5/orangeDependencies/src/orange/orange/Orange/distance/
    #import instances
    #measure = instances.MahalanobisConstructor(data)
    measure = None
    methodIdx = 1
    for method in methods:

        resultsFile = "CPresults_"+method+".txt"
        fid = open(resultsFile, "w")
        fid.write("ActualLabel\tLabel1\tLabel2\tPvalue1\tPvalue2\tConf1\tConf2\tPrediction\n")
        fid.close()

        # Run a 10 fold CV
        nFolds = 10
        ind = Orange.data.sample.SubsetIndicesCV(data, nFolds)
        for idx in range(nFolds):
            work = data.select(ind, idx)
            train = None
            for iidx in range(nFolds):
                if iidx != idx:
                    if not train:
                        train = data.select(ind, iidx)
                    else:
                        train.extend(data.select(ind, iidx))
            print "Length of train ", len(train)
            print "Length of work ", len(work)

            # Create results file and get the conformal predictions
            if cpMethod == "transductive":
                getConfPred(train, work, method, measure, resultsFile, verbose = True)
            elif cpMethod == "inductive":
                print "Please note, only kNNratio and probPred implemented for ICP!"
                getIndConfPred(train, work, method, measure, resultsFile, verbose = True)
            else:
                print "Valid cpMethod values are 'transductive' and 'inductive'"

            # Compare with RF predictions
            if methodIdx == 1:
                getRFAcc(train, work)
                prob = 0.1
                getRFprobAcc(train, work, prob)
                prob = 0.2
                getRFprobAcc(train, work, prob)
                prob = 0.3
                getRFprobAcc(train, work, prob)

            #print "Breaking after first fold"
            #if idx == 0: break
        methodIdx = methodIdx + 1      






