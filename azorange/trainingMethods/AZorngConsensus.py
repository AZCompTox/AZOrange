"""
AZorngonsensus
Module to build a Consensus model.
"""
import glob
import string
import os,sys
import AZBaseClasses
import orange
import sys
import re

from cStringIO import StringIO
from tokenize import generate_tokens

from AZutilities import dataUtilities
#from AZutilities import miscUtilities

import AZOrangeConfig as AZOC

class ConsensusLearner(AZBaseClasses.AZLearner):
    """
    Creates a Consensus as an Orange type of learner instance. 
    """

    def __new__(cls, trainingData = None, name = "Consensus learner", **kwds):
        self = AZBaseClasses.AZLearner.__new__(cls, **kwds)
        if trainingData:
            self.__init__(name, **kwds)
            return self.__call__(trainingData)
        else:
            self.__dict__.update(kwds)
            self.name = name
            return self

    def __init__(self, name = "Consensus learner", **kwds):
        """
        Set default values of the model parameters if they are not given as inputs
        """
        #NOTE: Only use one of the learnersNames  or learnersObj. learnersNames will have priority on learnersObj!
        self.learnersNames = None       # The Learners names that we want to include in out consensus model. Only <name> in AZorng<name> will be allowed
        self.learnersObj = None         # The objects handling the desired learners. Can be direcly specifyed or will be constructed from learnersNames
        self.rejectedLearners = []      # Learners that, for some reason, could not be used/trained
        self.regressionExpression = None# Expression used with the regression model
        self.discreteExpression = None  # Expression used with the discrete model
        self.learnerNameMap = None      # Connects a variable, used in an expression, with a learner name 
        self.learnerObjMap = None       # Connects a variable, used in an expression, with a learner object
        self.name = name
        self.verbose = 0 
        self.imputeData = None
        self.NTrainEx = 0
        self.basicStat  = None
        # Append arguments to the __dict__ member variable
        self.__dict__.update(kwds)


    def __call__(self, trainingData, weight = None):
        """Creates a Consensus model from the data in trainingData. """
        if not AZBaseClasses.AZLearner.__call__(self,trainingData, weight) or not trainingData:
            return None

        # Make sure the correct number of arguments are supplied
        if not self.learnersNames and not self.learnersObj:
            if not self.regressionExpression and not self.discreteExpression:
                return None
            
            if not self.learnerNameMap and not self.learnerObjMap:
                return None

        #Create the Learners
        if self.learnersNames:
            self.rejectedLearners = []
            self.learnersObj = [] 
            for learner in self.learnersNames:
                try:
                    exec("from trainingMethods import AZorng"+learner)
                    self.learnersObj.append(eval("AZorng"+learner+"."+learner+"Learner()"))
                except:
                    self.rejectedLearners.append(learner)

            for learner in self.rejectedLearners:
                self.learnersNames.remove(learner)
                
            if len(self.rejectedLearners) > 0:
                print "WARNING: There werte some learners that were rejected:\n"
                for l in self.rejectedLearners:
                    print "   ",l,"\n"

        if self.learnersObj:
            if len(self.learnersObj) <= 1:
                raise Exception("ERROR: The Consensus model needs at least 2 valid learners.\n"+\
                                "Learners: "+str(self.learnersObj))

        if self.learnerNameMap:
            self.rejectedLearners = []
            self.learnerObjMap = {}
            for learner in self.learnerNameMap:
                try:
                    exec("from trainingMethods import AZorng" + self.learnerNameMap[learner])
                    self.learnerObjMap[learner] =  eval("AZorng" + self.learnerNameMap[learner] + "." + self.learnerNameMap[learner] + "Learner()")
                except:
                    self.rejectedLearners.append(learner)

            for learner in self.rejectedLearners:
                del self.learnerNameMap[learner]

            if len(self.rejectedLearners) > 0:
                print "WARNING: There were some learners that were rejected:\n"
                for l in self.rejectedLearners:
                    print "   ",l,"\n"

        if self.learnerObjMap:
            if len(self.learnerObjMap) <= 1:
                raise Exception("ERROR: The Consensus model needs at least 2 valid learners.\n"+\
                                "Learners: "+str(self.learnerObjMap))

        
        if trainingData.domain.classVar.varType == orange.VarTypes.Discrete and len(trainingData.domain.classVar.values) != 2:
            raise Exception("ERROR: The Consensus model only supports binary classification or regression problems.")
       
        # Call the train method
        if not self.regressionExpression and not self.discreteExpression:
            # Default behaviour, no expression defined.
            classifiers = []
            for learner in self.learnersObj:
                classifiers.append(learner(trainingData))
                if not classifiers[-1]:
                    if self.verbose > 0:
                        print "ERROR: Could not create the model ",str(learner)

                    return None
                else:
                    
                    #Try to get the imputeData, basicStat from a model that have it!
                    if hasattr(classifiers[-1], "basicStat") and classifiers[-1].basicStat and not self.basicStat:
                        self.basicStat = classifiers[-1].basicStat
                    
                    if hasattr(classifiers[-1], "NTrainEx") and classifiers[-1].basicStat and not self.NTrainEx:
                        self.NTrainEx = len(trainingData)
                        
                    if hasattr(classifiers[-1], "imputeData") and classifiers[-1].imputeData and not self.imputeData:
                        self.imputeData = classifiers[-1].imputeData
                            
            return ConsensusClassifier(classifiers = classifiers,
                                       classVar = trainingData.domain.classVar,
                                       verbose = self.verbose,
                                       domain = trainingData.domain,
                                       varNames = [attr.name for attr in trainingData.domain.attributes],
                                       NTrainEx = self.NTrainEx,
                                       basicStat = self.basicStat,
                                       imputeData = self.imputeData)
        else:
            classifiers = {}
            for learner in self.learnerObjMap:
                newClassifier = self.learnerObjMap[learner](trainingData)

                if not newClassifier:
                    if self.verbose > 0:
                        print "ERROR: Could not create the model ",str(learner)

                    return None
                else:
                    classifiers[learner] = newClassifier
                    #Try to get the imputeData, basicStat from a model that have it!
                    if hasattr(newClassifier, "basicStat") and newClassifier.basicStat and not self.basicStat:
                        self.basicStat = newClassifier.basicStat
                    
                    if hasattr(newClassifier, "NTrainEx") and newClassifier.basicStat and not self.NTrainEx:
                        self.NTrainEx = len(trainingData)
                        
                    if hasattr(newClassifier, "imputeData") and newClassifier.imputeData and not self.imputeData:
                        self.imputeData = newClassifier.imputeData
                            
            return ConsensusClassifier(classifiers = classifiers,
                                       discreteExpression = self.discreteExpression,
                                       regressionExpression = self.regressionExpression,
                                       classVar = trainingData.domain.classVar,
                                       verbose = self.verbose,
                                       domain = trainingData.domain,
                                       varNames = [attr.name for attr in trainingData.domain.attributes],
                                       NTrainEx = self.NTrainEx,
                                       basicStat = self.basicStat,
                                       imputeData = self.imputeData)


class ConsensusClassifier(AZBaseClasses.AZClassifier):
    def __new__(cls, name = "Consensus classifier", **kwds):
        self = AZBaseClasses.AZClassifier.__new__(cls, name = name,  **kwds)        
        return self
    def __init__(self, name = "Consensus classifier", **kwds):
        self.verbose = 0
        varNames = None
        self.domain = None
        self.classVar = None
        self.__dict__.update(kwds)
        self._isRealProb = False
        self.name = name
        self.status = ""
        if not self.classVar or not self.domain: self.setDomainAndClass()

    def __call__(self, origExample = None, resultType = orange.GetValue, returnDFV = False):
        """
        orange.GetBoth -          <type 'tuple'>                     ->    (<orange.Value 'Act'='3.44158792'>, <3.442: 1.000>)
        orange.GetValue -         <type 'orange.Value'>              ->    <orange.Value 'Act'='3.44158792'>
        orange.GetProbabilities - <type 'orange.DiscDistribution'>   ->    <0.000, 0.000>
        returnDFV - Flag indicating to return the Decision Function Value. If set to True, it will encapsulate the original result asked by the keyword resultType and the DFV into a tuple:
                ((<orange.Value 'Act'='3.44158792'>, <3.442: 1.000>), 0.34443)
                (<orange.Value 'Act'='3.44158792'>, 0.34443) 
                (<0.000, 0.000>, 0.34443)
                If it is not a binary classifier, DFV will be equal to None
                DFV will be a value from -0.5 to 0.5
        """
        self.status = ""
        if not self.classVar or not self.domain:
            self.setDomainAndClass()

        if type(self.classifiers).__name__ == 'list':
            if origExample == None:
                return self.classifiers[0](None, resultType)
            else:
                # Predict using the models  
                predictions = []   #The individual predictions for each repective classifier in classifiers
                predicted = None
                probabilities = None
                DFV = None
                if self.classVar.varType == orange.VarTypes.Discrete:
                    if (len(self.classifiers) % 2) != 0:     #Odd number of Classifiers, Using the Majority
                        self.status = "Using Majority (Odd number of classifiers)"
                        if self.verbose:
                            print self.status
                        votes = {self.classVar.values[0]:0, self.classVar.values[1]:0} # We already assured that if var is discrete, it is a binary classifyer
                        for c in self.classifiers:
                            predictions.append(c(origExample))
                            votes[predictions[-1].value] += 1
                        if self.verbose: 
                            print str([cf.name for cf in self.classifiers])[1:-1].replace(",","\t")
                            print str(predictions)[1:-1].replace(",","\t")
                        invVotes = dict(map(lambda item: (item[1],item[0]),votes.items()))
                        maxVoted = invVotes[max(invVotes.keys())]
                        predicted = self.classVar[self.classVar.values.index(maxVoted)]
                        # generate probabilities based on voting. We assured already this is a binary classifier
                        self._isRealProb = True
                        probOf1 = votes[self.classVar.values[1]]/len(self.classifiers)
                        DFV = self.convert2DFV(probOf1)
                        probabilities = self.__getProbabilities(probOf1)
                    else:     #Even number of Classifiers: Use the average of the N probabilities.
                              #if at least one of the Classifiers do not support true probabilities, STOP!       
                        self.status = "Using probabilities average (Even number of classifiers)"
                        if self.verbose: print self.status                     
                        for c in self.classifiers:
                            FailedPredicting = False
                            try:
                                predictions.append(c(origExample, orange.GetBoth))
                            except:
                                print "Learner ",c,"couldn't predict: ",sys.exc_info()
                                FailedPredicting = True
                            if FailedPredicting or (hasattr(c,"isRealProb")  and not c.isRealProb()):
                                msg = "ERROR: Cannot predict because the learner "+(hasattr(c,"name") and str(c.name) or str(c))+" does not return real probabilities.\n"+\
                                "       In order to use the Consensus model with this dataset, do one of the following:\n"+\
                                "         a) Train the Consensus model with an odd number of learners so that voting can be done instead.\n"+\
                                "         b) Train the consensus model with learners that support real probabilities."
                                raise Exception(msg)
                        self._isRealProb = True
                        overallProb = {self.classVar.values[0]:0, self.classVar.values[1]:0} 
                        for p in predictions:
                            overallProb[p[0].value] += max(p[1]) #For each predicted value of a learner, add the repective probability
                        for prob in overallProb:
                            overallProb[prob] = overallProb[prob]/len(predictions)
                        if overallProb[self.classVar.values[0]] > overallProb[self.classVar.values[1]]:
                            predicted = self.classVar[self.classVar.values.index(self.classVar.values[0])]
                        else:
                            predicted = self.classVar[self.classVar.values.index(self.classVar.values[1])]
                        if self.verbose: 
                            print str([cf.name for cf in self.classifiers])[1:-1].replace(",","\t")
                            print str(predictions)[1:-1].replace(",","\t")
                        probOf1 = overallProb[self.classVar.values[1]]
                        DFV = self.convert2DFV(probOf1)
                        probabilities = self.__getProbabilities(probOf1)
                else: # With Regression, use always the average fo the N Classifiers
                    self.status = "Using average of N classifiers (Regression)"
                    if self.verbose: print self.status
                    for c in self.classifiers:
                        predictions.append(c(origExample))
                    DFV = predicted = sum(predictions,0.0) / len(predictions)
                    probabilities = None 

                if resultType == orange.GetBoth:
                    if predicted:
                        orangePrediction = orange.Value(self.classVar, predicted)
                    else:
                        orangePrediction = None
                    res = orangePrediction, probabilities
                elif resultType == orange.GetProbabilities:
                    res = probabilities
                else: 
                    if predicted:
                        orangePrediction = orange.Value(self.classVar, predicted)
                    else:
                        orangePrediction = None
                    res = orangePrediction

                self.nPredictions += 1
                if returnDFV:
                    return (res,DFV)
                else:
                    return res
        else:
            # expression specified
            if origExample == None:
                return self.classifiers[0](None, resultType)
            else:
                # Predict using the models  
                predictions = {}   #The individual predictions for each repective classifier in classifiers
                predicted = None
                probabilities = None
                DFV = None
                if self.classVar.varType == orange.VarTypes.Discrete:
                    # Discrete expression
                    self.status = "Using supplied discrete expression (Discrete)"
                    if self.verbose:
                        print self.status

                    # Init votes
                    votes = {}
                    for pv in self.classVar.values:
                        votes[pv] = 0

                    # Do prediction and construct distribution of votes
                    for c in self.classifiers:
                        predictions[c] = self.classifiers[c](origExample)
                        votes[predictions[c].value] += 1

                    result = None
                    for exp in self.discreteExpression:
                        logicalExp, logicalRes = exp.split('->')
                        logicalExp = logicalExp.strip()
                        logicalRes = logicalRes.strip()
                        
                        if len(logicalExp) == 0:
                            predicted = logicalRes
                            break

                        rawParseTree = self._lexLogicalExp(logicalExp)
                        modParseTree = self._parseLogicalTree(rawParseTree, predictions, self.classVar.values)
                        result = self._interpretLogicalTree(modParseTree)
                        if result:
                            if self.verbose:
                                print "Logical Expression is True: ", ''.join(modParseTree)
                            predicted = logicalRes
                            break

                    self._isRealProb = True
                    #probOf1 = votes[self.classVar.values[1]]/len(self.classifiers)
                    #DFV = self.convert2DFV(probOf1)
                    #probabilities = self.__getProbabilities(probOf1)
                else:
                    self.status = "Using supplied regression expression (Regression)"
                    if self.verbose:
                        print self.status
                    
                    for c in self.classifiers:
                        predictions[c] = self.classifiers[c](origExample)

                    rawParseTree = self._lexRegressionExp(self.regressionExpression)
                    modParseTree = self._parseRegressionTree(rawParseTree, predictions)
                    result = self._interpretRegressionTree(modParseTree)
                
                    DFV = predicted = result
                    probabilities = None 
                    # Regression expression

                if resultType == orange.GetBoth:
                    if predicted:
                        orangePrediction = orange.Value(self.classVar, predicted)
                    else:
                        orangePrediction = None
                    
                    res = orangePrediction, probabilities
                
                elif resultType == orange.GetProbabilities:
                    res = probabilities
                else: 
                    if predicted:
                        orangePrediction = orange.Value(self.classVar, predicted)
                    else:
                        orangePrediction = None
                        
                    res = orangePrediction
                        
                    self.nPredictions += 1
                    if returnDFV:
                        return (res,DFV)
                    else:
                        return res
                
    def _parseRegressionTree(self, tree, predictionResults):
        """ Replace the variables with results from the classifiers """
        for k in predictionResults.keys():
            for n,i in enumerate(tree):
                if k == i:
                    tree[n] = str(predictionResults[k])
                    
        return tree
            
    def _interpretRegressionTree(self, tree):
        #print tree
        return eval(''.join(tree))

    def _lexRegressionExp(self, exp):
        STRING = 1
        exprList = list(token[STRING] for token
                        in generate_tokens(StringIO(exp).readline)
                        if token[STRING])
        return exprList

    def _parseLogicalTree(self, tree, predictionResults, predictionClasses):
        """ Replace the variables with results from the classifiers """
        for k in predictionResults.keys():
            for n,i in enumerate(tree):
                if k == i:
                    tree[n] = str(predictionResults[k])
                    
        """ Modify the tree such that predictionClasses are actual strings """
        for i, item in enumerate(tree):
            if item in predictionClasses:
                tree[i] = '"%s"' % tree[i]
                
        return tree
            
    def _interpretLogicalTree(self, tree):
        return eval(''.join(tree))

    def _lexLogicalExp(self, exp):
        STRING = 1
        exprList = re.split(r'[ ]| or | and ', exp)
        return exprList

    def convert2DFV(self,probOf1):
        # Subtract 0.5 so that the threshold is 0 and invert the signal as all learners have standard DFV:
        # Positive Values for the first element of the class attributes, and negatove values to the second
        DFV = -(probOf1-0.5)
        self._updateDFVExtremes(DFV)
        return DFV

    def __getProbabilities(self,ProbOf1):
        """Get the orange like output probabilities for the current predicted example"""
        #This is only valid for binary classifiers opencv limitation!
        # From openCv we know that the probability returned for this example represents the fraction of tree votes for class 1
        #Find the classValue string that is represented by the scalar 1 in opencvRF
        class1 = dataUtilities.CvMat2orangeResponse(1, self.classVar).value 
        dist = orange.DiscDistribution(self.classVar)
        dist[self.classVar.values.index(class1)] = ProbOf1
        dist[not self.classVar.values.index(class1)] = 1 - ProbOf1
        return dist


    def __generateProbabilities(self, prediction):
        # Method to artificialy generate a list the length of the number of classes and set the predicted class to 1
        dist = orange.DiscDistribution(self.classVar)
        dist[prediction]=1
        return dist

    def setDomainAndClass(self):
        #Find a possible Domain among the classifiers
        if not self.domain:
            for c in self.classifiers:
                if hasattr(c,"domain") and c.domain:
                    self.domain = c.domain
                    break
        #Find a possible classVar among the classifiers
        if not self.classVar:
            for c in self.classifiers:
                if hasattr(c,"classVar") and c.classVar:
                    self.classVar = c.classVar
                    break
        #If there wasn't found a classVar, try to get it from the Domain
        if not self.classVar and self.domain and self.domain.classVar:
            self.classVar = self.domain.classVar
 
        #If we were not success in getting a domain and a classVar, we cannot proceed!
        if not self.classVar or not self.domain:
            raise Exception("The classifiers are not compatible with the Consensus model. Try to use the respective Learners instead.")
            
                

    def write(self, dirPath):
        """ Save a Consensus model to disk including the domain used """
        if not self.classVar or not self.domain: self.setDomainAndClass() 
        try:
                #This removes any trailing '/'
                dirPath = os.path.realpath(str(dirPath))
                if os.path.isdir(dirPath):
                    modelFiles = glob.glob(os.path.join(dirPath,'C*.model'))
                    for Mfile in modelFiles:
                        os.system("rm -rf " + Mfile)     
                    os.system("rm -f " + os.path.join(dirPath,"trainDomain.tab"))
                    

                # This assures that all related files will be inside a folder
                os.system("mkdir -p " + dirPath) 
                # Save the models
                trainDomain = dataUtilities.DataTable(self.domain)
                #Save along with trainDomain file some dummy examples for compatibility
                ex = orange.Example(self.domain)
                for attr in self.domain:
                    if attr.varType == orange.VarTypes.Discrete:
                        ex[attr] = attr.values[0]
                    elif attr.varType == orange.VarTypes.Continuous:
                        ex[attr] = 0
                    elif attr.varType == orange.VarTypes.String:
                        ex[attr] = "NA"
                trainDomain.append(ex)
                trainDomain.save(os.path.join(dirPath,"trainDomain.tab"))
                for idx,c in enumerate(self.classifiers):
                    c.write(os.path.join(dirPath,"C"+str(idx)+".model"))
        except:            
                if self.verbose > 0: print "ERROR: Could not save the Consensus model to ", dirPath
                return False
        return True
        ##ecPA
             

def Consensusread(dirPath,verbose = 0):
    """Read a Consensus model from disk and return as a ConsensusClassifier instance. """
    # Read data from disk
    #This removes any trailing '/'
    dirPath = os.path.realpath(str(dirPath)) 
    basicStat = None
    NTrainEx = None
    imputeData = None
    # This assures that all related files will be inside a folder
    try:
        domainFile = dataUtilities.DataTable(os.path.join(dirPath,"trainDomain.tab"))

        #Load the models
        modelFiles = glob.glob(os.path.join(dirPath,'C*.model'))
        modelFiles.sort()
        if len(modelFiles) < 2:
                if verbose > 0: print "ERROR: Missing model files in ",dirPath    
                return None
        else:
                classifiers = []
                for mFile in modelFiles:
                    classifiers.append(AZBaseClasses.modelRead(mFile))
                    if not classifiers[-1]:
                        if verbose > 0: print "ERROR: Could not load the model ",mFile
                        return None
                    else:
                        #Try to load the imputeData, basicStat and NTrainEx from a model that saved it!
                        if hasattr(classifiers[-1], "basicStat") and classifiers[-1].basicStat and not basicStat:
                            basicStat = classifiers[-1].basicStat
                        if hasattr(classifiers[-1], "NTrainEx") and classifiers[-1].NTrainEx and not NTrainEx:
                            NTrainEx = classifiers[-1].NTrainEx
                        if hasattr(classifiers[-1], "imputeData") and classifiers[-1].imputeData and not imputeData:
                            imputeData = classifiers[-1].imputeData
                            domainFile = imputeData #This is needed for domain compatibilitu betwene imputer and domain var

    except:
        if verbose > 0: print "ERROR: It was not possible to load the Consensus model"
        return None
    return ConsensusClassifier(classifiers=classifiers ,varNames = [attr.name for attr in domainFile.domain.attributes],classVar = domainFile.domain.classVar, verbose = verbose, domain = domainFile.domain, basicStat = basicStat, NTrainEx = NTrainEx, imputeData = imputeData)


if __name__ == "__main__":
   
    Consensuslearner = ConsensusLearner()        
    
    # Read data set of orange type
    fileName = "train.tab"
    data = dataUtilities.DataTable(fileName)
    
    # Create a Consensus model from data
    ConsensusClassifier = Ccclearner(data)


