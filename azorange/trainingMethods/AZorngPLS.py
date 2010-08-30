"""
AZorngPLS
Module for the incoperation of the PLearn-PLS package into Orange. 
The aim is to be able to use PLS just like any other learning algorithm in Orange.
"""
import string
import os
import time

import pls
import orange
import orngTest
import orngBayes
import orngStat
import random
import time
import types

import AZOrangeConfig as AZOC
import AZBaseClasses
from AZutilities import dataUtilities
from AZutilities import miscUtilities

class PLSLearner(AZBaseClasses.AZLearner):

    def __new__(cls, trainingData = None, name = "PLS learner", **kwds):
        self = AZBaseClasses.AZLearner.__new__(cls, **kwds)
        if trainingData:
            self.__init__(name, **kwds)
            return self.__call__(trainingData)
        else:
	    # Append arguments to the __dict__ member variable
            self.__dict__.update(kwds)
	    self.name = name
            return self
   
    #def setattr(self, name, value):
        #self.__dict__[name] = value
	    
    def __init__(self, name = "PLS learner", **kwds):
        """
        Set default values of the PLS parameters if they are not given as inputs.
        """

	#Set verbosity flag
        self.verbose = 0

	self.imputer = None
        self.learner = None

	# Append arguments to the __dict__ member variable
        self.__dict__.update(kwds)
	self.name = name

	# Precision
	if not self.__dict__.has_key("precision"):
		self.__dict__["precision"] = AZOC.PLSDEFAULTDICT["precision"] 
	# Method
	if not self.__dict__.has_key("method"):
                self.__dict__["method"] = AZOC.PLSDEFAULTDICT["method"]
        
	# k - Number of components
	if not self.__dict__.has_key("k"):
		self.__dict__["k"] = AZOC.PLSDEFAULTDICT["k"]


    def __call__(self, trainingData, weight=None):
        """Creates an PLS model from the data in trainingData. """
        if not AZBaseClasses.AZLearner.__call__(self,trainingData, weight):
            return None
        #Remove from the domain any unused values of discrete attributes including class
        trainingData = dataUtilities.getDataWithoutUnusedValues(trainingData,True)
        # Create path for the Orange data
        scratchdir = miscUtilities.createScratchDir(desc="PLS")
        OrngFile = os.path.join(scratchdir,"OrngData.tab")

        # Remove meta attributes from training data to make the imputer work with examples without the meta attributes. 
        #dataUtilities.rmAllMeta(trainingData)
        if len(trainingData.domain.getmetas()) == 0:
            trainData = trainingData
        else:
            trainData = dataUtilities.getCopyWithoutMeta(trainingData)

	# Create the imputer
        self.imputer = orange.ImputerConstructor_average(trainData)
	# Impute the data 
	trainData = self.imputer(trainData)
        # Save the Data already imputed to an Orange formated file
	if self.verbose > 1: print time.asctime(), "Saving Orange Data to a tab file..."
        orange.saveTabDelimited(OrngFile,trainData)
	if self.verbose > 1: print time.asctime(), "done"

        # Create the PLS instance
	if self.verbose > 1: print time.asctime(), "Creating PLS Object..."
        learner = pls.PlsAPI()
	if self.verbose > 1: print time.asctime(), "done"

	# Assign the PLS parameters
	learner.SetParameter('v',str(self.verbose))
        learner.SetParameter('debug',str(int(self.verbose > 0)))
	learner.SetParameter('method',self.method)
        if types.IntType(self.k) > len(trainData.domain.attributes):
	    learner.SetParameter('k',str(len(trainData.domain.attributes)))
            if self.verbose > 0: print "Warning! The number of components were more than the number of attributes."
            if self.verbose > 0: print "   Components were set to ",len(trainData.domain.attributes)
        else:
	    learner.SetParameter('k',self.k)
	learner.SetParameter('precision',self.precision)	
	learner.SetParameter('sDir',AZOC.SCRATCHDIR)
	
        # Read the Orange Formated file and Train the Algorithm
	# TRAIN
	if self.verbose > 1: print time.asctime(), "Training..."
        learner.Train(OrngFile)
	if self.verbose > 1:
		print "Train finished at ", time.asctime()
		print "PLS trained in: " + str(learner.GetCPUTrainTime()) + " seconds";
		print "Method:     " +  learner.GetParameter("method")
		print "Components: " +  learner.GetParameter("k")
		print "Precision:  " +  learner.GetParameter("precision")

        # Remove the scratch file
        if self.verbose == 0:
	    miscUtilities.removeDir(scratchdir)
	else:
	    print "The directory " + scratchdir + " was not deleted because DEBUG flag is ON"
	del trainData
        impData=self.imputer.defaults
        return PLSClassifier(classifier = learner, name = "Classifier of " + self.name, classVar = trainingData.domain.classVar, imputeData=impData, verbose = self.verbose)#learner.GetClassVarName())#


class PLSClassifier(AZBaseClasses.AZClassifier):
    def __new__(cls, name = "PLS classifier", **kwds):
        self = AZBaseClasses.AZClassifier.__new__(cls, name = name,  **kwds)
        #self.__init__(name, **kwds)
	return self

    def __init__(self, name = "PLS classifier", **kwds): 
        self.verbose = 0
        self.__dict__.update(kwds)
	self.name = name
        self.imputer = None
        self.domain = None
        self.ExFix = dataUtilities.ExFix()

        if self.imputeData:
            '''Create the imputer: the imputer needs the imputeData to exists allong it's life time'''
            try :
                self.domain = self.imputeData.domain
                self.imputer = orange.Imputer_defaults(self.imputeData)
            except:
            	self.imputer = None
            	if self.verbose > 0: print "Unable to create the imputer"
	else:
	    if self.verbose > 0: print "Warning! - No impute data defined"


    def __call__(self, origExamples = None, resultType = orange.GetValue):
        if origExamples == None:
            return self.classifier(None, resultType)
        else:
            if len(origExamples.domain.getmetas()) == 0:
                examples = origExamples
            else:
                examples = dataUtilities.getCopyWithoutMeta(origExamples)
            #dataUtilities.rmAllMeta(examples) 
           
            #Check if the examples are compatible with the classifier (attributes order and varType compatibility)
            if self.imputer:
                dataUtilities.verbose = self.verbose
                if not self.ExFix.ready:
                    self.ExFix.set_domain(self.imputer.defaults.domain)
                    self.ExFix.set_examplesFixedLog(self.examplesFixedLog)
                inExamples = self.ExFix.fixExample(examples)
            else:
                inExamples=None
	    if not inExamples:
                if self.verbose > 1:
        	    print "No prediction made for example:"
                    try:
                        print str(examples)[0:str(examples).find(",",20)]+" ... "+str(examples)[str(examples).rfind(",")+1:]
                    except:
                        print examples
                    print "The example does not have the same variables as the model, or the varTypes are incompatible."
		return None

	    #Imput the examples if there are missing values	
	    if self.imputer:  
		examplesImp = self.imputer(inExamples)         
                # There is a problem with using the imputer when examples contain meta attributes.
                # Unable to remove meta attributes from the examples. OK to rm meta from ExampleTables, but not from Example objects.
                if not examplesImp:
                    if self.verbose > 0: print "Unable to predict with the PLS model."
                    if self.verbose > 0: print "Perhaps you need to remove meta attributes from your examples."
                    return None
	    else:
	        examplesImp = inExamples
		
            # Transform the orange data to the PLS prediction data format 
            PLSFeatureVector = self.getFeatureVector(examplesImp)
            # Return the result of the prediction for one feature vector
	    PLSOut = self.classifier.Run(PLSFeatureVector)
	    if PLSOut.find("ERROR")>=0:
                if self.verbose > 0: print PLSOut
                PLSOut = '?' #"ERROR"
	    orngOut=string.split(PLSOut,"\t")
	    #convert result to orange value
            try:
	        value=orange.Value(self.classVar,orngOut[len(orngOut)-1])
            except:
                print "Error converting the Class back to orange format:"
                print "Class:",str(self.classVar)  
                if self.classVar.varType == orange.VarTypes.Discrete:
                    print "values = ",str(self.classVar.values)
                else:
                    print "Numerical Variable"
                print "Returned by PLS:",str(PLSOut)
                print "Value in orange Format (Will be the last element of PLSout):",str(orngOut[len(orngOut)-1])
	    score = self.getProbabilities(value)
	    # Assure that large local variables are deleted
	    del examplesImp
	    del PLSFeatureVector

	    #Return the desired quantity	
            if resultType == orange.GetProbabilities:
		return score
	    else:
	 	if resultType == orange.GetBoth:
			return value, score
		else:
            		return value

    def getProbabilities(self, prediction):
	try:
            prob = [0.0]*len(self.classVar.values)
	    prob[self.classVar.values.index(str(prediction))]=1.0
	except:
	    return [0.0]
	return prob

    def getFeatureVector(self, orangeVector):
        """ Transforms one orange type of data example ([5.1, 3.5, POS, 0.2, 'Iris-setosa'], special list) 
            to a PLS type ("5.1 3.5 POS 0.2", string). """
        PLSStrVector = ""
        # The PLS expects the classification variable to be the last column
	# The PLS will only use the first n columns needed and ignore the last one(s), but just in case...
        if orangeVector.domain.classVar:
            orangeVector[orangeVector.domain.classVar] = '?'
        for idx in range(len(orangeVector)):
		if idx!=0: PLSStrVector+=" "
		PLSStrVector+=(str(orangeVector[idx]))
        return PLSStrVector

    def write(self, filePath, data=None):
        """Save a PLS classifier to disk"""
        try: 
                # Save classifier
                self.classifier.SavePLSModel(str(filePath))
                if not self.imputer:
                    if self.verbose > 0: print "ERROR: PLS model saved without impute data"
                    return False 
                # Save a data set with one row containing the impute values
                impData = dataUtilities.DataTable(self.imputer.defaults.domain)
                impData.append(self.imputer.defaults)
                # Remove the meta attributes from the imputer data. We don't need to store them along with the model
                impData = dataUtilities.getCopyWithoutMeta(impData)
                impData.save(str(filePath)+"/ImputeData.tab")
        except:            
                if self.verbose > 0: print "ERROR: Could not save model to ", path
                return False
        return True

def PLSread(filePath, verbose = 0):
    """Read a PLS classifier from disk and return as a PLSClassifier instance. """
    # Create a PLS instance
    PLSInstance=pls.PlsAPI()
    # Read PLS model from disk
    if PLSInstance.LoadPLSModel(str(filePath))!=1:
	if verbose > 0: print "PLS model not loaded!"
	return None
    
    try:
        impData = dataUtilities.DataTable(str(filePath)+"/ImputeData.tab",createNewOn=orange.Variable.MakeStatus.OK)
        classVar = impData.domain.classVar
    except:
	if verbose > 0: print "ERROR: It was not possible to load the impute data"
	return 
	 
    return PLSClassifier(classifier = PLSInstance, classVar=classVar, imputeData=impData[0], verbose = verbose)


