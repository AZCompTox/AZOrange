"""
AZorngRF
Module for interfacing the random forest models in OpenCV with AZOrange. 
The OpenCV RF models should be usable just like any other learner or classifier in the Orange package. 
"""
import string
import os
import pickle
#import time

##scPA
import AZBaseClasses
from math import sqrt
from opencv import ml
from opencv import cv

##ecPA
import orange

from AZutilities import dataUtilities
from AZutilities import miscUtilities

import AZOrangeConfig as AZOC

##scPA Now the RF derivates from base class AZBaseClasses.AZLearner  ##ecPA
class RFLearner(AZBaseClasses.AZLearner):
    """
    Creates an RF as an Orange type of learner instance. 
    """

    def __new__(cls, trainingData = None, name = "RF learner", **kwds):
        self = AZBaseClasses.AZLearner.__new__(cls, **kwds)
        if trainingData:
            self.__init__(name, **kwds)
            return self.__call__(trainingData)
        else:
##scPA
            self.__dict__.update(kwds)
            self.name = name
##ecPA
            return self

    def __init__(self, name = "RF learner", **kwds):
        """
        Set default values of the model parameters if they are not given as inputs
        """
        self.learner = None
        self.name = name
        self.priors = None
        ##scPA
        self.imputer = None
        self.verbose = 0 
        ##ecPA
        # Append arguments to the __dict__ member variable
        self.__dict__.update(kwds)

        # useBuiltInMissValHandling
        if not self.__dict__.has_key("useBuiltInMissValHandling"):
            self.__dict__["useBuiltInMissValHandling"] = AZOC.RFDEFAULTDICT["useBuiltInMissValHandling"]

        # NumThreads: Number of threads to be used by opencv. 
        # 0 means it will use as much threads as number of cores. Default 1.
        if not self.__dict__.has_key("NumThreads"):
            self.__dict__["NumThreads"] = AZOC.RFDEFAULTDICT["NumThreads"]



        # maxDepth
        if not self.__dict__.has_key("maxDepth"): 
            self.__dict__["maxDepth"] = AZOC.RFDEFAULTDICT["maxDepth"]

        # minSample
        if not self.__dict__.has_key("minSample"): 
            self.__dict__["minSample"] = AZOC.RFDEFAULTDICT["minSample"]

        # useSurrogates
        if not self.__dict__.has_key("useSurrogates"): 
            self.__dict__["useSurrogates"] = AZOC.RFDEFAULTDICT["useSurrogates"]

        # getVarVariance 
        if not self.__dict__.has_key("getVarVariance"): 
            self.__dict__["getVarVariance"] = AZOC.RFDEFAULTDICT["getVarVariance"]

        # nActVars
        if not self.__dict__.has_key("nActVars"): 
            self.__dict__["nActVars"] = AZOC.RFDEFAULTDICT["nActVars"]

        # nTrees
        if not self.__dict__.has_key("nTrees"): 
            self.__dict__["nTrees"] = AZOC.RFDEFAULTDICT["nTrees"]

        # forestAcc 
        if not self.__dict__.has_key("forestAcc"): 
            self.__dict__["forestAcc"] = AZOC.RFDEFAULTDICT["forestAcc"]

        # termCrit = 0 -> number of trees, termCrit = 1 -> oob error = forestAcc used to terminate training. 
        if not self.__dict__.has_key("termCrit"): 
            self.__dict__["termCrit"] = AZOC.RFDEFAULTDICT["termCrit"]

        # Stratified data sampling. 
        if not self.__dict__.has_key("stratify"): 
            self.__dict__["stratify"] = AZOC.RFDEFAULTDICT["stratify"]


    def __call__(self, trainingData, weight = None):
        """Creates an RF model from the data in trainingData. """
        if not AZBaseClasses.AZLearner.__call__(self,trainingData, weight):
            return None

        # Set the number of theatd to be used ny opencv
        cv.cvSetNumThreads(max(int(self.NumThreads),0))
        #Remove from the domain any unused values of discrete attributes including class
        trainingData = dataUtilities.getDataWithoutUnusedValues(trainingData,True)

        # Object holding the data req for predictions (model, domain, etc)
	#print time.asctime(), "=superRFmodel(trainingData.domain)"
        ##scPA
        # Remove meta attributes from training data
        #dataUtilities.rmAllMeta(trainingData)
        if len(trainingData.domain.getmetas()) == 0:
            trainData = trainingData
        else:
            trainData = dataUtilities.getCopyWithoutMeta(trainingData)
        # Impute the data and Convert the ExampleTable to CvMat 
        if self.useBuiltInMissValHandling:
            #Create the imputer empty since we will not be using it
            impData = dataUtilities.DataTable(trainData.domain)
            CvMatrices = dataUtilities.ExampleTable2CvMat(trainData)
        else:
            #Create the imputer
            self.imputer = orange.ImputerConstructor_average(trainData)
            impData=self.imputer.defaults
            trainData = self.imputer(trainData)
            CvMatrices = dataUtilities.ExampleTable2CvMat(trainData)
            CvMatrices["missing_data_mask"] = None
        ##ecPA
        self.learner = ml.CvRTrees()#superRFmodel(trainData.domain)    #This call creates a scratchDir

        # Set RF model parameter values
        #  when nActVars defined as 0, use the sqrt of number of attributes so the user knows what will be used
        # This would be done in the C level if left as 0
        if self.nActVars == "0" and len(trainData.domain.attributes)>0:
            self.nActVars =  str(int(sqrt(len(trainData.domain.attributes))))
	#print time.asctime(), "=self.setParameters"
        params = self.setParameters(trainData)
        # Print values of the parameters
        if self.verbose > 0: self.printOuts(params)
        #**************************************************************************************************//
        #                      Check for irrational input arguments
        #**************************************************************************************************//
        if params.min_sample_count >= len(trainingData):
            if self.verbose > 0: print "ERROR! Invalid minSample: ",params.min_sample_count
            if self.verbose > 0: print "minSample must be smaller than the number of examples."
            if self.verbose > 0: print "The number of examples is: ",len(trainingData)
            if len(trainingData) > 10:
                if self.verbose > 0: print "minSample assigned to default value: 10"
                params.min_sample_count = 10
            else:
                if self.verbose > 0: print "Too few examples!!"
                if self.verbose > 0: print "Terminating"
                if self.verbose > 0: print "No random forest model built"
                return None
        if params.nactive_vars > len(trainingData.domain.attributes):
            if self.verbose > 0: print "ERROR! Invalid nActVars: ",params.nactive_vars
            if self.verbose > 0: print "nActVars must be smaller than or equal to the number of variables."
            if self.verbose > 0: print "The number of variables is: ", len(trainingData.domain.attributes)
            if self.verbose > 0: print "nActVars assigned to default value: sqrt(nVars)=",sqrt(len(trainingData.domain.attributes))
            params.nactive_vars = 0;
        # Train RF model on data in openCVFile
	#print time.asctime(), "=Start Training"
        #Process the priors and Count the number of values in class var
        if  trainingData.domain.classVar.varType == orange.VarTypes.Discrete:
            cls_count = len(trainData.domain.classVar.values)
            priors = self.convertPriors(self.priors,trainingData.domain.classVar)
            if type(priors) == str: #If a string is returned, there was a failure, and it is the respective error mnessage.
                print priors
                return None 
        else:
            cls_count = 0
            priors = None
        # Call the train method
        self.learner.train( CvMatrices["matrix"],ml.CV_ROW_SAMPLE,CvMatrices["responses"],None,None,CvMatrices["varTypes"],CvMatrices["missing_data_mask"],params,cls_count,  priors and str(priors).replace(","," ") or None)
        if self.learner.get_var_importance():
            varImportanceList = self.learner.get_var_importance()
            varImportance = {}
            varName = []
            varImp = []
            for idx,attr in enumerate(CvMatrices["varNames"]):
                varImportance[attr] = varImportanceList[idx]
            #Uncomment next lines if needed the outpuit already ordered
            #============================= begin =================================
            #    varName.append(attr)
            #    varImp.append(varImportanceList[idx])
            #Order the vars in terms of importance
            # insertion sort algorithm
            #for i in range(1, len(varImp)):
            #    save = varImp[i]
            #    saveName = varName[i]
            #    j = i
            #    while j > 0 and varImp[j - 1] < save:
            #        varImp[j] = varImp[j - 1]
            #        varName[j] = varName[j - 1]
            #        j -= 1
            #    varImp[j] = save
            #    varName[j] = saveName
            #For debug: test if assign var importance was correct
            #for attr in varImportance:
            #    if varImportance[attr] != varImp[varName.index(attr)]:
            #        print "ERROR: Variable importance of ", attr, " is not correct!"
            #OrderedVarImportance = {"VarNames":varName, "VarImportance":varImp}
            #=============================  end  =================================
        else:
            varImportance = {}
        #print time.asctime(), "=Done"
        # Save info about the variables used in the model (used by the write method)
        #attributeInfo = dataUtilities.DataTable(trainData.domain)
        # place the impute data as the first example of this data
        #attributeInfo.append(self.imputer.defaults)
        return RFClassifier(classifier = self.learner, classVar = impData.domain.classVar, imputeData=impData, verbose = self.verbose, varNames = CvMatrices["varNames"],thisVer=True,useBuiltInMissValHandling = self.useBuiltInMissValHandling, varImportance = varImportance, basicStat = self.basicStat, NTrainEx = len(trainingData), parameters = self.parameters)


    def setParameters(self, trainingData):

        # Get all parameters for the RF in OpenCV
        self.nVars = str(len(trainingData.domain.attributes))
        self.nExamples = str(len(trainingData))
        self.maxDepth = self.__dict__["maxDepth"]
        self.minSample = self.__dict__["minSample"]
        self.useSurrogates = self.__dict__["useSurrogates"]
        self.getVarVariance = self.__dict__["getVarVariance"]
        self.nActVars = self.__dict__["nActVars"]
        self.nTrees = self.__dict__["nTrees"]
        self.forestAcc = self.__dict__["forestAcc"]
        self.termCrit = self.__dict__["termCrit"]


        #Configure RTrees params
        #ml.CvRTParams(max_depth, min_sample_count,  regression_accuracy, use_surrogates,   max_categories,  priors, calc_var_importance, __native_vars, __max_tree_count,  forestaccuracy, term_crit.type   )
        #params = ml.CvRTParams(20,  5,  0.0   ,False, 15,  None, False, 0,   50, 0.1, 1)

        params = ml.CvRTParams()
        if str(self.getVarVariance).lower() == "false":
            params.calc_var_importance = False
        else:
            params.calc_var_importance = True
        params.max_categories = 15
        params.max_depth = int(self.maxDepth)
        params.min_sample_count = int(self.minSample)
        #params.priors = None  # The priors are set in other way. 
        params.regression_accuracy = 0.0
        if str(self.useSurrogates).lower() == "false":
            params.use_surrogates = False
        else:
            params.use_surrogates = True
        params.nactive_vars = int(self.nActVars)
        #params.truncate_pruned_tree = False

        term_crit = cv.CvTermCriteria()
         #cv.CV_TERMCRIT_EPS (1 on Canvas)  #  or cv.CV_TERMCRIT_ITER (0 on Canvas)
        if int(self.termCrit) == 0:                      # CV_TERMCRIT_ITER
            term_crit.type = cv.CV_TERMCRIT_ITER
        else:                                           # CV_TERMCRIT_EPS
            term_crit.type = cv.CV_TERMCRIT_EPS                 
        term_crit.epsilon = float(self.forestAcc)     # OOB error
        term_crit.max_iter = int(self.nTrees)         # max_Tree_Count

        params.term_crit =  term_crit

        return params

    def printOuts(self,params):

        print("Random Forest model trained with the following parameters: ")
        print("maxDepth: "+str(params.max_depth))
        print("minSample: "+str(params.min_sample_count))
        print("useSurrogates: "+str(params.use_surrogates))
        print("getVarVariance: "+str(params.calc_var_importance))
        print("nActVars: "+str(params.nactive_vars))
        print("nTrees: "+str(params.term_crit.max_iter))
        print("forestAcc: "+str(params.term_crit.epsilon))
        if params.term_crit.type == cv.CV_TERMCRIT_ITER:
            print("termCrit: CV_TERMCRIT_ITER (NTrees)")
        elif params.term_crit.type == cv.CV_TERMCRIT_EPS:
            print("termCrit: CV_TERMCRIT_EPS (OOC error)")
        else:
            print("termCrit: unknown")

        print("maxCategories: "+str(params.max_categories))
        print("truncatePrunedTree: "+str(params.truncate_pruned_tree))
        print("priors: "+str(params.priors))
        print("regressionAccuracy: "+str(params.regression_accuracy))


##scPA Now the RF classifier derivates from base class AZBaseClasses.AZClassifier ##ecPA
class RFClassifier(AZBaseClasses.AZClassifier):
    def __new__(cls, name = "RF classifier", **kwds):
        self = AZBaseClasses.AZClassifier.__new__(cls, name = name,  **kwds)        
        #self.__init__(name, **kwds)
        return self
##scPA
    def __init__(self, name = "RF classifier", **kwds):
        self.verbose = 0
        self.varImportance =  {}
        self.__dict__.update(kwds)
        self._isRealProb = False
        self.name = name
        self.domain = None
        self.ExFix = dataUtilities.ExFix()
        if self.imputeData != None:
            '''Create the imputer: the imputer needs the imputeData to exists allong it's life time'''
            try :
                self.domain = self.imputeData.domain
                if self.useBuiltInMissValHandling:
                    self.imputer = None
                else:
                    self.imputer = orange.Imputer_defaults(self.imputeData)
            except:
                self.imputer = None
                if self.verbose > 0: print "Unable to create the imputer, or even the builtIn RF imputer."
                return None
        else:
            self.imputer = None
            if self.verbose > 0: print "Warning! - No impute data defined"
            return None
##ecPA
    def getVarImportance(self, var, data = None):
        """ retunrs the var variable importance 
            The data parameter is ignored. It'f just for compatibility with other estimators
        """
        if not self.varImportance:
            return None
        else:
            if type(var) == str:
                return self.varImportance[var]
            else:
                return self.varImportance[var.name] 


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
        if origExample == None:
            return self.classifier(None, resultType)
        else:
##scPA
            # Remove Meta attributes from example
            #dataUtilities.rmAllMeta(example)
            if len(origExample.domain.getmetas()) == 0:
                example = origExample
            else:
                example = dataUtilities.getCopyWithoutMeta(origExample)

            if not self.ExFix.ready:
                self.ExFix.set_domain(self.domain)
                self.ExFix.set_examplesFixedLog(self.examplesFixedLog)
            inExample = self.ExFix.fixExample(example) 
            if inExample: ##only procceds if the example was fixed or already ok, i.e.   inExample !=  None
##ecPA          
                ##scPA
                if  self.useBuiltInMissValHandling:
                    #compute the missing _mask
                    (exampleCvMat, missing_mask) = dataUtilities.Example2CvMat(inExample,self.varNames,self.thisVer,True) 
                else:
                    missing_mask = None
                    if self.imputer:
                         examplesImp = self.imputer(inExample)
                         # There is a problem with using the imputer when examples contain meta attributes.
                         # Unable to remove meta attributes from the examples. OK to rm meta from ExampleTables, but not from Example objects.
                         if not examplesImp:
                             if self.verbose > 0: print "Unable to predict with the RF model."
                             if self.verbose > 0: print "Perhaps you need to remove meta attributes from your examples."
                             return None
                    else:
                         examplesImp = inExample
                    ##ecPA
                    # Remove the response variable from the example to be predicted and transfrom the example to a tab sep string 
                    exampleCvMat = dataUtilities.Example2CvMat(examplesImp,self.varNames,self.thisVer)
                    del examplesImp
                if not exampleCvMat:
                    if self.verbose > 0: print "Could not convert the example to a valid CvMat objct for prediction"
                    return none
                # Predict using the RFmodel object
                prediction = self.classifier.predict(exampleCvMat,missing_mask)
	        probabilities = None
                DFV = None
                # Back transform the prediction to the original classes and calc probabilities
                prediction = dataUtilities.CvMat2orangeResponse(prediction, self.classVar)
                # Calculate artificial probabilities - not returned by the OpenCV RF algorithm
                if self.classVar.varType == orange.VarTypes.Discrete:
                    if resultType != orange.GetValue:
                        if len(self.classVar.values) == 2:
                            probOf1 = self.classifier.predict_prob(exampleCvMat,missing_mask)
                            probabilities = self.__getProbabilities(probOf1)
                            DFV = self.convert2DFV(probOf1)
                            self._isRealProb = True 
                        else:
                            #Need to make sure to return meanful probabilities to the cases where opencvRF does not support probabilities
                            # to be compatible with possible callers asking for probabilities. 
                            probabilities = self.__generateProbabilities(prediction)
                            self._isRealProb = False
                    elif len(self.classVar.values) == 2 and returnDFV:
                        DFV = self.convert2DFV(self.classifier.predict_prob(exampleCvMat,missing_mask))
                else:
                    #On Regression models assume the DVF as the value predicted
                    DFV = prediction
                    self._updateDFVExtremes(DFV)
                del exampleCvMat
                del inExample
            else:
                if self.verbose > 0:
                    print "No prediction made for example:"
                    print example
                    print "The example does not have the same variables as the model."
                prediction = None
                probabilities = None
                DFV = None
	    ##scPA
	    del example
	    ##ecPA
            if resultType == orange.GetBoth:
                if prediction: 
                    orangePrediction = orange.Value(self.classVar, prediction)
                else: 
                    orangePrediction = None
                res = orangePrediction, probabilities
            elif resultType == orange.GetProbabilities:
                res = probabilities
            else: 
                if prediction: 
                    orangePrediction = orange.Value(self.classVar, prediction)
                else: 
                    orangePrediction = None
                res = orangePrediction 

            self.nPredictions += 1
            if returnDFV:
                return (res,DFV)
            else:
                return res

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

        #probabilities = [0 , 0]
        #Set the classValue corresponding to 1  to the probability returned from openCv
        #probabilities[self.classVar.values.index(class1)] = ProbOf1
        #Set the other classValue to the remain of probability  
        #probabilities[not self.classVar.values.index(class1)] = 1 - ProbOf1
        #return probabilities 


    def __generateProbabilities(self, prediction):
        # Method to artificialy generate a list the length of the number of classes and set the predicted class to 1
        dist = orange.DiscDistribution(self.classVar)
        dist[prediction]=1
        return dist
        #probabilities = [0] * len(self.classVar.values)
        #probabilities[self.classVar.values.index(prediction.value)] = 1
        #return probabilities


    def write(self, dirPath):
        """Save a RF model to disk with the data used to train the model.
           It is imparative that the model is saved with the data used for training. Only the domain is used. """
         
        try:
                #This removes any trailing '/'
                dirPath = os.path.realpath(str(dirPath))

                # This assures that all related files will be inside a folder
                os.system("mkdir -p " + dirPath) 
                
                filePath = os.path.join(dirPath,"model.rf")

                # The impute Data was previously added to the self.attributeInfo
                # Remove the meta attributes from the imputer data. We don't need to store them along with the model

                if  self.useBuiltInMissValHandling:
                    impData = self.imputeData
                else:
                    # Save a data set with one row containing the impute values
                    impData = dataUtilities.DataTable(self.imputer.defaults.domain)
                    impData.append(self.imputer.defaults)
                # Remove the meta attributes from the imputer data. We don't need to store them along with the model
                impData = dataUtilities.getCopyWithoutMeta(impData)


                impData.save(os.path.join(dirPath,"ImputeData.tab"))


                #Save the info about the train data as:
                #    var names ordered the same way the Learner was trained
                #    NTrainEx
                #    basicStat 

                varNamesFile = open(os.path.join(dirPath,"varNames.txt"),"w")

                varNamesFile.write(str(self.varNames)+"\n") 
                varNamesFile.write(str(self.NTrainEx)+"\n") 
                varNamesFile.write(str(self.basicStat)+"\n") 
                varNamesFile.close()
                #Save the parameters
                self._saveParameters(os.path.join(dirPath,"parameters.pkl"))
                # Save the model
                self.classifier.save(filePath)
        except:            
                if self.verbose > 0: print "ERROR: Could not save model to ", path
                return False
        return True
        ##ecPA
             

def RFread(dirPath,verbose = 0):
    """Read a RF model from disk and return as a RFClassifier instance. """
    # Read data from disk
    ##scPA   
    #This removes any trailing '/'
    dirPath = os.path.realpath(str(dirPath)) 
    NTrainEx = 0
    basicStat = None 
    # This assures that all related files will be inside a folder
    loadedRFclassifier = ml.CvRTrees()
    # Read the parameters
    if os.path.isfile(os.path.join(dirPath,"parameters.pkl")):
        fileh = open(os.path.join(dirPath,"parameters.pkl"),"r")
        parameters = pickle.load(fileh)
        fileh.close()
    else:
        parameters = {}
    if os.path.isfile(os.path.join(dirPath,"model.rf")):
        #New format
        filePath = os.path.join(dirPath,"model.rf") 
        impDataPath = os.path.join(dirPath,"ImputeData.tab")
        varNamesPath = os.path.join(dirPath,"varNames.txt")
        loadedRFclassifier.load(filePath)
    else:
        #Old Format
        # For RF models we assume that the model file has the same name as the model folder
        #filePath = os.path.join(dirPath,os.path.split(dirPath)[1])
        #root, ext = os.path.splitext(filePath)
        files = os.listdir(dirPath)
        filePath = None
        impDataPath = None
        varNamesPath = None
        # Load the model when found
        for file in files:
            if len(file) >= 9 and file[-9:] == "Saved.tab":
                impDataPath = os.path.join(dirPath,file)
            elif  len(file) >= 12 and file[-12:] == "varNames.txt":
                varNamesPath = os.path.join(dirPath,file)
            elif filePath is None:
                # looking for opencv-ml-random-trees  in first 10 lines
                fh = open(os.path.join(dirPath,file),'r')
                for i in range(10):
                    if "opencv-ml-random-trees" in fh.readline():
                        filePath = os.path.join(dirPath,file)
                        break
                fh.close()
                if filePath:
                    try:
                        loadedRFclassifier.load(filePath)
                    except:
                        filePath = None
        if not filePath or not impDataPath or not varNamesPath:
            print "Error loading RF model: Missing files. Files found:",files
            return None
    ##scPA 
    try:
        impData = dataUtilities.DataTable(impDataPath,createNewOn=orange.Variable.MakeStatus.OK)
        classVar = impData.domain.classVar

        #Load the var names oredered the way it was used when training
        if (os.path.isfile(varNamesPath)):
            if len(impData) == 0:
                useBuiltInMissValHandling = True
            else:
                useBuiltInMissValHandling = False
                impData = impData[0]
            varNamesFile = open(varNamesPath,"r")
            lines = varNamesFile.readlines()
            varNames = eval(lines[0].strip())
            if len(lines) >= 3:
                NTrainEx = eval(lines[1].strip())
                basicStat = eval(lines[2].strip())
            varNamesFile.close()
            thisVer = True
        else:
            useBuiltInMissValHandling = False
            if verbose > 0: print "WARNING: The model loaded was probably saved with azorange version 0.2.1 or lower"
            varNames = [attr.name for attr in impData.domain.attributes]
            thisVer = False
    except:
        if verbose > 0: print "ERROR: It was not possible to load the impute data or the varNames."
        return None
    ##ecPA   also added , imputeData=impData to nexti call 
    return RFClassifier(classifier = loadedRFclassifier, classVar = classVar, imputeData=impData, verbose = verbose, varNames = varNames, thisVer=thisVer, useBuiltInMissValHandling = useBuiltInMissValHandling, NTrainEx = NTrainEx, basicStat = basicStat, parameters = parameters)


if __name__ == "__main__":
   
    RFlearner = RFLearner()        
    
    # Read data set of orange type
    fileName = "train.tab"
    data = dataUtilities.DataTable(fileName)
    
    # Create a RF model from data
    RFClassifier = RFlearner(data)


