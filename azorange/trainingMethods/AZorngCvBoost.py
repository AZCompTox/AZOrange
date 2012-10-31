import string,orange
import AZBaseClasses
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC
import os
from opencv import ml,cv
import pickle

class CvBoostLearner(AZBaseClasses.AZLearner):
    """
    Creates a opencv Boost learner derivated from AZBaseClasses.AZLearner. 
    """
    def __new__(cls, trainingData = None, name = "CvBoost learner", **kwds):
        self = AZBaseClasses.AZLearner.__new__(cls, **kwds)
        if trainingData:
            self.__init__(name, **kwds)
            return self.__call__(trainingData)
        else:
            self.__dict__.update(kwds)
            self.name = name
            return self

    def __init__(self, name = "CvBoost learner", **kwds):
        self.verbose = 0
        self.name = name
        self.trainData = None
        self.imputer = None
        self.compatibility = ["classification"]
	#Read default parameters from AZOrangeConfig.py file
        #CVBOOSTTYPE = { "DISCRETE":0, "REAL":1, "LOGIT":2, "GENTLE":3 }
        #CVBOOSTSPLITCRIT{ "DEFAULT":0, "GINI":1, "MISCLASS":3, "SQERR":4 }
        #CVBOOSTDEFAULTDICT = {"boost_type":"DISCRETE","weak_count":100,"split_criteria":"DEFAULT","weight_trim_rate":0.95, "max_depth":1, "use_surrogates":True, "priors":None}

        for par in ("boost_type","weak_count","split_criteria","weight_trim_rate", "max_depth", "use_surrogates","priors"): 
            setattr(self, par, AZOC.CVBOOSTDEFAULTDICT[par])
        self.__dict__.update(kwds)

    def isCompatible(self, classVar):
        """Checks if the learner is compatiblue with thw passed class variable""" 
        if classVar.varType != orange.VarTypes.Discrete or len(classVar.values) != 2:
            return False
        else:
            return True

    def printParams(self,params):
        print "Boost train parameters:"
        for attr in ("boost_type","weak_count","split_criteria","weight_trim_rate", "max_depth", "use_surrogates","priors"):
            print "%20s" % str(attr)+": ",
            print "%s" % str(getattr(params,attr))   

    def __call__(self, data, weight = None):
        """Creates a Boost model from the data in origTrainingData. """
        if not AZBaseClasses.AZLearner.__call__(self, data, weight):
            return None
        if data.domain.classVar.varType != orange.VarTypes.Discrete:
            print "AZorngCvBoost can only be used for binary classification."
            return None
        #Remove from the domain any unused values of discrete attributes including class
        data = dataUtilities.getDataWithoutUnusedValues(data,True)

        #dataUtilities.rmAllMeta(data) 
        if len(data.domain.getmetas()) == 0:
            trainingData = data
        else:
            trainingData = dataUtilities.getCopyWithoutMeta(data)
        # Create the imputer
        self.imputer = orange.ImputerConstructor_average(trainingData)
        # Impute the data 
        self.trainData = self.imputer(trainingData)

        impData=self.imputer.defaults
        #Convert the ExampleTable to CvMat
        CvMatrices = dataUtilities.ExampleTable2CvMat(self.trainData)
        mat = CvMatrices["matrix"]
        responses = CvMatrices["responses"]
        varTypes = CvMatrices["varTypes"]
        missingDataMask = CvMatrices["missing_data_mask"]

        #Configure Boost params
        #First, Correct any wrong parameters Combination:
        #   CVBOOSTTYPE = { "DISCRETE":0, "REAL":1, "LOGIT":2, "GENTLE":3 }
        #   CVBOOSTSPLITCRIT = { "DEFAULT":0, "GINI":1, "MISCLASS":3, "SQERR":4 }
        if self.boost_type not in AZOC.CVBOOSTTYPE:
            print "ERROR: Bad value for parameter boost_type. Possible values: " + string.join([x for x in AZOC.CVBOOSTTYPE],", ")
            return None
        if self.split_criteria not in AZOC.CVBOOSTSPLITCRIT:
            print "ERROR: Bad value for parameter split_criteria. Possible values: " + string.join([x for x in AZOC.AZOC.CVBOOSTSPLITCRIT],", ")  
            return None

        if self.boost_type == "DISCRETE":
            if self.split_criteria not in ["MISCLASS", "GINI"]:
                print "WARNING: For Discrete type, the split Criteria must be MISCLASS or GINI. MISCLASS was used by default."
                self.split_criteria = "MISCLASS"
        if self.boost_type == "REAL":
            if self.split_criteria not in ["MISCLASS", "GINI"]:
                print "WARNING: For REAL type, the split Criteria must be MISCLASS or GINI. GINI was used by default."
                self.split_criteria = "GINI"
        if self.boost_type in ["LOGIT","GENTLE"]:
            if self.split_criteria != "SQERR":
                print "WARNING: For LOGIT and GENTLE types, the split Criteria must be SQERR. SQERR was used by default."
                self.split_criteria = "SQERR"

        params = ml.CvBoostParams()
        params.boost_type = AZOC.CVBOOSTTYPE[self.boost_type]
        params.split_criteria = AZOC.CVBOOSTSPLITCRIT[self.split_criteria]
        params.weak_count = self.weak_count
        params.weight_trim_rate = self.weight_trim_rate
        params.max_depth = self.max_depth
        params.use_surrogates = self.use_surrogates

        #Create the model it MUST be created with the NON DEFAULT constructor or must call create
        classifier = ml.CvBoost()
        #Train the model
        #train(const CvMat* _train_data, int _tflag, const CvMat* _responses, const CvMat* _var_idx=0, const CvMat* _sample_idx=0, const CvMat* _var_type=0, const CvMat* _missing_mask=0, CvBoostParams params=CvBoostParams(), bool update=false)
        #sampleWeights = cv.cvCreateMat(1,len(self.trainData),cv.CV_32FC1)
        #cv.cvSet(sampleWeights,1.0)
        
        #compute priors (sample weights)
        priors = self.convertPriors(self.priors, self.trainData.domain.classVar)
        if type(priors) == str: #If a string is returned, there was a failure, and it is the respective error mnessage.
            print priors
            return None
        #Train the model
        if self.verbose: self.printParams(params)
        classifier.train(mat, ml.CV_ROW_SAMPLE, responses, None, None, varTypes, missingDataMask, params, False, priors and str(priors).replace(","," ") or None)
        return CvBoostClassifier(classifier = classifier, classVar = self.trainData.domain.classVar, imputeData=impData, verbose = self.verbose, varNames = CvMatrices["varNames"], nIter = None, basicStat = self.basicStat, NTrainEx = len(trainingData), parameters = self.parameters)

class CvBoostClassifier(AZBaseClasses.AZClassifier):
    def __new__(cls, name = "CvBoost classifier", **kwds):
        self = AZBaseClasses.AZClassifier.__new__(cls, name = name,  **kwds)
        #self.__init__(name, **kwds)
        return self

    def getTopImportantVars(self, inEx, nVars = 1, gradRef = None, absGradient = True, c_step = None, getGrad = False):
        return {"NA":"Not aplicable: No true DFV"}


    def __init__(self, name = "CvBoost classifier", **kwds):
        self.verbose = 0
        self.loadedModel = False
        self.__dict__.update(kwds)
        self._isRealProb = False
        self.name = name
        self.domain = None        
        self.ExFix = dataUtilities.ExFix()
        if self.imputeData:
            '''Create the imputer: the imputer needs the imputeData to exists allong it's life time'''
            try:
                self.domain = self.imputeData.domain
                self.imputer = orange.Imputer_defaults(self.imputeData)
            except:
                if self.verbose > 0: print "Unable to create the imputer"
                return None
        else:
            if self.verbose > 0: 
                print "Warning! - No impute data defined"
                return None


    def _singlePredict(self, origExamples = None, resultType = orange.GetValue, returnDFV = False):
        res = None
        """
        orange.GetBoth -          <type 'tuple'>                     ->    (<orange.Value 'Act'='3.44158792'>, <3.442: 1.000>)
        orange.GetValue -         <type 'orange.Value'>              ->    <orange.Value 'Act'='3.44158792'>
        orange.GetProbabilities - <type 'orange.DiscDistribution'>   ->    <0.000, 0.000> 
        """
        #dataUtilities.rmAllMeta(examples)
        if len(origExamples.domain.getmetas()) == 0:
            examples = origExamples
        else:
            examples = dataUtilities.getCopyWithoutMeta(origExamples)
        #Check if the examples are compatible with the classifier (attributes order and varType compatibility)
        dataUtilities.verbose = self.verbose
        if not self.ExFix.ready:
                self.ExFix.set_domain(self.imputer.defaults.domain)
                self.ExFix.set_examplesFixedLog(self.examplesFixedLog) 
        inExamples = self.ExFix.fixExample(examples)

        if not inExamples:
                return None

        #Imput the examples if there are missing values     
        examplesImp = self.imputer(inExamples)
        # There is a problem with using the imputer when examples contain meta attributes.
        # Unable to remove meta attributes from the examples. OK to rm meta from ExampleTables, but not from Example objects.
        if not examplesImp:
            if self.verbose > 0: print "Unable to predict with the Boost model."
            if self.verbose > 0: print "Perhaps you need to remove meta attributes from your examples."
            return None

        out = self.classifier.predict(dataUtilities.Example2CvMat(examplesImp,self.varNames))
        probabilities = None
        DFV = None
        # Back transform the prediction to the original classes and calc probabilities
        prediction = dataUtilities.CvMat2orangeResponse(out, self.classVar)
        # Calculate artificial probabilities - not returned by the OpenCV RF algorithm
        if self.classVar.varType == orange.VarTypes.Discrete:
                #Need to make sure to return meanful probabilities to the cases where opencvRF does not support probabilities
                # to be compatible with possible callers asking for probabilities. 
                probabilities = self.__generateProbabilities(prediction)
                self._isRealProb = False
                probOf1 = probabilities[self.classVar.values[1]]
                DFV = -(probOf1-0.5)
                self._updateDFVExtremes(DFV)

        else:
            #On Regression models assume the DVF as the value predicted
            DFV = prediction
            self._updateDFVExtremes(DFV)

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
  
   
    def __generateProbabilities(self, prediction):
        # Method to artificialy generate a list the length of the number of classes and set the predicted class to 1
        dist = orange.DiscDistribution(self.classVar)
        dist[prediction]=1
        return dist
 

    def write(self, path):
        '''Save a Boost classifier to disk'''
        thePath = str(path)
        try:
            if os.path.isdir(thePath):
                os.system("rm -f "+os.path.join(thePath,"ImputeData.tab"))
                os.system("rm -f "+os.path.join(thePath,"model.boost"))
                os.system("rm -f "+os.path.join(thePath,"varNames.txt"))
            else:
                os.mkdir(thePath)
            if not os.path.isdir(thePath):
                if self.verbose > 0: print "ERROR: Could not create ", path
                return False

            impData = dataUtilities.DataTable(self.imputer.defaults.domain)
            impData.append(self.imputer.defaults)
            # Remove the meta attributes from the imputer data. We don't need to store them along with the model
            impData = dataUtilities.getCopyWithoutMeta(impData)
            impData.save(os.path.join(thePath,"ImputeData.tab"))

            self.classifier.save(os.path.join(thePath,"model.boost"))
            #Save the var names orderes the same way the Learner was trained
            varNamesFile = open(os.path.join(thePath,"varNames.txt"),"w")
            varNamesFile.write(str(self.varNames)+"\n")
            varNamesFile.write(str(self.NTrainEx)+"\n")
            varNamesFile.write(str(self.basicStat)+"\n")
            varNamesFile.close()
            #Save the parameters
            self._saveParameters(os.path.join(thePath,"parameters.pkl"))
        except:
            if self.verbose > 0: print "ERROR: Could not save model to ", path
            return False
        return True
                        
def CvBoostread(path, verbose = 0):
    '''Read a CvBopost classifier from disk and return as a CvBoostClassifier instance. '''
    NTrainEx = 0
    basicStat = None
    thePath = str(path)
    try:
        if not os.path.isdir(thePath):
            if verbose > 0: print "ERROR: no such path:  ", path
            return None
        if not os.path.isfile(str(os.path.join(thePath,"ImputeData.tab"))) or not os.path.isfile(str(os.path.join(thePath,"model.boost"))):
            if verbose > 0: print "ERROR: Missing saved model/data files in:  ", path
            return None
       
        impData = dataUtilities.DataTable(str(os.path.join(thePath,"ImputeData.tab")),createNewOn=orange.Variable.MakeStatus.OK)
        loadedboost = ml.CvBoost()
        loadedboost.load(os.path.join(thePath,"model.boost"))
        
        #Load the var names oredered the way it was used when training
        varNamesFile = open(os.path.join(thePath,"varNames.txt"),"r")
        lines = varNamesFile.readlines()
        varNames = eval(lines[0].strip())
        if len(lines) >= 3:
            NTrainEx = eval(lines[1].strip())
            basicStat = eval(lines[2].strip())
        varNamesFile.close()
        # Read the parameters
        if os.path.isfile(os.path.join(thePath,"parameters.pkl")):
            fileh = open(os.path.join(thePath,"parameters.pkl"),"r")
            parameters = pickle.load(fileh)
            fileh.close()
        else:
            parameters = {} 

        return CvBoostClassifier(classifier = loadedboost, imputeData=impData[0], classVar = impData.domain.classVar, verbose = verbose, loadedModel = True, varNames = varNames, NTrainEx = NTrainEx, basicStat = basicStat, parameters = parameters)
    except:
        if verbose > 0: print "ERROR: Could not read model from ", path







if __name__ == "__main__":
    trainData = dataUtilities.DataTable("../../tests/source/data/BinClass_No_metas_Test.tab")
    learner = CvBoostLearner()
    learner.priors = {"POS":0.7,  "NEG":0.3}
    classifier = learner(trainData)
    preds = {}
    corrects = 0
    for ex in trainData[0:10]:
        pred = classifier(ex)
        print pred," <- ",ex.getclass()
        if pred.value not in preds:
            preds[pred.value] = 1
        else:
            preds[pred.value] +=1
        if pred == ex.getclass():
            corrects +=1
    print "Acc:",corrects/10.0
    for p in preds:
        print p,": ",preds[p]



 
