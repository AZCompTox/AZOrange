import pickle
import orange
import AZBaseClasses
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC
import os
from opencv import ml,cv

class CvBayesLearner(AZBaseClasses.AZLearner):
    """
    Creates an opencv Normal Bayes learner derivated from AZBaseClasses.AZLearner.
    It only support Classification
    All variables are used as ordered 
    It does not support missing measurements (imputation is needed)
    feature vectors from each class should be normally distributed
    """
    def __new__(cls, trainingData = None, name = "CvBayes learner", **kwds):
        self = AZBaseClasses.AZLearner.__new__(cls, **kwds)
        if trainingData:
            self.__init__(name, **kwds)
            return self.__call__(trainingData)
        else:
            self.__dict__.update(kwds)
            self.name = name
            return self

    def __init__(self, name = "CvBayes learner", **kwds):
        self.verbose = 0
        self.name = name
        self.trainData = None
        self.imputer = None
        self.scale = AZOC.CVBAYESDEFAULTDICT["scale"] 
        self.__dict__.update(kwds)


    def __call__(self, data, weight = None):
        """Creates a Bayes model from the data in origTrainingData. """
        if not AZBaseClasses.AZLearner.__call__(self, data, weight):
            return None
        if data.domain.classVar.varType != orange.VarTypes.Discrete:
            raise Exception("AZorngCvBayes can only be used for classification.")
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
        trainingData = self.imputer(trainingData)
        if self.scale:
            self.scalizer = dataUtilities.scalizer()
            self.scalizer.scaleClass = False
            self.scalizer.nMin = -1
            self.scalizer.nMax = 1 
            self.trainData = self.scalizer.scaleAndContinuizeData(trainingData)
        else:
            self.trainData = trainingData
            self.scalizer = None

        impData=self.imputer.defaults
        #Convert the ExampleTable to CvMat
        CvMatrices = dataUtilities.ExampleTable2CvMat(self.trainData)
        mat = CvMatrices["matrix"]
        responses = CvMatrices["responses"]
        varTypes = CvMatrices["varTypes"]
        missingDataMask = CvMatrices["missing_data_mask"]

        #Create the model it MUST be created with the NON DEFAULT constructor or must call create
        classifier = ml.CvNormalBayesClassifier()
        classifier.clear()
        #Train the model
        #CvNormalBayesClassifier::train(const CvMat* _train_data, const CvMat* _responses, const CvMat* _var_idx =0, const CvMat* _sample_idx=0, bool update=false)
        classifier.train(mat, responses, None, None, False)
        return CvBayesClassifier(classifier = classifier, classVar = trainingData.domain.classVar, imputeData=impData, verbose = self.verbose, varNames = CvMatrices["varNames"], nIter = None, basicStat = self.basicStat, NTrainEx = len(trainingData), scalizer = self.scalizer, parameters = self.parameters)

class CvBayesClassifier(AZBaseClasses.AZClassifier):
    def __new__(cls, name = "CvBayes classifier", **kwds):
        self = AZBaseClasses.AZClassifier.__new__(cls, name = name,  **kwds)
        #self.__init__(name, **kwds)
        return self

    def __init__(self, name = "CvBayes classifier", **kwds):
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


    def __call__(self, origExamples = None, resultType = orange.GetValue, returnDFV = False):
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
            if self.verbose > 0: print "Unable to predict with the Bayes model."
            if self.verbose > 0: print "Perhaps you need to remove meta attributes from your examples."
            return None

        if self.scalizer:
            ex = self.scalizer.scaleEx(examplesImp)
        else:
            ex = examplesImp        
        out = self.classifier.predict(dataUtilities.Example2CvMat(ex,self.varNames))
        #print "OUT:",out
        probabilities = None
        DFV = None
        # Back transform the prediction to the original classes and calc probabilities
        prediction = dataUtilities.CvMat2orangeResponse(out, self.classVar)
        #print "Prediction:",prediction
        # Calculate artificial probabilities - not returned by the OpenCV RF algorithm
        if self.classVar.varType == orange.VarTypes.Discrete:
            if resultType != orange.GetValue:
                #Need to make sure to return meanful probabilities to the cases where opencvRF does not support probabilities
                # to be compatible with possible callers asking for probabilities. 
                probabilities = self.__generateProbabilities(prediction)
                self._isRealProb = False
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
        '''Save a Bayes classifier to disk'''
        thePath = str(path)
        try:
            if os.path.isdir(thePath):
                os.system("rm -f "+os.path.join(thePath,"ImputeData.tab"))
                os.system("rm -f "+os.path.join(thePath,"model.bayes"))
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

            self.classifier.save(os.path.join(thePath,"model.bayes"))
            if self.scalizer != None:
                self.scalizer.saveScalingValues(os.path.join(thePath,"scalingValues"))
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
                        
def CvBayesread(path, verbose = 0):
    '''Read a CvBayes classifier from disk and return as a CvBayesClassifier instance. '''
    NTrainEx = 0
    basicStat = None
    thePath = str(path)
    try:
        if not os.path.isdir(thePath):
            if verbose > 0: print "ERROR: no such path:  ", path
            return None
        if not os.path.isfile(str(os.path.join(thePath,"ImputeData.tab"))) or not os.path.isfile(str(os.path.join(thePath,"model.bayes"))):
            if verbose > 0: print "ERROR: Missing saved model/data files in:  ", path
            return None
       
        impData = dataUtilities.DataTable(str(os.path.join(thePath,"ImputeData.tab")),createNewOn=orange.Variable.MakeStatus.OK)
        loadedbayes = ml.CvNormalBayesClassifier()
        loadedbayes.load(os.path.join(thePath,"model.bayes"))

        # Load Scaling Values        
        if os.path.isdir(str(os.path.join(thePath,"scalingValues"))):
            scalizer = dataUtilities.scalizer(file=str(os.path.join(thePath,"scalingValues")))
        else:
            scalizer = None

        
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
        return CvBayesClassifier(classifier = loadedbayes, imputeData=impData[0], classVar = impData.domain.classVar, verbose = verbose, loadedModel = True, varNames = varNames, NTrainEx = NTrainEx, basicStat = basicStat, scalizer = scalizer, parameters = parameters)
    except:
        if verbose > 0: print "ERROR: Could not read model from ", path


