
import orange
import AZBaseClasses
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC
import os
from opencv import ml,cv

class CvSVMLearner(AZBaseClasses.AZLearner):
    """
    Creates a opencv SVM learner derivated from AZBaseClasses.AZLearner. 
    """
    def __new__(cls, trainingData = None, name = "CvSVM learner", **kwds):
        self = AZBaseClasses.AZLearner.__new__(cls, **kwds)
        if trainingData:
            self.__init__(name, **kwds)
            return self.__call__(trainingData)
        else:
            self.__dict__.update(kwds)
            self.name = name
            return self

    def __init__(self, name = "CvSVM learner", **kwds):
        self.verbose = 0
        self.name = name
        self.scalizer = None
        self.trainData = None
        self.nMin = -1
        self.nMax = 1
        self.scaleClass = False
        self.nClassMin = None
        self.nClassMax = None
        self.imputer = None
        self.params = None

	#Read default parameters from AZOrangeConfig.py file
        for par in ("kernel_type", "svm_type","gamma","C","p","epsC","epsR","maxIter",\
	            "stopCrit","nu","scaleData","scaleClass","coef0","degree","priors"):
            setattr(self, par, AZOC.CVSVMDEFAULTDICT[par])
	# if SVMType is set for Classification problems, use epsC else if it is set to regression use epsR
        if self.svm_type in (103,104):
	    self.eps = self.epsR
        else:
            self.eps = self.epsC
        

        self.__dict__.update(kwds)



    def __call__(self, data, weight = None):
        """Creates an SVM model from the data in origTrainingData. """
        if not AZBaseClasses.AZLearner.__call__(self, data, weight):
            if self.verbose > 0: print "Could not create base class instance"
            return None
        dataUtilities.verbose = self.verbose
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
        if self.scaleData:
            self.scalizer = dataUtilities.scalizer()
            for attr in ("nMin","nMax","nClassMin","nClassMax"):
                setattr(self.scalizer, attr, getattr(self, attr))
            #Only scale the class in regression. On classification, set scaleClass to False
            self.scalizer.scaleClass = self.scaleClass  and trainingData.domain.classVar.varType == orange.VarTypes.Continuous or False
            self.scalizer.nClassMin = self.nClassMin
            self.scalizer.nClassMax = self.nClassMax
            self.trainData = self.scalizer.scaleAndContinuizeData(trainingData)
        else:
            self.trainData = trainingData
            self.scalizer = None

        impData=self.imputer.defaults
        #Adjust the svm type according to the problem (regression or classification)
        if self.svm_type != 102:
            if trainingData.domain.classVar.varType == orange.VarTypes.Continuous:
                if self.svm_type in (100,101):
                    self.svm_type += 3
                    self.eps = self.epsR    #Regression eps
            else:
                if self.svm_type in (103,104):
                    self.svm_type -= 3
                    self.eps = self.epsC    #Classification eps
        #Convert the ExampleTable to CvMat
        CvMatices = dataUtilities.ExampleTable2CvMat(self.trainData)
        mat = CvMatices["matrix"]
        responses = CvMatices["responses"]
        varTypes = CvMatices["varTypes"]

        #Configure SVM self.params
        self.params = ml.CvSVMParams()
        self.params.svm_type = self.svm_type
        self.params.kernel_type = self.kernel_type
        self.params.degree = self.degree
        self.params.gamma = self.gamma
        self.params.coef0 = self.coef0
        self.params.C = self.C
        self.params.nu = self.nu
        self.params.p = self.p
        #Process the priors from a str, list or dict to  a valid list 
        priors = self.convertPriors(self.priors,trainingData.domain.classVar)
        if type(priors) == str: #If a string is returned, there was a failure, and it is the respective error mnessage.
            print priors
            return None

        if priors and self.params.svm_type != ml.CvSVM.C_SVC:
            priors = None
            if self.verbose > 0: print "WARNING: The priors will not have any effect. They can only be used with C_SVC SVM-Type."
        elif priors:
            priors = dataUtilities. List2CvMat(priors)

        self.params.class_weights = priors

        term_crit = cv.CvTermCriteria()
        term_crit.type = self.stopCrit #cv.CV_TERMCRIT_EPS  #  or CV_TERMCRIT_ITER
        term_crit.epsilon = self.eps           #Or use:  term_crit.max_iter = x
        term_crit.max_iter = self.maxIter           #Or use:  term_crit.max_iter = x
        self.params.term_crit =  term_crit

        #Create the model
        classifier = ml.CvSVM()
        #Train the model
        #train(trainData, responses, varIdx, SampleIdx, Params)
        classifier.train(mat,responses,None,None,self.params)
        if classifier.get_support_vector_count() < 1:
            print "WARNING: The number of support vectors is 0." 
            print "This could be becasue the margin between the hyper plane and the support vectors has become zero."
            print "Try to modify the parameters controlling the margin. "
            print "For example decrease C or p(regression only)."
            print "No SVM model returned!"
            return None
        else:
            return CvSVMClassifier(classifier = classifier, classVar = data.domain.classVar, scalizer = self.scalizer, imputeData=impData, verbose = self.verbose, varNames = CvMatices["varNames"], basicStat = self.basicStat, NTrainEx = len(trainingData))

    def printParams(self):
        if not self.params:
            print "Paramaters will only be known after training!"
            return
        print "svm_type: ",self.params.svm_type
        print "kernek_type: ", self.params.kernel_type
        print "Degree: ", self.params.degree
        print "Gamma: ", self.params.gamma
        print "coef0: ", self.params.coef0
        print "C: ", self.params.C
        print "Nu: ", self.params.nu
        print "P: ", self.params.p
        print "Priors: ",self.params.class_weights
        print "TermCrit: ",self.params.term_crit.type
        print "    EPS: ",self.params.term_crit.epsilon
        print "    ITER:",self.params.term_crit.max_iter


class CvSVMClassifier(AZBaseClasses.AZClassifier):
    def __new__(cls, name = "CvSVM classifier", **kwds):
        self = AZBaseClasses.AZClassifier.__new__(cls, name = name,  **kwds)
        #self.__init__(name, **kwds)
        return self

    def __init__(self, name = "CvSVM classifier", **kwds):
        self.verbose = 0
        self.loadedModel = False
        self.__dict__.update(kwds)
        self.name = name
        self.domain = None
        self.ExFix = dataUtilities.ExFix()
        if self.imputeData:
            '''Create the imputer: the imputer needs the imputeData to exists allong it's life time'''
            try:
                self.domain = self.imputeData.domain
                self.imputer = orange.Imputer_defaults(self.imputeData)
            except:
                self.imputer = None
                if self.verbose > 0: print "Unable to create the imputer"
        else:
            if self.verbose > 0: print "Warning! - No impute data defined"


    def __call__(self, origExamples = None, resultType = orange.GetValue, returnDFV = False):
        """
        orange.GetBoth -          <type 'tuple'>                     ->    (<orange.Value 'Act'='3.44158792'>, <3.442: 1.000>)
        orange.GetValue -         <type 'orange.Value'>              ->    <orange.Value 'Act'='3.44158792'>
        orange.GetProbabilities - <type 'orange.DiscDistribution'>   ->    <0.000, 0.000>
        returnDFV - Flag indicating to return the Decision Function Value. If set to True, it will encapsulate the original result asked by the keyword resultType and the DFV into a tuple:
                ((<orange.Value 'Act'='3.44158792'>, <3.442: 1.000>), 2.34443)
                (<orange.Value 'Act'='3.44158792'>, 2.34443) 
                (<0.000, 0.000>, 2.34443)
                If it is not a binary classifier, DFV will be equal to None
                DFV will be a value from greater or equal to 0  
        """
        res = None
        #dataUtilities.rmAllMeta(examples)
        if len(origExamples.domain.getmetas()) == 0:
            examples = origExamples
        else:
            examples = dataUtilities.getCopyWithoutMeta(origExamples)

        #Check if the examples are compatible with the classifier (attributes order and varType compatibility)
        if self.imputer:
            dataUtilities.verbose = self.verbose
            if not self.ExFix.ready:
                self.ExFix.set_domain(self.imputer.defaults.domain)
                self.ExFix.set_examplesFixedLog(self.examplesFixedLog)
            inExamples = self.ExFix.fixExample(examples)

            if not inExamples:
                if self.verbose > 0: print "Warning no example. Returning None prediction"
                return None

            #Imput the examples if there are missing values     
            examplesImp = self.imputer(inExamples)
            # There is a problem with using the imputer when examples contain meta attributes.
            # Unable to remove meta attributes from the examples. OK to rm meta from ExampleTables, but not from Example objects.
            if not examplesImp:
                if self.verbose > 0: print "Unable to predict with the SVM model."
                if self.verbose > 0: print "Perhaps you need to remove meta attributes from your examples."
                return None
        else:
            if self.verbose > 0: print "Warning: No Imputer in SVM Classifier"
            examplesImp = examples

        if self.classifier.get_support_vector_count() ==0:
            if self.verbose > 0: print "WARNING:  Support Vectors count is 0 (zero)" 
        DFV = None
        if examplesImp: 
            if self.scalizer:
                exToPredict = dataUtilities.Example2CvMat(self.scalizer.scaleEx(examplesImp,True), self.varNames)
                res = self.classifier.predict(exToPredict)
                res = self.scalizer.convertClass(res)
                if self.classVar.varType != orange.VarTypes.Continuous and len(self.classVar.values) == 2 and returnDFV:
                    DFV = self.classifier.predict(exToPredict, True)
                else:
                    #On Regression models assume the DVF as the value predicted
                    DFV = res 
                self._updateDFVExtremes(DFV)
                res = dataUtilities.CvMat2orangeResponse(res,self.classVar)
            else:
                exToPredict = dataUtilities.Example2CvMat(examplesImp,self.varNames)
                res = self.classifier.predict(exToPredict)
                if self.classVar.varType != orange.VarTypes.Continuous and len(self.classVar.values) == 2 and returnDFV:
                    DFV = self.classifier.predict(exToPredict, True)
                else:
                    #On Regression models assume the DVF as the value predicted
                    DFV = res
                self._updateDFVExtremes(DFV)
                res = dataUtilities.CvMat2orangeResponse(res,self.classVar)
             
            if resultType!=orange.GetValue:
                if examplesImp.domain.classVar.varType != orange.VarTypes.Continuous:
                    dist = orange.DiscDistribution(examplesImp.domain.classVar)
                    dist[res]=1
                else:
                    dist = res
                if resultType==orange.GetProbabilities:
                    res = dist
                else:
                    res = (res,dist)
                    
            if returnDFV:
                res = (res,DFV)
                
        self.nPredictions += 1
        return res

    def write(self, path):
        '''Save an SVM classifier to disk'''
        thePath = str(path)
        try:
            if os.path.isdir(thePath):
                os.system("rm -f "+os.path.join(thePath,"ImputeData.tab"))
                os.system("rm -f "+os.path.join(thePath,"model.svm"))
                os.system("rm -Rf "+os.path.join(thePath,"scalingValues"))
                os.system("rm -f "+os.path.join(thePath,"varNames.txt"))
                #if os.path.isdir(thePath):
                    #print "ERROR: Cannot overwrite ", path
                    #return False
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

            self.classifier.save(os.path.join(thePath,"model.svm"))
            if self.scalizer != None:
                self.scalizer.saveScalingValues(os.path.join(thePath,"scalingValues")) 
            #Save the var names orderes the same way the Learner was trained
            varNamesFile = open(os.path.join(thePath,"varNames.txt"),"w")
            varNamesFile.write(str(self.varNames)+"\n")
            varNamesFile.write(str(self.NTrainEx)+"\n")
            varNamesFile.write(str(self.basicStat)+"\n")
            varNamesFile.close()

        except:
            if self.verbose > 0: print "ERROR: Could not save model to ", path
            return False
        return True
                        
def CvSVMread(path, verbose = 0):
    '''Read a PLS classifier from disk and return as a PLSClassifier instance. '''
    thePath = str(path)
    NTrainEx = 0
    basicStat = None 
    try:
        if not os.path.isdir(thePath):
            if verbose > 0: print "ERROR: no such path:  ", path
            return None
        if not os.path.isfile(str(os.path.join(thePath,"ImputeData.tab"))) or not os.path.isfile(str(os.path.join(thePath,"model.svm"))):
            if verbose > 0: print "ERROR: Missing saved model/data files in:  ", path
            return None
        if os.path.isdir(str(os.path.join(thePath,"scalingValues"))):
            scalizer = dataUtilities.scalizer(file=str(os.path.join(thePath,"scalingValues"))) 
        else:
            scalizer = None
       
        impData = dataUtilities.DataTable(str(os.path.join(thePath,"ImputeData.tab")),createNewOn=orange.Variable.MakeStatus.OK)
        loadedsvm = ml.CvSVM()
        loadedsvm.load(os.path.join(thePath,"model.svm"))
        
        #Load the var names oredered the way it was used when training
        if (os.path.isfile(os.path.join(thePath,"varNames.txt"))):
            varNamesFile = open(os.path.join(thePath,"varNames.txt"),"r")
            lines = varNamesFile.readlines()
            varNames = eval(lines[0].strip())
            if len(lines) >= 3:
                NTrainEx = eval(lines[1].strip())
                basicStat = eval(lines[2].strip())
            varNamesFile.close()
            #thisVer = True
        else:
            if verbose > 0: print "WARNING: The model loaded was probably saved with azorange version 0.2.1 or lower"
            varNames = [attr.name for attr in impData.domain.attributes]
            #thisVer = False

        return CvSVMClassifier(classifier = loadedsvm, scalizer = scalizer,imputeData=impData[0], classVar = impData.domain.classVar, verbose = verbose, loadedModel = True, varNames = varNames, NTrainEx = NTrainEx, basicStat = basicStat)
    except:
        if verbose > 0: print "ERROR: Could not read model from ", path


