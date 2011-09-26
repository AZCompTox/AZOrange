import pickle
import orange
import AZBaseClasses
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC
import os
from opencv import ml,cv
from AZutilities import evalUtilities

class CvANNLearner(AZBaseClasses.AZLearner):
    """
    Creates a opencv ANN learner derivated from AZBaseClasses.AZLearner. 
    """
    def __new__(cls, trainingData = None, name = "CvANN learner", **kwds):
        self = AZBaseClasses.AZLearner.__new__(cls, **kwds)
        if trainingData:
            self.__init__(name, **kwds)
            return self.__call__(trainingData)
        else:
            self.__dict__.update(kwds)
            self.name = name
            return self

    def __init__(self, name = "CvANN learner", **kwds):
        self.verbose = 0
        self.name = name
        self.trainData = None
        self.imputer = None
	#Read default parameters from AZOrangeConfig.py file
        for par in ("activationFunction","sigmoidAlpha","sigmoidBeta","nHidden", "scaleData",
                    "scaleClass", "optAlg", "bp_dw_scale", 
                    "bp_moment_scale", "rp_dw0", "rp_dw_plus",
                    "rp_dw_minus", "rp_dw_max",'stopCrit', 'maxIter','eps','priors',
                    "nDiffIniWeights","stopUPs"): 
            setattr(self, par, AZOC.CVANNDEFAULTDICT[par])
        self.__dict__.update(kwds)

    def __call__(self, data, weight = None):
        bestSeed = None
        bestAcc = None
        bestNiter = None
        bestModel = None
        #fix self.nDiffIniWeights for the disabled mode
        if self.nDiffIniWeights <= 1:
            self.nDiffIniWeights = 1 #loop over n different initial weights Disabled
        #Fix self.stopUPs for the disabled mode
        if self.stopUPs <=0:
            self.stopUPs = 0  # Optimization of nIter will be disabled

        #Remove from the domain any unused values of discrete attributes including class
        data = dataUtilities.getDataWithoutUnusedValues(data,True)
        #dataUtilities.rmAllMeta(data) 
        if len(data.domain.getmetas()) == 0:
            cleanedData = data
        else:
            cleanedData = dataUtilities.getCopyWithoutMeta(data)
        # Create the imputer
        self.imputer = orange.ImputerConstructor_average(cleanedData)
        # Impute the data 
        self.trainData = self.imputer(cleanedData)
         # If we are not seetin neither weights init optimization or nEphocs optimization (opencvLayer), the do nto split the data
        if self.stopUPs != 0 or self.nDiffIniWeights > 1:
            #Define train-80% and validation set-20% of the input data
            indices = orange.MakeRandomIndices2(p0=0.2, stratified = orange.MakeRandomIndices.StratifiedIfPossible)
            ind = indices(cleanedData)
            self.trainData = cleanedData.select(ind,1)
            validationSet = cleanedData.select(ind,0)
        else:
            validationSet = None

        if self.verbose and self.nDiffIniWeights>1: print "=========== Training ",self.nDiffIniWeights," times with different initial weights =============="
        for n in range(self.nDiffIniWeights):
            if self.nDiffIniWeights <=1:
                seed=0  #in opencv  mmlann seed=0 means the seed is disabled, and original seed will be used
            else:
                seed = len(cleanedData) * len(cleanedData.domain) * (n+1)  #seed can be any integer
            #Create a model with a specific seed for training opencv ANN. 
            #Also passing the step for the nIter optimization (self.stopUPs=0 - disable nIter optimization)
            #Also passing the validation set to be used in internal opencv implemented nEphocs optimization.
            model = self.__train__(weight = None, seed = seed, validationSet = validationSet)
            #Skip evaluation if the weights loop is disabled
            if self.nDiffIniWeights <=1:
                return model
                break
            if cleanedData.domain.classVar.varType == orange.VarTypes.Discrete:
                Acc = evalUtilities.getClassificationAccuracy(validationSet, model)
            else:
                Acc = -evalUtilities.getRMSE(validationSet, model)
            if bestModel == None or (Acc > bestAcc) or (Acc == bestAcc and model.nIter < bestNiter):
                bestSeed = seed
                bestAcc = Acc
                bestNiter = model.nIter
                bestModel = model
            if self.verbose:  print "nIter:%-7s  Acc:%-20s  seed: %s" % (model.nIter,Acc,seed)

        if self.verbose: print "================ Best model Found: ==================="
        if self.verbose: print "nIter:%-7s  Acc:%-20s  seed: %s" % (bestNiter,bestAcc,bestSeed)

        # DEBUG for check if the returned model is indeed the best model, and not the last trainted
        #if cleanedData.domain.classVar.varType == orange.VarTypes.Discrete:
        #    Acc = evalUtilities.getClassificationAccuracy(validationSet, bestModel)
        #else:
        #    Acc = -evalUtilities.getRMSE(validationSet, bestModel)
        #if self.verbose: print "================ Best model returned: ==================="
        #if self.verbose:  print "nIter:%-7s  Acc:%-20s  seed: %s" % (bestModel.nIter,Acc,bestModel.seed)

        return bestModel

    def __train__(self, weight = None, seed=0, validationSet=None):
        if not AZBaseClasses.AZLearner.__call__(self, self.trainData, weight):
            return None
        """Creates an ANN model from the data in origTrainingData. """


        #Convert the ExampleTable to CvMat
        CvMatices = dataUtilities.ExampleTable2CvMat(self.trainData, True)
        mat = CvMatices["matrix"]
        responses = CvMatices["responses"]
        varTypes = CvMatices["varTypes"]

        #Configure ANN params
        params = ml.CvANN_MLP_TrainParams()
        params.train_method = self.optAlg 
        params.bp_dw_scale = self.bp_dw_scale
        params.bp_moment_scale = self.bp_moment_scale
        params.rp_dw0 = self.rp_dw0
        params.rp_dw_plus = self.rp_dw_plus
        params.rp_dw_minus = self.rp_dw_minus
        #params.rp_dw_min = ##default is the minimum float value 
        params.rp_dw_max = self.rp_dw_max

        term_crit = cv.CvTermCriteria()
        term_crit.type = self.stopCrit #cv.CV_TERMCRIT_EPS  #  or CV_TERMCRIT_ITER
        term_crit.epsilon = self.eps           #Or use:  term_crit.max_iter = x
        term_crit.max_iter = self.maxIter           #Or use:  term_crit.max_iter = x
        params.term_crit =  term_crit

        #Create the model it MUST be created with the NON DEFAULT constructor or must call create
        classifier = ml.CvANN_MLP()
        if self.trainData.domain.classVar.varType == orange.VarTypes.Discrete:
            Nout = len(self.trainData.domain.classVar.values)
        else:
            Nout = 1
        if type(self.nHidden) != list: 
            nHidden = [self.nHidden]
        else:
            nHidden = self.nHidden
        layers = [len(self.trainData.domain.attributes)] + nHidden + [Nout]
        layerSizes = dataUtilities.List2CvMat(layers,"CV_32SC1")
        classifier.create(layerSizes, self.activationFunction, self.sigmoidAlpha, self.sigmoidBeta)
        #Train the model
        #train(trainData, responses, sampleWeights (RPROP only), sampleIdx, TrainParams, flags for scaling)
        #sampleWeights = cv.cvCreateMat(1,len(self.trainData),cv.CV_32FC1)
        #cv.cvSet(sampleWeights,1.0)
        
        scaleFlag = 0
        if not self.scaleData:
            scaleFlag = scaleFlag | ml.CvANN_MLP.NO_INPUT_SCALE
        if not self.scaleClass:
            scaleFlag = scaleFlag |  ml.CvANN_MLP.NO_OUTPUT_SCALE
       
        #compute priors (sample weights)
        priors = self.convertPriors(self.priors, self.trainData.domain.classVar,getDict = True)
        if type(priors) == str: #If a string is returned, there was a failure, and it is the respective error mnessage.
            print priors
            return None
 
        if priors and self.optAlg == 1:
            #scale priors
            pSum=sum(priors.values())
            if pSum==0:
                print "ERROR: The priors cannot be all 0!"
                return None
            map(lambda k,v:priors.update({k: (v+0.0)/pSum}),priors.keys(),priors.values())
            #Apply the priors to each respective sample
            sample_weights = [1] * len(self.trainData)
            for idx,sw in enumerate(sample_weights):
                actualClass = str(self.trainData[idx].getclass().value)
                if actualClass in priors:
                    sample_weights[idx] = sample_weights[idx] * priors[actualClass]
            CV_sample_weights = dataUtilities.List2CvMat(sample_weights,"CV_32FC1")
        else:
            CV_sample_weights = None
        #Train the model
        if self.stopUPs > 0:
            #Convert the ExampleTable to CvMat
            CvMatices = dataUtilities.ExampleTable2CvMat(validationSet, True)
            valMat = CvMatices["matrix"]
            valResponses = CvMatices["responses"]
        else:
            valMat = None
            valResponses = None
        nIter = classifier.train(mat, responses, CV_sample_weights, None, params, scaleFlag, seed,self.stopUPs,valMat,valResponses)
        model = CvANNClassifier(seed = seed, classifier = classifier, classVar = self.trainData.domain.classVar,
                        imputeData=self.imputer.defaults, verbose = self.verbose, varNames = CvMatices["varNames"],
                        nIter = nIter, basicStat = self.basicStat, NTrainEx = len(self.trainData), parameters = self.parameters)
        return model


class CvANNClassifier(AZBaseClasses.AZClassifier):
    def __new__(cls, name = "CvANN classifier", **kwds):
        self = AZBaseClasses.AZClassifier.__new__(cls, name = name,  **kwds)
        #self.__init__(name, **kwds)
        return self

    def __init__(self, name = "CvANN classifier", **kwds):
        self.verbose = 0
        self.loadedModel = False
        self.__dict__.update(kwds)
        self._isRealProb = True     #This learner always return real probabilities
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
        if self.verbose > 1: dataUtilities.verbose = self.verbose
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
            if self.verbose > 0: print "Unable to predict with the ANN model."
            if self.verbose > 0: print "Perhaps you need to remove meta attributes from your examples."
            return None

        res = None
        if self.classVar.varType == orange.VarTypes.Continuous: 
                Nout = 1
        else:
                Nout = len(self.classVar.values)
        out = cv.cvCreateMat(1,Nout,cv.CV_32FC1)
        self.classifier.predict(dataUtilities.Example2CvMat(examplesImp,self.varNames),out)
        #print "OUT = ",out
        #print out,"->",dataUtilities.CvMat2orangeResponse(out,self.classVar,True),":",origExamples[self.classVar.name].value
        res = dataUtilities.CvMat2orangeResponse(out,self.classVar,True)
        #print "RES=",res
        DFV = None
        if out.cols > 1:
            fannOutVector = dataUtilities.CvMat2List(out)[0]
            probabilities = self.__getProbabilities(fannOutVector)
            #Compute the DFV
            if self.classVar.varType == orange.VarTypes.Discrete and len(self.classVar.values) == 2:
                DFV = probabilities[0]
                # Subtract 0.5 so that the threshold is 0 as all learners DFV
                DFV -= 0.5
                self._updateDFVExtremes(DFV)    
        
            # Retrun the desired quantity
            if resultType == orange.GetProbabilities:
                    res = probabilities
            else:
                    if resultType == orange.GetBoth:
                        res = (res, probabilities)
        else:
            #On Regression models, assume the DFV as the value predicted
            DFV = res.value
            self._updateDFVExtremes(DFV)
            if resultType == orange.GetProbabilities:
                res  =  [0.0]
            else:
                if resultType==orange.GetBoth:
                    res = (res,[0.0])
        self.nPredictions += 1
        
        if returnDFV:
            return (res,DFV)
        else:
            return res
    
    def __getProbabilities(self, fannOutVector):
        """Get the orange like output probabilities for the current predicted example"""
        dist = orange.DiscDistribution(self.classVar)
        vectorSum = sum(fannOutVector)
        #fix the probabilities so that values are between 0 and 1
        OutVector = [p/vectorSum for p in fannOutVector]
        subtract = abs(sum([x for x in OutVector if x < 0]))
        for idx,p in enumerate(OutVector):
            if p > 1:
                dist[self.classVar.values[idx]] = p - subtract
            elif p <= 0:
                dist[self.classVar.values[idx]] = 0
            else:
                dist[self.classVar.values[idx]] = p
        return dist



    def write(self, path):
        '''Save an ANN classifier to disk'''
        thePath = str(path)
        try:
            if os.path.isdir(thePath):
                os.system("rm -f "+os.path.join(thePath,"ImputeData.tab"))
                os.system("rm -f "+os.path.join(thePath,"model.ann"))
                os.system("rm -f "+os.path.join(thePath,"varNames.txt"))
                #os.system("rm -R "+thePath)
                #if os.path.isdir(thePath):
                    #if self.verbose > 0: print "ERROR: Could not overwrite ", path
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

            self.classifier.save(os.path.join(thePath,"model.ann"))
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
                        
def CvANNread(path, verbose = 0):
    '''Read a PLS classifier from disk and return as a PLSClassifier instance. '''
    NTrainEx = 0
    basicStat = None
    thePath = str(path)
    try:
        if not os.path.isdir(thePath):
            if verbose > 0: print "ERROR: no such path:  ", path
            return None
        if not os.path.isfile(str(os.path.join(thePath,"ImputeData.tab"))) or not os.path.isfile(str(os.path.join(thePath,"model.ann"))):
            if verbose > 0: print "ERROR: Missing saved model/data files in:  ", path
            return None
       
        impData = dataUtilities.DataTable(str(os.path.join(thePath,"ImputeData.tab")),createNewOn=orange.Variable.MakeStatus.OK)
        loadedann = ml.CvANN_MLP()
        loadedann.load(os.path.join(thePath,"model.ann"))
        
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
        # Read the parameters
        if os.path.isfile(os.path.join(thePath,"parameters.pkl")):
            fileh = open(os.path.join(thePath,"parameters.pkl"),"r")
            parameters = pickle.load(fileh)
            fileh.close()
        else:
            parameters = {} 

        return CvANNClassifier(classifier = loadedann, imputeData=impData[0], classVar = impData.domain.classVar, verbose = verbose, loadedModel = True, varNames = varNames, NTrainEx = NTrainEx, basicStat = basicStat, parameters = parameters)
    except:
        if verbose > 0: print "ERROR: Could not read model from ", path


