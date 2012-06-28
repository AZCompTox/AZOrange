"""
AZBaseClasses
Base Classes for AZ methods.
"""
import imp
import orange
import types,os
from AZutilities import dataUtilities
import pickle
#from trainingMethods import AZorngConsensus 

def getCorrespondingLearner(modelPath, getParameters = True):
    """ Determines what is the learner used to build the model in modelPath.
        if the Learner is to be optimized, one can set getParameters = False so that the model does not have to be loaded (if it would be needed when there is no parameters.pkl for loading)
        It allows the caller to do not have to load the necessary training method(s)
        It loads the necessary learner module(s) and returns the Learner with parameters set like the corresponding one.
        
        Internal getCorrespondingLearner variables:
                if single model in modelPath:
                                ex.:  "learners":  {"RF": <RFLearner 'RF learner'>}   
                                      "parameters" : {"RF":{"nActVars":5, "maxDepth":20}}
                                    
                if consensus model in modelPath and expression is present:
                               ex.:  "learners":  {"RF": <RFLearner 'RF learner'> , 
                                                    "SVM":<RFLearner 'CvSVM learner'>, 
                                                    "Consensus":<ConsensusLearner 'Consensus learner'>, 
                                                    ...} 
                                      "Cexpression" : ['RF == POS and SVM == POS->POS', '->NEG']
                                      "parameters" : {"RF":{"nActVars":5, "maxDepth":20},
                                                      "CvSVM": {"C":2, "svm_type":"RBF"},
                                                      ...}  
                                    
                if consensus model in modelPath and NO expression is present  (for Backcompatibility):
                               ex.:  "learners":  {"RF": <RFLearner 'RF learner'> , 
                                                    "CvSVM":<RFLearner 'CvSVM learner'>, 
                                                    "Consensus":<ConsensusLearner 'Consensus learner'>, 
                                                    ...} 
                                      "Cexpression" : []
                                      "parameters" : {"RF":{"nActVars":5, "maxDepth":20},
                                                      "CvSVM": {"C":2, "svm_type":"RBF"},
                                                      ...}
                                     

    """
    learners = {}
    Cexpression = []
    parameters = {}

    learnersDict = {}
    theLearner = None

    modelType = modelRead(modelFile=modelPath, retrunClassifier = False)
    exec "from trainingMethods import AZorng" + modelType
    if modelType != "Consensus":
        learners[modelType] = eval("AZorng"+modelType+"."+modelType+"Learner()")        
        if getParameters:
            if os.path.isfile(os.path.join(modelPath,"parameters.pkl")):
                fileh = open(os.path.join(modelPath,"parameters.pkl"))
                parameters[modelType] = pickle.load(fileh)
                fileh.close()
            else:
                model = modelRead(modelPath)
                parameters[modelType] = model.parameters
        else:
            parameters[modelType] = {}
        for par in parameters[modelType]:
            learners[modelType].setattr(par, parameters[modelType][par])
        theLearner = learners[modelType]
    else:
        if os.path.isfile(os.path.join(modelPath,"learnerDict.pkl")) and (os.path.isfile(os.path.join(modelPath,"expression.pkl")) or os.path.isfile(os.path.join(modelPath,"expressionList.pkl"))):
            if os.path.isfile(os.path.join(modelPath,"expression.pkl")):
                fh = open(os.path.join(modelPath,"expression.pkl"))
            else:
                fh = open(os.path.join(modelPath,"expressionList.pkl"))
            Cexpression = pickle.load(fh)
            fh.close()

            fh = open(os.path.join(modelPath,"learnerDict.pkl"))
            learnersDict = pickle.load(fh)
            fh.close()
        if not Cexpression or not learnersDict:  # Trivial Consensus Model
            modelPaths = [dir for dir in os.listdir(modelPath) if "C" in dir and ".model" in dir]
            for modelP in modelPaths:
                mType = modelRead(modelFile=os.path.join(modelPath,modelP), retrunClassifier = False)
                exec "from trainingMethods import AZorng" + mType
                learners[mType] = eval("AZorng"+mType+"."+mType+"Learner()") 
                if getParameters:
                    if os.path.isfile(os.path.join(modelPath,modelP,"parameters.pkl")):
                        fileh = open(os.path.join(modelPath,modelP,"parameters.pkl"))
                        parameters[mType] = pickle.load(fileh)
                        fileh.close()
                    else:
                        model = modelRead(os.path.join(modelPath,modelP))
                        parameters[mType] = model.parameters
                else:
                    parameters[mType] = {}
            for l in learners:
                for par in parameters[l]:
                    learners[l].setattr(par, parameters[l][par])
            exec "from trainingMethods import AZorngConsensus"
            theLearner = AZorngConsensus.ConsensusLearner(learners = [learners[l] for l in learners])
        else:                                   # advanced Consensus model
            for modelName in learnersDict:
                modelP = "C"+str(learnersDict[modelName])+".model"
                mType = modelRead(modelFile=os.path.join(modelPath,modelP), retrunClassifier = False)
                exec "from trainingMethods import AZorng" + mType
                learners[modelName] = eval("AZorng"+mType+"."+mType+"Learner()")            
                if getParameters:
                    if os.path.isfile(os.path.join(modelPath,modelP,"parameters.pkl")):
                        fileh = open(os.path.join(modelPath,modelP,"parameters.pkl"))
                        parameters[mType] = pickle.load(fileh)
                        fileh.close()
                    else:
                        model = modelRead(os.path.join(modelPath,modelP))
                        parameters[modelName] = model.parameters
                else:
                    parameters[modelName] = {}
            for l in learners:
                for par in parameters[l]:
                    learners[l].setattr(par, parameters[l][par])
            exec "from trainingMethods import AZorngConsensus"
            theLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = Cexpression)


    return theLearner

def getModelDomain(modelPath):
    """
    Looks for the domain used to train the model. 
    If looks for the file ImputeData.tab where it extracts the model domain
    If it does not exist, it loads the model and get the parameters from the variable domain
    if can't do it, returns None
    """
    try:
        filePath = os.path.join(modelPath, "ImputeData.tab")
        if os.path.isfile(filePath):
            impData = dataUtilities.DataTable(filePath)
            return impData.domain
        else:
            model = modelRead(modelPath)
            return model.domain
    except:
        print "ERROR: Can't find the domain of the model in ", modelPath
        return None

 


def modelRead(modelFile=None,verbose = 0,retrunClassifier = True):
    """Get the type of model saved in 'modelPath' and loads the respective model
       Returns the Classifier saved in the respective model path
       If called without parameters, it returns a list of known classifier types
       It can returns the classifier, or just a string with the Type

            modelRead (modelFile [, verbose = 0] [, retrunClassifier = True] )"""

    if not modelFile:
        return ("CvSVM", "CvANN", "PLS", "CvRF", "CvBoost", "CvBayes", "Consensus")

    modelType = None
    loadedModel = None
    if os.path.isfile(os.path.join(modelFile,"model.svm")):
        modelType =  "CvSVM"
        if not retrunClassifier: return modelType
        from trainingMethods import AZorngCvSVM    
        loadedModel = AZorngCvSVM.CvSVMread(modelFile,verbose)
    elif os.path.isfile(os.path.join(modelFile,"model.ann")):
        modelType =  "CvANN"
        if not retrunClassifier: return modelType
        from trainingMethods import AZorngCvANN
        loadedModel = AZorngCvANN.CvANNread(modelFile,verbose)
    elif os.path.isfile(os.path.join(modelFile,"Model.pls")):
        modelType =  "PLS"
        if not retrunClassifier: return modelType
        from trainingMethods import AZorngPLS
        loadedModel = AZorngPLS.PLSread(modelFile,verbose)
    elif os.path.isfile(os.path.join(modelFile,"model.rf")):
        modelType =  "RF"
        if not retrunClassifier: return modelType
        from trainingMethods import AZorngRF
        loadedModel = AZorngRF.RFread(modelFile,verbose)
    elif os.path.isdir(os.path.join(modelFile,"C0.model")):
        modelType =  "Consensus"
        if not retrunClassifier: return modelType
        from trainingMethods import AZorngConsensus
        loadedModel = AZorngConsensus.Consensusread(modelFile,verbose)
    elif os.path.isfile(os.path.join(modelFile,"model.boost")):
        modelType =  "CvBoost"
        if not retrunClassifier: return modelType
        from trainingMethods import AZorngCvBoost
        loadedModel = AZorngCvBoost.CvBoostread(modelFile,verbose)
    elif os.path.isfile(os.path.join(modelFile,"model.bayes")):
        modelType =  "CvBayes"
        if not retrunClassifier: return modelType
        from trainingMethods import AZorngCvBayes
        loadedModel = AZorngCvBayes.CvBayesread(modelFile,verbose)
    else:   # Assuming an RF old format for backcompatibility
        try:
            if os.path.isdir(modelFile):
                modelType =  "RF"
                if not retrunClassifier: return modelType
                from trainingMethods import AZorngRF
                loadedModel = AZorngRF.RFread(modelFile,verbose)
            else:
                modelType = None
                loadedModel = None
        except:
            modelType = None
            loadedModel = None

    return loadedModel
 


class AZLearner(orange.Learner):
    """
    Base Class for AZLearners
    It assures that all AZLearners have al teast the setattr method
    Here we can add more methods commun to all Learners which derivates from this class
    """
    def __new__(cls, trainingData = None, name = "AZ learner", **kwds):
        self = orange.Learner.__new__(cls, **kwds)
        self.__dict__.update(kwds)
        self.name = name
        self.basicStat = None
        self.parameters = {}
        return self
      
    def isCompatible(self, classVar):
        """Checks if the learner is compatiblue with thw passed class variable"""
        return True

 
    def __call__(self, trainingData = None, weight = None): 
        self.basicStat = None
        if not trainingData:
            print "AZBaseClasses ERROR: Missing training data!"
            return False
        elif dataUtilities.findDuplicatedNames(trainingData.domain):
            print "AZBaseClasses ERROR: Duplicated names found in the training data. Please use the method dataUtilities.DataTable() when loading a dataset in order to fix the duplicated names and avoid this error."
            return False
        elif not trainingData.domain.classVar:
            print "AZBaseClasses ERROR: No class attribute found in training data!"
            return False
        elif not len(trainingData):
            print "AZBaseClasses ERROR: No examples in training data!"
            return False
        elif not len(trainingData.domain.attributes):
            print "AZBaseClasses ERROR: No attributes in training data!"
            return False



        possibleMetas = dataUtilities.getPossibleMetas(trainingData, checkIndividuality = True)
        if possibleMetas:
            msg="\nAZBaseClasses ERROR: Detected attributes that should be considered meta-attributes:"
            for attr in possibleMetas:
                msg += "\n    "+attr
            raise Exception(msg)
            #return False
        #Get the Domain basic statistics and save only the desired info in self.basicStat
        basicStat = orange.DomainBasicAttrStat(trainingData)
        self.basicStat = {}
        for attr in trainingData.domain:
            if attr.varType == orange.VarTypes.Discrete:
                self.basicStat[attr.name] = None
            else:       
                self.basicStat[attr.name] = {"min":basicStat[attr].min, "max":basicStat[attr].max, "avg":basicStat[attr].avg}
        # Gather all the learner parameters to be stored along with the classifier 
        # Find the name of the Learner
        learnerName = str(self.__class__)[:str(self.__class__).rfind("'")].split(".")[-1] 
        self.parameters = {}
        if learnerName != "ConsensusLearner":
            # Load the AZLearnersParamsConfig.py from the AZORANGEHOME!
            AZOLearnersConfig = imp.load_source("AZLearnersParamsConfig", os.path.join(os.environ["AZORANGEHOME"],'azorange',"AZLearnersParamsConfig.py"))
            pars = AZOLearnersConfig.API(learnerName)
            if pars:
                for par in pars.getParameterNames():
                    self.parameters[par] = getattr(self,par)
        return True

    def setattr(self, name, value):
        self.__dict__[name] = value


    def convertPriors(self,priors,classVar,getDict = False):
        """Converts the passed priors to a list according to the classVar
              a) returns a list if success convertion
              b) returns None if no priors are to be used, or if the class is not discrete
              b) returns a string with the error message if failed
        """
        if not priors or not classVar or classVar.varType != orange.VarTypes.Discrete:
            return None
        else:
            if type(priors) == str:
                try:
                    InPriors = eval(priors)
                    if InPriors != None and (type(InPriors) not in (dict,list)):
                        raise Exception("ERROR: Priors were specifyed incorrectly! Use a string defining None, dict or a list in python syntax")
                except:
                    InPriors = self.__convertStrPriors2Dict(priors)
                    if not InPriors:
                        raise Exception("ERROR: Priors were specifyed incorrectly. Use a formated sting. Ex: 'POS:2, NEG:1'")
            else:
                InPriors = priors
            #If priors are defined, they are now a List or a Dict
            if type(InPriors) == dict:
                if len(InPriors) != len(classVar.values) or [x for x in [str(v) for v in classVar.values] if x in InPriors] != [str(v) for v in classVar.values]:
                    raise Exception("Error: Wrong priors: "+str(InPriors.keys())+"\n"+\
                           "       Acepted priors: "+str(classVar.values))
                # Create an empty list of priors
                priorsList = [None]*len(classVar.values)
                for cls_idx,cls_v in enumerate(classVar.values):
                        priorsList[cls_idx] = InPriors[cls_v]
                if getDict:
                    return InPriors
                else:
                    return priorsList
            elif type(InPriors) == list:
                if len(InPriors) != len(classVar.values):
                    raise Exception("ERROR: The number of priors specified are not according to the number of class values.\n"+\
                           "       Available Class values :" + str(classVar.values))
                else:
                    if getDict:
                        priorsDict = {}
                        map(lambda k,v: priorsDict.update({k: v}),[str(x) for x in classVar.values],InPriors)
                        return priorsDict
                    else:
                        return InPriors
            else:
                return None

    def __convertStrPriors2Dict(self,priorStr):
        """Convert a string with priors ex: "POS:2, NEG:1" to the correspondent dictionary.
           Returnes the dictionary or None if there was an error"""
        #Validate and convert the priors
        try:
            priorsDict = {}
            priorsList = priorStr.split(",")
            for prior in priorsList:
                priorEl = prior.strip().split(":")
                if len(priorEl) != 2:
                    return None
                else:
                    priorsDict[priorEl[0]] = float(priorEl[1])
            if priorsDict:
                return priorsDict
            else:
                return None
        except:
            return None


class AZClassifier(object):
    """
    Base Class for all AZClassifiers
    Here we can add methods that are commun to all Classifiers that derivates from this class
    """
    def __new__(cls,name = "AZ classifier", **kwds):
        self = object.__new__(cls)
        self.parameters = {}
        self.__dict__.update(kwds)
        self.examplesFixedLog = {}
        self.name = name
        self.classifier = None
        self.domain = None
        self._isRealProb = False
        self._DFVExtremes = {"min":0, "max":0}  # Store for the Extremes of DFV
        self.nPredictions = 0                    # Number of predictions made with this Classifier
        self.basicStat = None
        self.NTrainEx = 0
        return self


    def _saveParameters(self, path):
        fileh = open(path, 'w') 
        pickle.dump(self.parameters,fileh)
        fileh.close()
 
    def _updateDFVExtremes(self, DFV):
        if DFV < self._DFVExtremes["min"]:
            self._DFVExtremes["min"] = DFV
        if DFV > self._DFVExtremes["max"]:
            self._DFVExtremes["max"] = DFV

    def isRealProb(self):
        return self._isRealProb

    def getDFVExtremes(self):
        #If the extremes were never set, return None
        if self._DFVExtremes["min"]==0 and self._DFVExtremes["max"] == 0:
            return None
        else:
            return self._DFVExtremes

    def resetCounters(self):
        self._DFVExtremes = {"min":0, "max":0}
        self.nPredictions = 0

    def getNTrainEx(self):
        """Return the number of examples used at Train time"""
        return self.NTrainEx


    def getTopImportantVars(self, inEx, nVars = 1, gradRef = None, absGradient = True):
        """Return the n top important variables (n = nVars) for the given example
            if nVars is 0, it returns all variables ordered by importance
        """
        varGradValUP = []
        varGradValDOWN = []
        varGradNameUP = []
        varGradNameDOWN = []


        ExFix = dataUtilities.ExFix()
        ExFix.set_domain(self.domain)
        ex = ExFix.fixExample(inEx)
        if self.basicStat == None or not self.NTrainEx or (self.domain.classVar.varType == orange.VarTypes.Discrete and len(self.domain.classVar.values)!=2):
            return None
        if gradRef == None:
            gradRef = self(ex,returnDFV = True)[1]
        # the calcDiscVarGrad and calcContVarGrad will return 
        def calcDiscVarGrad(var,ex,gradRef):
            step = 1
            if ex[var].isSpecial():
                return (gradRef,step)
            localVarGrad = 0
            localMaxGrad = gradRef
            #Uncomment next line to skip discrete variables
            #return localMaxGrad
            for val in self.domain[var].values:
                localEx = orange.Example(ex)
                localEx[var] = val
                pred = self(localEx,returnDFV = True)[1]
                if abs(gradRef-pred) > localVarGrad:
                    localVarGrad = abs(gradRef-pred)
                    localMaxGrad = pred
            return (localMaxGrad, step)


        def calcContVarGrad(var,ex,gradRef):
            localEx = orange.Example(ex)
            step = 1 #((self.basicStat[var]["max"]-self.basicStat[var]["min"])/(self.getNTrainEx()+1))
            if ex[var].isSpecial():
                return (gradRef, step)
            localEx[var] = ex[var] + step
            localPred = self(localEx,returnDFV = True)[1]
            ResUp = (localPred, step)
            # Print used for testing significance
            #print "%10s%10s" % ( round(step,2),round(localPred,2)),
            #single direction step:
            return ResUp
            # double direction steps: 
            #localEx[var] = ex[var] - step 
            #ResDown = (self(localEx,returnDFV = True)[1], step)
            #if abs(ResUp[0]-gradRef) > abs(ResDown[0]-gradRef):
            #    return ResUp
            #else:
            #    return ResDown

        def calcVarGrad(var,ex,gradRef):
            if attr.varType == orange.VarTypes.Discrete:
                return calcDiscVarGrad(attr.name,ex,gradRef)
            else:
                return calcContVarGrad(attr.name,ex,gradRef)
 
        def insSort(list, IDlist):
            # Sorts both lists simultaniously with reference to list
            # insertion sort algorithm descendent
            for i in range(1, len(list)):
              save = list[i]
              saveName = IDlist[i]
              j = i
              while j > 0 and list[j - 1] < save:
                  list[j] = list[j - 1]
                  IDlist[j] = IDlist[j - 1]
                  j -= 1
              list[j] = save
              IDlist[j] = saveName

        for attr in self.domain.attributes:
            res,step = calcVarGrad(attr.name,ex,gradRef)

            #if self.domain.classVar.varType == orange.VarTypes.Discrete:
                # IMPORTANT: A positive gradient component means the predicted value will be "more strong"
                #            A Nagative one means the predicted value will be weak and may lead towrds the other class value
            #    grad = (abs(res) - abs(gradRef))/step
            #else:
            grad = (res-gradRef)/step
            # Print used for testing significance
            #print  "%10s" % (round(grad,2)),
            if grad == 0:
                continue
            if grad > 0:
                varGradNameUP.append(attr.name)
                varGradValUP.append(abs(grad))
            else:
                varGradNameDOWN.append(attr.name)
                varGradValDOWN.append(abs(grad))
        #print "%6s" % (round(gradRef,2)),
        #Order the vars in terms of importance
        nRet = min(int(nVars),len(self.domain.attributes))
        if absGradient:
            varGradVal = varGradValUP + varGradValDOWN 
            varGradName = varGradNameUP + varGradNameDOWN
            insSort(varGradVal ,varGradName)
            if nVars == 0:
                return varGradName
            else:
                return varGradName[0:nRet]
        else:  
            insSort(varGradValUP ,varGradNameUP)
            insSort(varGradValDOWN ,varGradNameDOWN)
            if nVars == 0:
                return { "UP":varGradNameUP,  "DOWN":varGradNameDOWN}
            else:
                return { "UP":varGradNameUP[0:nRet], "DOWN":varGradNameDOWN[0:nRet] }




