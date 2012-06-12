import orange, time, pickle, copy
import sys, traceback
import AZOrangeConfig as AZOC
from AZutilities import paramOptUtilities
from AZutilities import dataUtilities
from trainingMethods import AZorngConsensus
import AZLearnersParamsConfig
from AZutilities import evalUtilities
from AZutilities import miscUtilities
import orngTest, orngStat
import os,random
from pprint import pprint
import statc



class UnbiasedAccuracyGetter():
    def __init__(self, **kwds):
        self.verbose = 0
        self.logFile = None
        self.resultsFile = None
        self.nExtFolds = 5
        self.nInnerFolds = 5
        self.data = None
        self.learner = None
        self.paramList = None
        self.queueType = "NoSGE"
        self.responseType = None
        self.fixedParams = {} 
        self.testAttrFilter = None
        self.testFilterVal = None
        self.sampler = dataUtilities.SeedDataSampler
        # Append arguments to the __dict__ member variable 
        self.__dict__.update(kwds)
        self.learnerName = ""

        self.trainBiasMask = []
        if self.testFilterVal and type(self.testFilterVal) == list and self.testAttrFilter and type(self.testAttrFilter) == str and self.testAttrFilter in self.data.domain:
            self.useVarCtrlCV = True
            for ex in self.data:
                if ex[self.testAttrFilter].value in self.testFilterVal: # Compound selected to be allowed in the test set
                    self.trainBiasMask.append(0)
                else:                                                 # Compound to not include in the test set. Always to be shifted to the train
                    self.trainBiasMask.append(1)
            self.data = dataUtilities.attributeDeselectionData(self.data, [self.testAttrFilter]) 
        else:
            self.useVarCtrlCV = False
            self.testAttrFilter = None
            self.testFilterVal = None



    def __checkTrainData(self, trainData, stopIfError = True):
        """
            Checks if the train data complies with rules:
                -Have at least 20 examples
                -If classifier:
                        Have at least 10 examples for each class value
        """
        error = False
        logStr = "Too few compounds to build a QSAR model"
        if len(trainData) < 20:
            error = True
        elif self.responseType == "Classification":
            valuesCount = dict.fromkeys([str(val) for val in trainData.domain.classVar.values], 0)
            for ex in trainData:
                valuesCount[str(ex.getclass().value)] +=1
            for val in valuesCount:
                if valuesCount[val] < 10:
                    logStr += "\n    " + val 
                    error = True
        if error:
            if stopIfError:
                self.__log(logStr)
                raise Exception(logStr)
            else:
                return False
        else:
            return True
        


    def __writeResults(self, statObj):
        if self.resultsFile and os.path.isdir(os.path.split(self.resultsFile)[0]):
            file = open(self.resultsFile, "w")
            pickle.dump(statObj, file)
            file.close()


    def __log(self, text):
        """Adds a new line (what's in text) to the logFile"""
        textOut = str(time.asctime()) + ": " +text
        if self.logFile and os.path.isdir(os.path.split(self.logFile)[0]):
            file = open(self.logFile, "a")
            file.write(textOut+"\n")
            file.close()
        else:
            print textOut

    def __areInputsOK(self):
        if not self.learner or (not self.paramList and type(self.learner)!=dict) or not self.nExtFolds or not self.nInnerFolds or not self.data or not self.sampler:
            self.__log("   Missing configuration in UnbiasedAccuracyGetter object")
            return False
        if not self.data.domain.classVar:
            self.__log("   The data has no Class!")
            return False
        if self.queueType not in ["NoSGE", "batch.q", "quick.q"]:
            self.__log("   Invalid queueType")
            return False
        if not len(self.data):
            self.__log("   Data is empty")
            return False
        if len(self.data)/self.nExtFolds < 1:
            self.__log("   Too few examples for " + str(self.nExtFolds) + "folds.")
            return False
        if type(self.learner)==dict and self.paramList:
            self.__log("   WARNING: A set of learners was provided, and therefore the paramList will be ignored. Default paramneters will be optimized instead.")
            self.paramList = None
        elif self.paramList:
            try:
                # Find the name of the Learner
                self.learnerName = str(self.learner.__class__)[:str(self.learner.__class__).rfind("'")].split(".")[-1]
            except:
                self.__log("   Couldn't find the Learner Name of: "+ str(a))
                return False

            if not hasattr(AZLearnersParamsConfig,self.learnerName):
                self.__log("   The learner '"+str(self.learnerName)+"' is not compatible with the optimizer")
                return False
            parsAPI = AZLearnersParamsConfig.API(self.learnerName)
            for par in self.paramList: 
                if par not in parsAPI.getParameterNames():
                    self.__log("   Parameter "+str(par)+" does not exist for the learner "+str(self.learnerName))
                    return False
        return True

    def createStatObj(self, results=None, exp_pred=None, nTrainCmpds=None, nTestCmpds=None, responseType=None, nExtFolds=None, userAlert = "", labels = None):
        #Initialize res (statObj) for statistic results
        res = {}
        # Classification
        res["CA"] = None
        res["CM"] = None
        res["MCC"] = None
        res["labels"] = None
        #Regression
        res["Q2"] = None
        res["RMSE"] = None
        #Both
        res["StabilityValue"] = None
        res["userAlert"] = userAlert
        res["selected"] = False
        res["stable"] = False
        res["responseType"] = False
        res["foldStat"] = {
                "nTrainCmpds": None,
                "nTestCmpds": None,
                #Regression
                "Q2"   : None,
                "RMSE" : None,
                #Classification
                "CM"   : None,
                "CA"   : None,
                "MCC"  : None }
        if results is None or exp_pred is None or responseType is None or nExtFolds is None or nTestCmpds is None or nTrainCmpds is None:
            return res 
        res["responseType"] = responseType
        #Calculate the (Q2, RMSE) or (CM, CA) results depending on Classification or regression
        if responseType == "Classification":
            #Compute CA
            res["CA"] = sum(r[0] for r in results) / nExtFolds
            #Compute CM
            res["CM"] = copy.deepcopy(results[0][1])                      # Get the first ConfMat
            for r in results[1:]:
                for Lidx,line in enumerate(r[1]):
                    for idx,val in enumerate(line):
                        res["CM"][Lidx][idx] = res["CM"][Lidx][idx] + val   #Add each same ConfMat position
            #Compute MCC 
            res["MCC"] = evalUtilities.calcMCC(res["CM"])
            res["labels"] = labels
            #Compute foldStat
            res["foldStat"]["nTrainCmpds"] = [n for n in nTrainCmpds]
            res["foldStat"]["nTestCmpds"] = [n for n in nTestCmpds]
            res["foldStat"]["CA"] = [r[0] for r in results]
            res["foldStat"]["CM"] = [r[1] for r in results]
            res["foldStat"]["MCC"] = [evalUtilities.calcMCC(r[1]) for r in results]
            #Compute Stability
            res["StabilityValue"] = evalUtilities.stability(res["foldStat"]["CA"])
        else:
            #compute Q2
            res["Q2"] = evalUtilities.calcRsqrt(exp_pred)
            #compute RMSE
            res["RMSE"] = evalUtilities.calcRMSE(exp_pred)
            #Compute foldStat
            res["foldStat"]["nTrainCmpds"] = [n for n in nTrainCmpds]
            res["foldStat"]["nTestCmpds"] = [n for n in nTestCmpds]
            res["foldStat"]["RMSE"] = [r[0] for r in results]
            res["foldStat"]["Q2"] = [r[1] for r in results]
            #Compute Stability value
            res["StabilityValue"] = evalUtilities.stability(res["foldStat"]["Q2"])
        #Evaluate stability of ML
        StabilityValue = res["StabilityValue"]
        if StabilityValue is not None:
            if responseType == "Classification":
                if statc.mean(res["foldStat"]["nTestCmpds"]) > 50:
                    stableTH = AZOC.QSARSTABILITYTHRESHOLD_CLASS_L
                else:
                    stableTH = AZOC.QSARSTABILITYTHRESHOLD_CLASS_H
            else:
                if statc.mean(res["foldStat"]["nTestCmpds"]) > 50:
                    stableTH = AZOC.QSARSTABILITYTHRESHOLD_REG_L
                else:
                    stableTH = AZOC.QSARSTABILITYTHRESHOLD_REG_H
            if StabilityValue < stableTH:   # Select only stable models
                res["stable"] = True

        return res
        
    def getAcc(self, callBack = None, callBackWithFoldModel = None):
        """ For regression problems, it returns the RMSE and the Q2 
            For Classification problems, it returns CA and the ConfMat
            The return is made in a Dict: {"RMSE":0.2,"Q2":0.1,"CA":0.98,"CM":[[TP, FP],[FN,TN]]}
            For the EvalResults not supported for a specific learner/datase, the respective result will be None

            if the learner is a dict {"LearnerName":learner, ...} the results will be a dict with results for all Learners and for a consensus
                made out of those that were stable

            It some error occurred, the respective values in the Dict will be None
        """
        self.__log("Starting Calculating MLStatistics")
        statistics = {}
        if not self.__areInputsOK():
            return None
        # Set the response type
        self.responseType =  self.data.domain.classVar.varType == orange.VarTypes.Discrete and "Classification"  or "Regression"
        self.__log("  "+str(self.responseType))

        #Create the Train and test sets
        DataIdxs = dataUtilities.SeedDataSampler(self.data, self.nExtFolds) 
        
        #Fix the Indexes based on self.trainBiasMask
        # (0s) represents the train set  (1s) represents the test set
        if self.useVarCtrlCV:
            nShifted = [0] * len(DataIdxs)
            self.__log("Number of examples marked for testing: "+str(len(self.trainBiasMask) - sum(self.trainBiasMask)))
            for idx,isTrainBias in enumerate(self.trainBiasMask):
                if isTrainBias:
                    for foldIdx in range(len(DataIdxs)):
                        if DataIdxs[foldIdx][idx]:
                            DataIdxs[foldIdx][idx] = 0  
                            nShifted[foldIdx] +=1
            for idx,shift in enumerate(nShifted):
                self.__log("In fold "+str(idx)+", "+str(shift)+" examples were shifted to the train set.")

        #Var for saving each Fols result
        optAcc = {}
        results = {}
        exp_pred = {}
        nTrainEx = {}
        nTestEx = {}
        
        #Set a dict of learners
        MLmethods = {}
        if type(self.learner) == dict:
            for ml in self.learner:
                MLmethods[ml] = self.learner[ml]
        else:
            MLmethods[self.learner.name] = self.learner

        models={}
        self.__log("Calculating Statistics for MLmethods:")
        self.__log("  "+str([x for x in MLmethods]))

        #Check data in advance so that, by chance, it will not faill at the last fold!
        for foldN in range(self.nExtFolds):
            trainData = self.data.select(DataIdxs[foldN],negate=1)
            self.__checkTrainData(trainData)

        #Optional!!
        # Order Learners so that PLS is the first
        sortedML = [ml for ml in MLmethods]
        if "PLS" in sortedML:
            sortedML.remove("PLS")
            sortedML.insert(0,"PLS")

        stepsDone = 0
        nTotalSteps = len(sortedML) * self.nExtFolds  
        for ml in sortedML:
          self.__log("    > "+str(ml)+"...")
          try:
            #Var for saving each Fols result
            results[ml] = []
            exp_pred[ml] = []
            models[ml] = []
            nTrainEx[ml] = []
            nTestEx[ml] = []
            optAcc[ml] = []
            logTxt = "" 
            for foldN in range(self.nExtFolds):
                if type(self.learner) == dict:
                    self.paramList = None

                trainData = self.data.select(DataIdxs[foldN],negate=1)
                testData = self.data.select(DataIdxs[foldN])
                nTrainEx[ml].append(len(trainData))
                nTestEx[ml].append(len(testData))
                #Test if trainsets inside optimizer will respect dataSize criterias.
                #  if not, don't optimize, but still train the model
                dontOptimize = False
                if self.responseType != "Classification" and (len(trainData)*(1-1.0/self.nInnerFolds) < 20):
                    dontOptimize = True
                else:                      
                    tmpDataIdxs = dataUtilities.SeedDataSampler(trainData, self.nInnerFolds)
                    tmpTrainData = trainData.select(tmpDataIdxs[0],negate=1)
                    if not self.__checkTrainData(tmpTrainData, False):
                        dontOptimize = True

                if dontOptimize:
                    logTxt += "       Fold "+str(foldN)+": Too few compounds to optimize model hyper-parameters\n"
                    self.__log(logTxt)
                    if trainData.domain.classVar.varType == orange.VarTypes.Discrete:
                        res = orngTest.crossValidation([MLmethods[ml]], trainData, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible, randomGenerator = random.randint(0, 100))
                        CA = evalUtilities.CA(res)[0]
                        optAcc[ml].append(CA)
                    else:
                        res = orngTest.crossValidation([MLmethods[ml]], trainData, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible, randomGenerator = random.randint(0, 100))
                        R2 = evalUtilities.R2(res)[0]
                        optAcc[ml].append(R2)
                else:
                    runPath = miscUtilities.createScratchDir(baseDir = AZOC.NFS_SCRATCHDIR, desc = "AccWOptParam", seed = id(trainData))
                    trainData.save(os.path.join(runPath,"trainData.tab"))

                    tunedPars = paramOptUtilities.getOptParam(
                        learner = MLmethods[ml], 
                        trainDataFile = os.path.join(runPath,"trainData.tab"), 
                        paramList = self.paramList, 
                        useGrid = False, 
                        verbose = self.verbose, 
                        queueType = self.queueType, 
                        runPath = runPath, 
                        nExtFolds = None, 
                        nFolds = self.nInnerFolds,
                        logFile = self.logFile,
                        getTunedPars = True,
                        fixedParams = self.fixedParams)
                    if not MLmethods[ml] or not MLmethods[ml].optimized:
                        self.__log("       WARNING: GETACCWOPTPARAM: The learner "+str(ml)+" was not optimized.")
                        self.__log("                It will be ignored")
                        #self.__log("                It will be set to default parameters")
                        self.__log("                    DEBUG can be done in: "+runPath)
                        #Set learner back to default 
                        #MLmethods[ml] = MLmethods[ml].__class__()
                        raise Exception("The learner "+str(ml)+" was not optimized.")
                    else:
                        if trainData.domain.classVar.varType == orange.VarTypes.Discrete:
                            optAcc[ml].append(tunedPars[0])
                        else:
                            res = orngTest.crossValidation([MLmethods[ml]], trainData, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible, randomGenerator = random.randint(0, 100))
                            R2 = evalUtilities.R2(res)[0]
                            optAcc[ml].append(R2)

                        miscUtilities.removeDir(runPath) 
                #Train the model
                model = MLmethods[ml](trainData)
                models[ml].append(model)
                #Test the model
                if self.responseType == "Classification":
                    results[ml].append((evalUtilities.getClassificationAccuracy(testData, model), evalUtilities.getConfMat(testData, model) ) )
                else:
                    local_exp_pred = []
                    for ex in testData:
                        local_exp_pred.append((ex.getclass(), model(ex)))
                    results[ml].append((evalUtilities.calcRMSE(local_exp_pred), evalUtilities.calcRsqrt(local_exp_pred) ) )
                    #Save the experimental value and correspondent predicted value
                    exp_pred[ml] += local_exp_pred
                if callBack:
                     stepsDone += 1
                     if not callBack((100*stepsDone)/nTotalSteps): return None
                if callBackWithFoldModel:
                    callBackWithFoldModel(model) 

            res = self.createStatObj(results[ml], exp_pred[ml], nTrainEx[ml], nTestEx[ml],self.responseType, self.nExtFolds, logTxt, labels = hasattr(self.data.domain.classVar,"values") and list(self.data.domain.classVar.values) or None )
            if self.verbose > 0: 
                print "UnbiasedAccuracyGetter!Results  "+ml+":\n"
                pprint(res)
            if not res:
                raise Exception("No results available!")
            statistics[ml] = copy.deepcopy(res)
            self.__writeResults(statistics)
            self.__log("       OK")
          except:
            self.__log("       Learner "+str(ml)+" failed to create/optimize the model!")
            error = str(sys.exc_info()[0]) +" "+\
                        str(sys.exc_info()[1]) +" "+\
                        str(traceback.extract_tb(sys.exc_info()[2]))
            self.__log(error)
 
            res = self.createStatObj()
            statistics[ml] = copy.deepcopy(res)
            self.__writeResults(statistics)

        if not statistics or len(statistics) < 1:
            self.__log("ERROR: No statistics to return!")
            return None
        elif len(statistics) > 1:
            #We still need to build a consensus model out of the stable models 
            #   ONLY if there are more that one model stable!
            #   When only one or no stable models, build a consensus based on all models
            consensusMLs={}
            for modelName in statistics:
                StabilityValue = statistics[modelName]["StabilityValue"]
                if StabilityValue is not None and statistics[modelName]["stable"]:
                    consensusMLs[modelName] = copy.deepcopy(statistics[modelName])

            self.__log("Found "+str(len(consensusMLs))+" stable MLmethods out of "+str(len(statistics))+" MLmethods.")

            if len(consensusMLs) <= 1:   # we need more models to build a consensus!
                consensusMLs={}
                for modelName in statistics:
                    consensusMLs[modelName] = copy.deepcopy(statistics[modelName])
 
            if len(consensusMLs) >= 2:
                #Var for saving each Fols result
                Cresults = []
                Cexp_pred = []
                CnTrainEx = []
                CnTestEx = []
                self.__log("Calculating the statistics for a Consensus model based on "+str([ml for ml in consensusMLs]))
                for foldN in range(self.nExtFolds):
                    if self.responseType == "Classification":
                        CLASS0 = str(self.data.domain.classVar.values[0])
                        CLASS1 = str(self.data.domain.classVar.values[1])
                        # exprTest0
                        exprTest0 = "(0"
                        for ml in consensusMLs:
                            exprTest0 += "+( "+ml+" == "+CLASS0+" )*"+str(optAcc[ml][foldN])+" "
                        exprTest0 += ")/IF0(sum([False"
                        for ml in consensusMLs:
                            exprTest0 += ", "+ml+" == "+CLASS0+" "
                        exprTest0 += "]),1)"
                        # exprTest1
                        exprTest1 = "(0"
                        for ml in consensusMLs:
                            exprTest1 += "+( "+ml+" == "+CLASS1+" )*"+str(optAcc[ml][foldN])+" "
                        exprTest1 += ")/IF0(sum([False"
                        for ml in consensusMLs:
                            exprTest1 += ", "+ml+" == "+CLASS1+" "
                        exprTest1 += "]),1)"
                        # Expression
                        expression = [exprTest0+" >= "+exprTest1+" -> "+CLASS0," -> "+CLASS1]
                    else:
                        Q2sum = sum([optAcc[ml][foldN] for ml in consensusMLs])
                        expression = "(1 / "+str(Q2sum)+") * (0"
                        for ml in consensusMLs:
                            expression += " + "+str(optAcc[ml][foldN])+" * "+ml+" "
                        expression += ")"

                    testData = self.data.select(DataIdxs[foldN])
                    CnTestEx.append(len(testData))
                    consensusClassifiers = {}
                    for learnerName in consensusMLs:
                        consensusClassifiers[learnerName] = models[learnerName][foldN]

                    model = AZorngConsensus.ConsensusClassifier(classifiers = consensusClassifiers, expression = expression)     
                    CnTrainEx.append(model.NTrainEx)
                    #Test the model
                    if self.responseType == "Classification":
                        Cresults.append((evalUtilities.getClassificationAccuracy(testData, model), evalUtilities.getConfMat(testData, model) ) )
                    else:
                        local_exp_pred = []
                        for ex in testData:
                            local_exp_pred.append((ex.getclass(), model(ex)))
                        Cresults.append((evalUtilities.calcRMSE(local_exp_pred), evalUtilities.calcRsqrt(local_exp_pred) ) )
                        #Save the experimental value and correspondent predicted value
                        Cexp_pred += local_exp_pred

                res = self.createStatObj(Cresults, Cexp_pred, CnTrainEx, CnTestEx, self.responseType, self.nExtFolds, labels = hasattr(self.data.domain.classVar,"values") and list(self.data.domain.classVar.values) or None )
                statistics["Consensus"] = copy.deepcopy(res)
                statistics["Consensus"]["IndividualStatistics"] = copy.deepcopy(consensusMLs)
                self.__writeResults(statistics)
            self.__log("Returned multiple ML methods statistics.")
            return statistics
                 
        #By default return the only existing statistics!
        self.__writeResults(statistics)
        self.__log("Returned only one ML method statistics.")
        return statistics[statistics.keys()[0]]




