import orange
import AZOrangeConfig as AZOC
from AZutilities import paramOptUtilities
from AZutilities import dataUtilities
from trainingMethods import AZorngConsensus
import AZLearnersParamsConfig
from AZutilities import evalUtilities
from AZutilities import miscUtilities
import orngTest, orngStat
import os
from pprint import pprint

class AccWOptParamGetter():
    def __init__(self, **kwds):
        self.verbose = 0
        self.nExtFolds = 5
        self.nInnerFolds = 5
        self.data = None
        self.learner = None
        self.paramList = None
        self.queueType = "NoSGE"
        self.sampler = dataUtilities.SeedDataSampler
        # Append arguments to the __dict__ member variable 
        self.__dict__.update(kwds)
        self.learnerName = ""

    def __areInputsOK(self):
        if not self.learner or (not self.paramList and type(self.learner)!=dict) or not self.nExtFolds or not self.nInnerFolds or not self.data or not self.sampler:
            print "Missing configuration in AccWOptParamGetter object"
            return False
        if not self.data.domain.classVar:
            print "The data has no Class!"
            return False
        if self.queueType not in ["NoSGE", "batch.q", "quick.q"]
            print "Invalid queueType"
            return False
        if not len(self.data):
            print "Data is empty"
            return False
        if len(self.data)/self.nExtFolds < 1:
            print "Too few examples for ", self.nExtFolds, "folds."
            return False
        if type(self.learner)==dict and self.paramList:
            print "WARNING: A set of learners was provided, and therefore the paramList will be ignored. Default paramneters will be oprtimized instead."
            self.paramList = None
        elif self.paramList:
            try:
                # Find the name of the Learner
                self.learnerName = str(self.learner.__class__)[:str(self.learner.__class__).rfind("'")].split(".")[-1]
            except:
                print "Couldn't find the Learner Name of: ", str(a)
                return False

            if not hasattr(AZLearnersParamsConfig,self.learnerName):
                print "The learner '"+self.learnerName+"' is not compatible with the optimizer"
                return False
            parsAPI = AZLearnersParamsConfig.API(self.learnerName)
            for par in self.paramList: 
                if par not in parsAPI.getParameterNames():
                    print "Parameter "+par+" does not exist for the learner "+self.learnerName
                    return False
        return True

    def createStatObj(self, results=None, exp_pred=None, responseType=None, nExtFolds=None):
        #Initialize res (statObj) for statistic results
        res = {}
        # Classification
        res["CA"] = None
        res["CM"] = None
        res["MCC"] = None
        #Regression
        res["R2"] = None
        res["RMSE"] = None
        #Both
        res["StabilityValue"] = None
        res["foldStat"] = {
                #Regression
                "R2"   : None,
                "RMSE" : None,
                #Classification
                "CM"   : None,
                "CA"   : None,
                "MCC"  : None }
        if results is None or exp_pred is None or responseType is None or nExtFolds is None:
            return res 
        #Calculate the (R2, RMSE) or (CM, CA) results depending on Classification or regression
        if responseType == "Classification":
            #Compute CA
            res["CA"] = sum(r[0] for r in results) / self.nExtFolds
            #Compute CM
            res["CM"] = results[0][1]                      # Get the first ConfMat
            for r in results[1:]:
                for Lidx,line in enumerate(r[1]):
                    for idx,val in enumerate(line):
                        res["CM"][Lidx][idx] = res["CM"][Lidx][idx] + val   #Add each same ConfMat position
            #Compute MCC 
            res["MCC"] = evalUtilities.calcMCC(res["CM"])
            #Compute foldStat
            res["foldStat"]["CA"] = [r[0] for r in results]
            res["foldStat"]["CM"] = [r[1] for r in results]
            res["foldStat"]["MCC"] = [evalUtilities.calcMCC(r[1]) for r in results]
            #Compute Stability
            res["StabilityValue"] = evalUtilities.stability(res["foldStat"]["CA"])
        else:
            #compute R2
            res["R2"] = evalUtilities.calcRsqrt(exp_pred)
            #compute RMSE
            res["RMSE"] = evalUtilities.calcRMSE(exp_pred)
            #Compute foldStat
            res["foldStat"]["RMSE"] = [r[0] for r in results]
            res["foldStat"]["R2"] = [r[1] for r in results]
            #Compute Stability
            res["StabilityValue"] = evalUtilities.stability(res["foldStat"]["R2"])
        return res
        
    def getAcc(self):
        """ For regression problems, it returns the RMSE and the R2 
            For Classification problems, it returns CA and the ConfMat
            The return is made in a Dict: {"RMSE":0.2,"R2":0.1,"CA":0.98,"CM":[[TP, FP],[FN,TN]]}
            For the EvalResults not supported for a specific learner/datase, the respective result will be None

            if the learner is a dict {"LearnerName":learner, ...} the results will be a dict with results for all Learners and for a consensus
                made out of those that were stable

            It some error occurred, the respective values in the Dict will be None
        """
        statistics = {}
        if not self.__areInputsOK():
            return None
        # Set the response type
        responseType =  self.data.domain.classVar.varType == orange.VarTypes.Discrete and "Classification"  or "Regression"
        

        #Create the Train and test sets
        DataIdxs = dataUtilities.SeedDataSampler(self.data, self.nExtFolds) 
        
        #Var for saving each Fols result
        results = {}
        exp_pred = {}
        
        #Set a dict of learners
        MLmethods = {}
        if type(self.learner) == dict:
            for ml in self.learner:
                MLmethods[ml] = self.learner[ml]
        else:
            MLmethods[self.learner.name] = self.learner

        models={}
        for ml in MLmethods:
          try:
            #Var for saving each Fols result
            results[ml] = []
            exp_pred[ml] = []
            models[ml] = []
            for foldN in range(self.nExtFolds):
                if type(self.learner) == dict:
                    self.paramList = None

                trainData = self.data.select(DataIdxs[foldN],negate=1)
                runPath = miscUtilities.createScratchDir(desc = "AccWOptParam")
                trainData.save(os.path.join(runPath,"trainData.tab"))
                testData = self.data.select(DataIdxs[foldN])

                paramOptUtilities.getOptParam(
                    learner = MLmethods[ml], 
                    trainDataFile = os.path.join(runPath,"trainData.tab"), 
                    paramList = self.paramList, 
                    useGrid = False, 
                    verbose = self.verbose, 
                    queueType = self.queueType, 
                    runPath = runPath, 
                    nExtFolds = None, 
                    nFolds = self.nInnerFolds)
                if not MLmethods[ml].optimized:
                    print "The learner "+str(ml)+" was not optimized."
                    raise Exception("The learner "+str(ml)+" was not optimized.")
                miscUtilities.removeDir(runPath) 
                #Train the model
                model = MLmethods[ml](trainData)
                models[ml].append(model)
                #Test the model
                if responseType == "Classification":
                    results[ml].append((evalUtilities.getClassificationAccuracy(testData, model), evalUtilities.getConfMat(testData, model) ) )
                else:
                    local_exp_pred = []
                    for ex in testData:
                        local_exp_pred.append((ex.getclass(), model(ex)))
                    results[ml].append((evalUtilities.calcRMSE(local_exp_pred), evalUtilities.calcRsqrt(local_exp_pred) ) )
                    #Save the experimental value and correspondent predicted value
                    exp_pred[ml] += local_exp_pred
   
            res = self.createStatObj(results[ml], exp_pred[ml], responseType, self.nExtFolds)
            if self.verbose > 0: 
                print "AccWOptParamGetter!Results  "+ml+":\n"
                pprint(res)
            if not res:
                raise Exception("No results available!")
            statistics[ml] = res.copy()
          except:
            print "Learner "+ml+" failed to optimize!"
            res = self.createStatObj()
            statistics[ml] = res.copy()

        if len(statistics) < 1:
            print "ERROR: No statistics to return!"
        elif len(statistics) > 1:
            #We still need to build a consensus model out of the stable models 
            #   ONLY if there are more that one model stable!
            stableML={}
            for modelName in statistics:
                if statistics[modelName]["StabilityValue"] < AZOC.QSARSTABILITYTHRESHOLD:   # Select only stable models
                    stableML[modelName] = statistics[modelName].copy()
            if len(stableML) >= 2:
                if responseType == "Classification":
                    CLASS0 = str(self.data.domain.classVar.values[0])
                    CLASS1 = str(self.data.domain.classVar.values[1])
                    exprTest0 = "(0"
                    for ml in stableML:
                        exprTest0 += "+( "+ml+" == "+CLASS0+" )*"+str(stableML[ml]["CA"])+" "
                    exprTest0 += ")/IF0(sum([False"
                    for ml in stableML:
                        exprTest0 += ", "+ml+" == "+CLASS0+" "
                    exprTest0 += "]),1)"
                    exprTest1 = exprTest0.replace(CLASS0,CLASS1)
                    expression = [exprTest0+" >= "+exprTest1+" -> "+CLASS0," -> "+CLASS1]
                else:
                    R2sum = sum([stableML[ml]["R2"] for ml in stableML])
                    expression = "(1 / "+str(R2sum)+") * (0"
                    for ml in stableML:
                        expression += " + "+str(stableML[ml]["R2"])+" * "+ml+" "
                    expression += ")"

                #Var for saving each Fols result
                Cresults = []
                Cexp_pred = []
                for foldN in range(self.nExtFolds):
                    testData = self.data.select(DataIdxs[foldN])
                    consensusClassifiers = {}
                    for learnerName in stableML:
                        consensusClassifiers[learnerName] = models[learnerName][foldN]

                    model = AZorngConsensus.ConsensusClassifier(classifiers = consensusClassifiers, expression = expression)     
                    #Test the model
                    if responseType == "Classification":
                        Cresults.append((evalUtilities.getClassificationAccuracy(testData, model), evalUtilities.getConfMat(testData, model) ) )
                    else:
                        local_exp_pred = []
                        for ex in testData:
                            local_exp_pred.append((ex.getclass(), model(ex)))
                        Cresults.append((evalUtilities.calcRMSE(local_exp_pred), evalUtilities.calcRsqrt(local_exp_pred) ) )
                        #Save the experimental value and correspondent predicted value
                        Cexp_pred += local_exp_pred

                res = self.createStatObj(Cresults, Cexp_pred, responseType, self.nExtFolds)
                statistics["Consensus"] = res.copy()
                statistics["Consensus"]["IndividualStatistics"] = stableML.copy()
            return statistics
                 
        #By default return the only existing statistics!
        return statistics[statistics.keys()[0]]




