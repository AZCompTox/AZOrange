import orange
from AZutilities import paramOptUtilities
from AZutilities import dataUtilities
import AZLearnersParamsConfig
from AZutilities import evalUtilities
from AZutilities import miscUtilities
import orngTest, orngStat
import os

class AccWOptParamGetter():
    def __init__(self, **kwds):
        self.verbose = 0
        self.nExtFolds = 5
        self.nInnerFolds = 5
        self.data = None
        self.learner = None
        self.paramList = None
        self.sampler = dataUtilities.SeedDataSampler
        # Append arguments to the __dict__ member variable 
        self.__dict__.update(kwds)
        self.learnerName = ""

    def __areInputsOK(self):
        if not self.paramList or not self.nExtFolds or not self.nInnerFolds or not self.data or not self.learner or not self.sampler:
            print "Missing configuration in AccWOptParamGetter object"
            return False
        if not self.data.domain.classVar:
            print "The data has no Class!"
            return False
        if not len(self.data):
            print "Data is empty"
            return False
        if len(self.data)/self.nExtFolds < 1:
            print "Too few examples for ", self.nExtFolds, "folds."
            return False
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
                print "Parameter "+par+" does not exist for the lraener "+self.learnerName
                return False
        return True

    def getAcc(self):
        """ For regression problems, it returns the RMSE and the R2 
            For Classification problems, it returns CA and the ConfMat
            The return is made in a Dict: {"RMSE":0.2,"R2":0.1,"CA":0.98,"CM":[[TP, FP],[FN,TN]]}
            For the EvalResults not supported for a specific learner/datase, the respective result will be None
            It some error occurred, the respective values in the Dict will be None
        """
        if not self.__areInputsOK():
            return None
        res = {"RMSE":None,"R2":None,"CA":None,"CM":None}
        # Set the response type
        responseType =  self.data.domain.classVar.varType == orange.VarTypes.Discrete and "Classification"  or "Regression"
        

        #Create the Train and test sets
        DataIdxs = dataUtilities.SeedDataSampler(self.data, self.nExtFolds) 
        
        #Var for saving each Fols result
        results = []

        for foldN in range(self.nExtFolds):
            trainData = self.data.select(DataIdxs[foldN],negate=1)
            runPath = miscUtilities.createScratchDir(desc = "AccWOptParam")
            trainData.save(os.path.join(runPath,"trainData.tab"))
            testData = self.data.select(DataIdxs[foldN])

            paramOptUtilities.optimizeSelectedParam(
                learner = self.learner, 
                learnerName = self.learnerName,
                trainDataFile = os.path.join(runPath,"trainData.tab"), 
                paramList = self.paramList, 
                responseType = responseType, 
                grid = False, 
                useGrid = False, 
                verbose = 0, 
                queueType = "batch.q", 
                runPath = runPath, 
                nExtFolds = None, 
                nFolds = self.nInnerFolds)
            if not self.learner.optimized:
                print "The learner was not optimized."
                return None
            miscUtilities.removeDir(runPath) 
            #Train the model
            model = self.learner(trainData)
            #Test teh model
            if responseType == "Classification":
                results.append((evalUtilities.getClassificationAccuracy(testData, model), evalUtilities.getConfMat(testData, model) ) )
            else:
                results.append((evalUtilities.getRMSE(testData, model), evalUtilities.getRsqrt(testData, model) ) )

        #Calculate the average of results
        #Compute the first result (CA or RMSE)
        if responseType == "Classification":
            resName = "CA"
        else:
            resName = "RMSE"
        res[resName] = 0.0
        for r in results:
            res[resName] += r[0]
        res[resName] = res[resName] / self.nExtFolds
        #Compute the second result (ConfMat or R2)
        if responseType == "Classification":
            res["CM"] = results[0][1]                      # Get the first ConfMat
            for r in results[1:]:
                for Lidx,line in enumerate(r[1]):
                    for idx,val in enumerate(line):
                        res["CM"][Lidx][idx] = res["CM"][Lidx][idx] + val   #Add each same ConfMat position
        else:
            res["R2"] = 0.0
            for r in results:
                res["R2"] += r[1]
            res["R2"] = res["R2"] / self.nExtFolds

        if self.verbose > 0: print "AccWOptParamGetter!Results: ",results, "\n res = ",res
        return res


