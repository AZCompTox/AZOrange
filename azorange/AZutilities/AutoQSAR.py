import os,time,pickle
import AZOrangeConfig as AZOC
import AZLearnersParamsConfig as OPTconf
import orange
import getAccWOptParam
from AZutilities import miscUtilities
from AZutilities import paramOptUtilities
from trainingMethods import AZorngConsensus
from pprint import pprint
#Create the ML methods (learners)
MLMETHODS = {}
for ML in AZOC.MLMETHODS:
    exec("import "+AZOC.MLMETHODS[ML]["module"])
    MLMETHODS[ML] = eval(AZOC.MLMETHODS[ML]["module"]+"."+AZOC.MLMETHODS[ML]["object"])

print "Available MLMETHODS:",[ml for ml in MLMETHODS]


def getMLStatistics(trainData,savePath = None, queueType = "NoSGE", verbose = 0, logFile = None):
        """
        Loop over all MLMETHODS to get their statistics
        Write to disk the full MLStatistics including the consensus model:
                  Consensus model statistics will be calculated out of the a Consensus model based on MLmethods that are stable (beased on StabilityValue)
        """
        MLStatistics = {}
        learners = {}
        for ml in MLMETHODS:
            learner = MLMETHODS[ml](name = ml)
            if not learner.isCompatible(trainData.domain.classVar):
                print "Ignored learner ",ml," since it's not compatible with this class."
                continue
            learners[ml] = learner
        evaluator = getAccWOptParam.AccWOptParamGetter(data = trainData, learner = learners, paramList = None, nExtFolds = AZOC.QSARNEXTFOLDS, nInnerFolds = AZOC.QSARNINNERFOLDS, queueType = queueType, verbose = verbose, logFile = logFile, resultsFile = None)
        MLStatistics = evaluator.getAcc()

        if savePath and os.path.isdir(os.path.split(savePath)[0]):
            if os.path.isfile(savePath):
                savePath = os.path.splitext(savePath)[0] + "_"+str(time.time()).replace(".","")+os.path.splitext(savePath)[1]
                print "File alerady exists. Saving in:",savePath
            file = open(savePath, "w")
            pickle.dump(MLStatistics, file)
            file.close()
        return MLStatistics



def selectModel(MLStatistics):
        """Return the model with highest R2/CA amongst methods with a stability less than 0.1.
           If no methods is considered stable, select the method with the greatest R2/CA
        """
        bestModelName = None
        bestRes = None
        for modelName in MLStatistics:
            if MLStatistics[modelName]["StabilityValue"] < AZOC.QSARSTABILITYTHRESHOLD:
                valRes = max( MLStatistics[modelName]["R2"], MLStatistics[modelName]["CA"])  # One of them is always None
                if bestRes is None or valRes > bestRes:
                    bestRes = valRes
                    bestModelName = modelName

        if bestModelName is None:
            for modelName in MLStatistics:
                valRes = max( MLStatistics[modelName]["R2"], MLStatistics[modelName]["CA"])  # One of them is always None
                if bestRes is None or valRes > bestRes:
                    bestRes = valRes
                    bestModelName = modelName
        MLMethod = MLStatistics[bestModelName].copy()
        MLMethod["MLMethod"] = bestModelName
        return MLMethod


def buildConsensus(trainData, learners, MLMethods):
        if trainData.domain.classVar.varType == orange.VarTypes.Discrete:
            #Expression: If  CAavg_{POS} ge CAavg_{NEG} -> POS  else -> NEG 
            #    where CAavg_{POS} is the average of classification accuracies of all models predicting POS.
            CLASS0 = str(trainData.domain.classVar.values[0])
            CLASS1 = str(trainData.domain.classVar.values[1])
            exprTest0 = "(0"
            for ml in MLMethods:
                exprTest0 += "+( "+ml+" == "+CLASS0+" )*"+str(MLMethods[ml]["CA"])+" "
            exprTest0 += ")/IF0(sum([False"
            for ml in MLMethods:
                exprTest0 += ", "+ml+" == "+CLASS0+" "
            exprTest0 += "]),1)"
            exprTest1 = exprTest0.replace(CLASS0,CLASS1)
            expression = [exprTest0+" >= "+exprTest1+" -> "+CLASS0," -> "+CLASS1]
        else:
            R2sum = sum([MLMethods[ml]["R2"] for ml in MLMethods])
            expression = "(1 / "+str(R2sum)+") * (0"
            for ml in MLMethods:
                expression += " + "+str(MLMethods[ml]["R2"])+" * " + ml +" "
            expression += ")" 

        consensusLearners = {}
        for learnerName in learners:
            consensusLearners[learnerName] = learners[learnerName]
        
        learner = AZorngConsensus.ConsensusLearner(learners = consensusLearners, expression = expression)
        return learner(trainData)

                    
def buildModel(trainData, MLMethod, queueType = "NoSGE", verbose = 0):
        """
        Buld the methods passed in MLMethods and optimize (if len(MLMethods) == 1)
        if MLMethods is the Consensus (len(MLMethods) > 1) , build each and optimize first all models and after build the consensus!
        """
        learners = {}
        MLMethods = {}
        if "IndividualStatistics"  in MLMethod:                        #It is a consensus
            for ML in MLMethod["IndividualStatistics"]:
                MLMethods[ML] = MLMethod["IndividualStatistics"][ML]
        else:
            MLMethods[MLMethod["MLMethod"]] = MLMethod

        # optimize all MLMethods
        for ML in MLMethods:
            learners[ML] = MLMETHODS[ML](name = ML)

            runPath = miscUtilities.createScratchDir(baseDir = AZOC.NFS_SCRATCHDIR, desc = "AutoQSAR")
            trainData.save(os.path.join(runPath,"trainData.tab"))

            paramOptUtilities.getOptParam(
                learner = learners[ML],
                trainDataFile = os.path.join(runPath,"trainData.tab"),
                useGrid = False,
                verbose = verbose,
                queueType = queueType,
                runPath = runPath,
                nExtFolds = None)

            if not learners[ML].optimized:
                print "ERROR: AutoQSAR: The learner was not optimized."
                return None
            else:
                print "Optimized learner ",learners[ML]           
            miscUtilities.removeDir(runPath)

        #Train the model
        if len(learners) == 1:
            model = learners[learners.keys()[0]](trainData)
        elif len(learners) >= 1:
            model = buildConsensus(trainData,learners,MLMethods)
        else:
            print "ERROR: No Learners were selected!"
            return None

        return model


def getModel(trainData, savePath = None,queueType = "NoSGE", verbose = 0):
        MLStatistics = getMLStatistics(trainData, savePath, queueType = queueType, verbose = verbose)
        MLMethod = selectModel(MLStatistics)
        model = buildModel(trainData, MLMethod, queueType = queueType, verbose = verbose)
        return model


def getStatistics(dataset):
        pass

if __name__ == "__main__":
        data = orange.ExampleTable("./dataReg.tab")
        model = getModel(data, verbose = 3, savePath = "./MLStat_reg.txt")
        print model

        data = orange.ExampleTable("./dataClass.tab")
        model = getModel(data, savePath = "./MLStat_class.txt")
        print model

