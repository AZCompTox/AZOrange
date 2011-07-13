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

def log(logFile, text):
        """Adds a new line (what's in text) to the logFile"""
        textOut = str(time.asctime()) + ": " +text
        if logFile and os.path.isdir(os.path.split(logFile)[0]):
            file = open(logFile, "a")
            file.write(textOut+"\n")
            file.close()
        else:
            print textOut


def getMLStatistics(trainData, savePath = None, queueType = "NoSGE", verbose = 0, logFile = None):
        """
        Loop over all MLMETHODS to get their statistics
        Write to disk the full MLStatistics including the consensus model:
                  Consensus model statistics will be calculated out of the a Consensus model based on MLmethods that are stable (beased on StabilityValue)
        """
        log(logFile, "Running getMLStatistics...")
        MLStatistics = {}
        learners = {}
        for ml in MLMETHODS:
            learner = MLMETHODS[ml](name = ml)
            if not learner.isCompatible(trainData.domain.classVar):
                print "Ignored learner ",ml," since it's not compatible with this class."
                continue
            learners[ml] = learner
        evaluator = getAccWOptParam.AccWOptParamGetter(data = trainData, learner = learners, paramList = None, nExtFolds = AZOC.QSARNEXTFOLDS, nInnerFolds = AZOC.QSARNINNERFOLDS, queueType = queueType, verbose = verbose, logFile = logFile, resultsFile = savePath)
        MLStatistics = evaluator.getAcc()

        if savePath and os.path.isdir(os.path.split(savePath)[0]):
            file = open(savePath, "w")
            pickle.dump(MLStatistics, file)
            file.close()
            log(logFile, "MLStatistics saved to: "+savePath)
        return MLStatistics



def selectModel(MLStatistics, logFile = None):
        """Return the model with highest R2/CA amongst methods with a stability less than 0.1.
           If no methods is considered stable, select the method with the greatest R2/CA
        """
        log(logFile, "Selecting MLmethod...")
        bestModelName = None
        bestRes = None
        bestStableVal = None
        for modelName in MLStatistics:
            if MLStatistics[modelName]["StabilityValue"] < AZOC.QSARSTABILITYTHRESHOLD:
                valRes = max( MLStatistics[modelName]["R2"], MLStatistics[modelName]["CA"])  # One of them is always None
                if bestRes is None or valRes > bestRes:
                    bestRes = valRes
                    bestModelName = modelName
                    bestStableVal = MLStatistics[modelName]["StabilityValue"]
                elif valRes == bestRes and MLStatistics[modelName]["StabilityValue"] > bestStableVal:
                    bestRes = valRes
                    bestModelName = modelName
                    bestStableVal = MLStatistics[modelName]["StabilityValue"]
                    


        if bestModelName is None:
            log(logFile, "  No stable models found! Selecting the one with best result still...")
            for modelName in MLStatistics:
                valRes = max( MLStatistics[modelName]["R2"], MLStatistics[modelName]["CA"])  # One of them is always None
                if bestRes is None or valRes > bestRes:
                    bestRes = valRes
                    bestModelName = modelName
            log(logFile, "  Selected the non-stable MLmethod: " + bestModelName)
        else:
            log(logFile, "  Selected the stable MLmethod: " + bestModelName)
        MLMethod = MLStatistics[bestModelName].copy()
        MLMethod["MLMethod"] = bestModelName
        return MLMethod


def buildConsensus(trainData, learners, MLMethods, logFile = None):
        log(logFile, "Building a consensus model based on optimized MLmethods: "+str([ml for ml in MLMethods])+"...")
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
        log(logFile, "  Training Consensus Learner")
        return learner(trainData)

                    
def buildModel(trainData, MLMethod, queueType = "NoSGE", verbose = 0, logFile = None):
        """
        Buld the method passed in MLMethod and optimize ( "IndividualStatistics"  not in MLMethod)
        if MLMethod is a Consensus ("individualStatistics"  in MLMethod) , build each and optimize first all models and after build the consensus!
        """
        log(logFile, "Building and optimizing learner: "+MLMethod["MLMethod"]+"...")
        learners = {}
        MLMethods = {}
        if "IndividualStatistics"  in MLMethod:                        #It is a consensus
            for ML in MLMethod["IndividualStatistics"]:
                MLMethods[ML] = MLMethod["IndividualStatistics"][ML]
        else:
            MLMethods[MLMethod["MLMethod"]] = MLMethod

        # optimize all MLMethods
        for ML in MLMethods:
            log(logFile, "  Optimizing MLmethod: "+ML)
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
            log(logFile, "  Building the optimized learner:"+learners.keys()[0])
            model = learners[learners.keys()[0]](trainData)
        elif len(learners) >= 1:
            model = buildConsensus(trainData,learners,MLMethods)
        else:
            print "ERROR: No Learners were selected!"
            return None

        return model


def getModel(trainData, savePath = None, queueType = "NoSGE", verbose = 0):
        """
            Chooses the best model based on calculated MLStatistics
            trainData           Data for calculating the MLStatistics and finaly for training the selected model
            savePath            A path to a file for saving the MLStatistics (using pickle)
                                A log will be created in same place as savePath but with extension .log
            queueType              'NoSGE'   (without access to the distributed environment) - default
                                   'batch.q'
                                   'quick.q' (jobs start immediatly but are terminated after 30 min)
            verbose             Define a verbose level (default = 0)
 
        """
        if savePath:
            logFile = os.path.splitext(savePath)[0]+".log"
            file = open(logFile,"w")
            file.write("== AutoQSAR pipelie: getModel ==\n\n")
            file.close()
        else:
            logFile = None
        MLStatistics = getMLStatistics(trainData, savePath, queueType = queueType, verbose = verbose, logFile = logFile)
        MLMethod = selectModel(MLStatistics, logFile = logFile)
        model = buildModel(trainData, MLMethod, queueType = queueType, verbose = verbose, logFile = logFile)
        log(logFile, "-"*20)
        log(logFile, "getModel is returning: "+str(model)+"\n\n")
        return model


def getStatistics(dataset):
        pass

if __name__ == "__main__":
        data = orange.ExampleTable("./dataReg.tab")
        model = getModel(data, savePath = "./MLStat_reg.txt")
        print model

        data = orange.ExampleTable("./dataClass.tab")
        model = getModel(data, savePath = "./MLStat_class.txt")
        print model

