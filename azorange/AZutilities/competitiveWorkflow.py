import commands, os,time,pickle,statc, sys,copy
import AZOrangeConfig as AZOC
import AZLearnersParamsConfig as OPTconf
import orange,orngTest
import getUnbiasedAccuracy
from AZutilities import miscUtilities
from AZutilities import evalUtilities
from AZutilities import dataUtilities
from AZutilities import paramOptUtilities
from trainingMethods import AZBaseClasses
from trainingMethods import AZorngConsensus
from pprint import pprint
import random
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

def saveMLStatistics(savePath, MLStatistics, logFile=None):
         if savePath and os.path.isdir(os.path.split(savePath)[0]):
            file = open(savePath, "w")
            pickle.dump(MLStatistics, file)
            file.close()
            log(logFile, "MLStatistics saved to: "+savePath)


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
        evaluator = getUnbiasedAccuracy.UnbiasedAccuracyGetter(data = trainData, learner = learners, paramList = None, nExtFolds = AZOC.QSARNINNERFOLDS, nInnerFolds = AZOC.QSARNCVFOLDS, queueType = queueType, verbose = verbose, logFile = logFile, resultsFile = savePath)
        MLStatistics = evaluator.getAcc()

        saveMLStatistics(savePath, MLStatistics, logFile)
        return MLStatistics



def selectModel(MLStatistics, logFile = None):
        """Return the model with highest Q2/CA amongst methods with a stability less than 0.1.
           If no methods is considered stable, select the method with the greatest Q2/CA
        """
        log(logFile, "Selecting MLmethod...")
        bestModelName = None
        bestRes = None
        bestStableVal = None
        #Select only from stable models
        for modelName in MLStatistics:
            StabilityValue = MLStatistics[modelName]["StabilityValue"]
            if StabilityValue is not None:
                if MLStatistics[modelName]["responseType"] == "Classification": 
                    if statc.mean(MLStatistics[modelName]["foldStat"]["nTestCmpds"]) > 50:
                        stableTH = AZOC.QSARSTABILITYTHRESHOLD_CLASS_L
                    else:
                        stableTH = AZOC.QSARSTABILITYTHRESHOLD_CLASS_H
                elif MLStatistics[modelName]["responseType"] == "Regression":
                    if statc.mean(MLStatistics[modelName]["foldStat"]["nTestCmpds"]) > 50:
                       stableTH = AZOC.QSARSTABILITYTHRESHOLD_REG_L
                    else:
                        stableTH = AZOC.QSARSTABILITYTHRESHOLD_REG_H
                if StabilityValue < stableTH:
                    valRes = max( MLStatistics[modelName]["Q2"], MLStatistics[modelName]["CA"])  # One of them is always None
                    if bestRes is None or valRes > bestRes:
                        bestRes = valRes
                        bestModelName = modelName
                        bestStableVal = StabilityValue
                    elif valRes == bestRes and StabilityValue < bestStableVal:
                        bestRes = valRes
                        bestModelName = modelName
                        bestStableVal = StabilityValue
                    

        # No stable models found! Selecting the one with best result still... 
        if bestModelName is None:
            log(logFile, "  No stable models found! Selecting the one with best result still...")
            for modelName in MLStatistics:
                valRes = max( MLStatistics[modelName]["Q2"], MLStatistics[modelName]["CA"])  # One of them is always None
                if bestRes is None or valRes > bestRes:
                    bestRes = valRes
                    bestModelName = modelName
            log(logFile, "  Selected the non-stable MLmethod: " + bestModelName)
        else:
            log(logFile, "  Selected the stable MLmethod: " + bestModelName)
        MLMethod = copy.deepcopy(MLStatistics[bestModelName])
        MLMethod["MLMethod"] = bestModelName
        MLStatistics[bestModelName]["selected"] = True
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
                exprTest0 += "+( "+ml+" == "+CLASS0+" )*"+str(MLMethods[ml]["optAcc"])+" "
            exprTest0 += ")/IF0(sum([False"
            for ml in MLMethods:
                exprTest0 += ", "+ml+" == "+CLASS0+" "
            exprTest0 += "]),1)"
            exprTest1 = exprTest0.replace(CLASS0,CLASS1)
            expression = [exprTest0+" >= "+exprTest1+" -> "+CLASS0," -> "+CLASS1]
        else:
            Q2sum = sum([MLMethods[ml]["optAcc"] for ml in MLMethods])
            expression = "(1 / "+str(Q2sum)+") * (0"
            for ml in MLMethods:
                expression += " + "+str(MLMethods[ml]["optAcc"])+" * " + ml +" "
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
                MLMethods[ML] = copy.deepcopy(MLMethod["IndividualStatistics"][ML])
        else:
            MLMethods[MLMethod["MLMethod"]] = MLMethod

        # optimize all MLMethods
        for ML in MLMethods:
            log(logFile, "  Optimizing MLmethod: "+ML)
            learners[ML] = MLMETHODS[ML](name = ML)

            runPath = miscUtilities.createScratchDir(baseDir = AZOC.NFS_SCRATCHDIR, desc = "competitiveWorkflow")
            trainData.save(os.path.join(runPath,"trainData.tab"))

            tunedPars = paramOptUtilities.getOptParam(
                learner = learners[ML],
                trainDataFile = os.path.join(runPath,"trainData.tab"),
                useGrid = False,
                verbose = verbose,
                queueType = queueType,
                runPath = runPath,
                nExtFolds = None,
                logFile = logFile,
                getTunedPars = True)

            
            if not learners[ML].optimized:
                print "WARNING: competitiveWorkflow: The learner "+str(learners[ML])+" was not optimized."
                #print "         Using default parameters"
                print "         The "+str(learners[ML])+" will not be included"
                #print "         Returning None"
                print "             DEBUG can be made in: "+runPath 
                #Setting default parameters
                #learners[ML] = learners[ML].__class__()   
                #return None
                learners.pop(ML)
                continue
            else:
                print "Optimized learner ",learners[ML]      
                if trainData.domain.classVar.varType == orange.VarTypes.Discrete:
                    MLMethods[ML]["optAcc"] = tunedPars[0] 
                else:
                    res = orngTest.crossValidation([learners[ML]], trainData, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible, randomGenerator = random.randint(0, 100))
                    R2 = evalUtilities.R2(res)[0]  
                    MLMethods[ML]["optAcc"] = R2
                miscUtilities.removeDir(runPath)
        #Train the model
        if len(learners) == 1:
            log(logFile, "  Building the learner:"+learners.keys()[0])
            model = learners[learners.keys()[0]](trainData)
        elif len(learners) >= 1:
            model = buildConsensus(trainData,learners,MLMethods)
        else:
            print "ERROR: No Learners were selected!"
            return None

        return model


def getModel(trainData, savePath = None, queueType = "NoSGE", verbose = 0, getAllModels = False):
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
            file.write("== competitiveWorkflow pipelie: getModel ==\n\n")
            file.close()
        else:
            logFile = None
        MLStatistics = getMLStatistics(trainData, savePath, queueType = queueType, verbose = verbose, logFile = logFile)
        MLMethod = selectModel(MLStatistics, logFile = logFile)
        #Save again the MLStatistics to update the selected flag
        saveMLStatistics(savePath, MLStatistics, logFile) 
        if getAllModels:
            models = {}
            for ml in MLStatistics:
                MLStatistics[ml]["MLMethod"] = ml
                models[ml] = buildModel(trainData, MLStatistics[ml], queueType = queueType, verbose = verbose, logFile = logFile)
            log(logFile, "-"*20)
            log(logFile, "getModel is returning all models: "+str(models)+"\n\n")
            return models
        else:
            model = buildModel(trainData, MLMethod, queueType = queueType, verbose = verbose, logFile = logFile)
            log(logFile, "-"*20)
            log(logFile, "getModel is returning the selected model: "+str(model)+"\n\n")
            return {MLMethod["MLMethod"]:model}


def createStatObj(results=None, exp_pred=None, nTrainCmpds=None, nTestCmpds=None, responseType=None, nExtFolds=None, userAlert = ""):
        #Initialize res (statObj) for statistic results
        res = {}
        # Classification
        res["CA"] = None
        res["CM"] = None
        res["MCC"] = None
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
        if not results or results is None or exp_pred is None or responseType is None or nExtFolds is None or nTestCmpds is None or nTrainCmpds is None:
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



def getStatistics(dataset, runningDir, resultsFile, queueType = "batch.q", verbose = 0, getAllModels = False):
        """
                runningDir           (An existing dir for creating one job dir per fold)
                    |
                    +---- status     (The overall status:   "started", "finished" or the progress "1/10", "2/10", ...)
                    |
                    +---- fold_1
                    |
                    +---- fold_2
                    |
                    .
                    .
                    .
               
            The running will be monitorixed by this method.
            Whenever a MLMethod fails the respective fold job is restarted 
        """
        if dataset.domain.classVar.varType == orange.VarTypes.Discrete: 
            responseType = "Classification"
        else:
            responseType = "Regression"
        #Create the Train and test sets
        DataIdxs = dataUtilities.SeedDataSampler(dataset, AZOC.QSARNEXTFOLDS )
        #Check data in advance so that, by chance, it will not faill at the last fold!
        #for foldN in range(AZOC.QSARNEXTFOLDS):
            #trainData = dataset.select(DataIdxs[foldN],negate=1)
            #checkTrainData(trainData)

        jobs = {}
        thisDir = os.getcwd()
        os.chdir(runningDir)
        #PID = os.getpid() 
        #print "Started getStatistics in Process with PID: "+str(PID)
        #os.system('echo "'+str(PID)+'" > '+os.path.join(runningDir,"PID"))
        os.system('echo "started" > '+os.path.join(runningDir,"status"))
        # Start  all Fold jobs
        for fold in range(AZOC.QSARNEXTFOLDS):
            job = str(fold)
            print "Starting job for fold ",job
            trainData = dataset.select(DataIdxs[fold],negate=1)
            jobs[job] = {"job":job,"path":os.path.join(runningDir, "fold_"+job), "running":False, "failed":False, "finished":False}

            # Uncomment next 3 lines for running in finished jobs dirs
            #st, jID = commands.getstatusoutput("cat "+os.path.join(runningDir, "fold_"+job,"jID"))
            #jobs[job]["jID"] = jID
            #continue




            os.system("rm -rf "+jobs[job]["path"])
            os.system("mkdir -p "+jobs[job]["path"])
            trainData.save(os.path.join(jobs[job]["path"],"trainData.tab"))
            file_h = open(os.path.join(jobs[job]["path"],"run.sh"),"w")
            file_h.write("#!/bin/tcsh\n")
            file_h.write("source /home/palmeida/dev/AZOrange/templateProfile\n")
            file_h.write("python "+os.path.join(jobs[job]["path"],"QsubScript.py")+"\n")
            file_h.close()

            file_h = open(os.path.join(jobs[job]["path"],"QsubScript.py"),"w")
            file_h.write("import os\n")
            file_h.write("from AZutilities import dataUtilities\n")
            file_h.write("from AZutilities import competitiveWorkflow\n")
            file_h.write("data = dataUtilities.DataTable('"+os.path.join(jobs[job]["path"],"trainData.tab")+"')\n")
            file_h.write('os.system(\'echo "running" > '+os.path.join(jobs[job]["path"],"status")+' \')\n')
            file_h.write("models = competitiveWorkflow.getModel(data, savePath = '"+os.path.join(jobs[job]["path"],"results.txt")+"', queueType = '"+queueType+"', getAllModels = "+str(getAllModels)+")\n")
            file_h.write("nModelsSaved = 0\n")
            file_h.write("for model in models:\n")
            file_h.write("    if not models[model] is None:\n")
            file_h.write("        models[model].write('"+os.path.join(jobs[job]["path"],"model")+"'+'_'+model)\n")
            file_h.write('        nModelsSaved += 1\n')
            file_h.write('if nModelsSaved == len([m for m in models if not models[m] is None ]):\n')
            file_h.write('    os.system(\'echo "finished" > '+os.path.join(jobs[job]["path"],"status")+' \')\n')
            file_h.write('else:\n')
            file_h.write('    os.system(\'echo "failed" > '+os.path.join(jobs[job]["path"],"status")+' \')\n')
            file_h.close()

            os.chdir(os.path.join(jobs[job]["path"]))
            status, out = commands.getstatusoutput("qsub -cwd -q batch.q " + os.path.join(jobs[job]["path"],"run.sh"))
            if status:
                print "ERROR on Job "+str(job)+" (will be skipped)"
                print out
                #raise Exception("ERROR starting job for folder "+str(job))
            # Your job 955801 ("template_run.sh") has been submitted
            jID = out.strip().split(" ")[2]
            print "    jID: ",jID
            os.system('echo "'+jID+'" > '+os.path.join(jobs[job]["path"], "jID"))
            jobs[job]["running"] = True
            jobs[job]["jID"] = jID
            os.chdir(runningDir)
        os.chdir(thisDir)

        #Monitor Fold jobs
        finished = []
        updateJobsStatus(jobs)
        for job in jobs:
            if jobs[job]["finished"]:
                finished.append(job)
        print "Jobs already finished: ",finished
        os.system(' echo "'+str(len(finished))+'/'+str(AZOC.QSARNEXTFOLDS)+'" > '+os.path.join(runningDir,"status"))
        while len(finished) < AZOC.QSARNEXTFOLDS:
            print ".",
            sys.stdout.flush() 
            updateJobsStatus(jobs)
            for job in jobs:
                if jobs[job]["finished"] and job not in finished:
                    finished.append(job)
                    print time.asctime()+": Finished job "+str(job)
            os.system(' echo "'+str(len(finished))+'/'+str(AZOC.QSARNEXTFOLDS)+'" > '+os.path.join(runningDir,"status"))
            for job in [j for j in jobs if jobs[j]["failed"]]:
                jobs[job] = restartJob(jobs[job]) 
            time.sleep(5)                

        print "All fold jobs finished!"
        # Gather the results
        print "Gathering results..."
        #Var for saving each Fols result
        results = {}
        exp_pred = {}
        nTrainEx = {}
        nTestEx = {}
        # Var for saving the statistics results
        statistics = {}

        mlMethods = [ml for ml in AZOC.MLMETHODS] + ["Consensus"] 
        sortedJobs = [job for job in jobs]
        sortedJobs.sort(cmp = lambda x,y:int(x)>int(y) and 1 or -1)
        # Place for storing the selected models results
        results["selectedML"] = []
        exp_pred["selectedML"] = []
        nTrainEx["selectedML"] = []
        nTestEx["selectedML"] = []

        for ml in mlMethods:   # Loop over each MLMethod
            try:
                #Var for saving each Fols result
                results[ml] = []
                exp_pred[ml] = []
                nTrainEx[ml] = []
                nTestEx[ml] = []
                logTxt = ""

                
                for job in sortedJobs:   #loop over each fold
                    modelPath = os.path.join(jobs[job]["path"], "model_"+ml)
                    if not os.path.isdir(modelPath):
                        print "MLMethod "+ml+" not available in fold "+job
                        continue

                    resFile = os.path.join(jobs[job]["path"], "results.txt")
                    statFile_h = open(resFile)
                    foldStat = pickle.load(statFile_h)
                    statFile_h.close()

                    #load model
                    model = AZBaseClasses.modelRead(modelPath)
                    #Test the model
                    testData = dataset.select(DataIdxs[int(job)])
                    nTrainEx[ml].append(model.NTrainEx)
                    nTestEx[ml].append(len(testData))
                    if foldStat[ml]["selected"]:
                        nTrainEx["selectedML"].append(model.NTrainEx)
                        nTestEx["selectedML"].append(len(testData))

                    if responseType == "Classification":
                        results[ml].append((evalUtilities.getClassificationAccuracy(testData, model), evalUtilities.getConfMat(testData, model) ) )
                        if foldStat[ml]["selected"]:
                            results["selectedML"].append(results[ml][-1])
                    else:
                        local_exp_pred = []
                        for ex in testData:
                            local_exp_pred.append((ex.getclass(), model(ex)))
                        results[ml].append((evalUtilities.calcRMSE(local_exp_pred), evalUtilities.calcRsqrt(local_exp_pred) ) )
                        #Save the experimental value and correspondent predicted value
                        exp_pred[ml] += local_exp_pred
                        if foldStat[ml]["selected"]:
                            results["selectedML"].append(results[ml][-1])
                            exp_pred["selectedML"]+= local_exp_pred
                res = createStatObj(results[ml], exp_pred[ml], nTrainEx[ml], nTestEx[ml],responseType, len(sortedJobs), logTxt)
                if not res:
                    raise Exception("No results available!")
                statistics[ml] = copy.deepcopy(res)
                writeResults(statistics, resultsFile)
                print "       OK",ml
            except:
                print "Error on MLmethod "+ml+". It will be skipped"
        ml = "selectedML"
        res = createStatObj(results[ml], exp_pred[ml], nTrainEx[ml], nTestEx[ml],responseType, len(sortedJobs), logTxt)
        if not res:
            raise Exception("No results available!")
        statistics[ml] = copy.deepcopy(res)
        writeResults(statistics, resultsFile)
        os.system(' echo "finished" > '+os.path.join(runningDir,"status"))
        return statistics

def writeResults(statObj, resultsFile):
        if resultsFile and os.path.isdir(os.path.split(resultsFile)[0]):
            file = open(resultsFile, "w")
            pickle.dump(statObj, file)
            file.close()

def updateJobsStatus(jobs):
        # read the jobs qsub  status
        #  job-ID     prior      name       user    state        submit/start                     queue         slots ja-task-ID 
        #['955826', '0.50586', 'run.sh', 'palmeida', 'r', '07/14/2011', '19:18:17', 'batch.q@semldxwolf.seml.astraz', '1']
        status, out = commands.getstatusoutput("qstat")
        qstat = {}
        for job in out.split('\n')[2:]:
            status = job.split()
            qstat[status[0]] = status[4]

        for job in [j for j in jobs if not jobs[j]["finished"]]:
            statusFile = os.path.join(jobs[job]["path"],"status")
            if os.path.isfile(statusFile):
                st, status = commands.getstatusoutput("cat "+statusFile)
            else:
                print "WARNING: Missing status file"
                status = None
            if not status:
                print "WARNING! job "+job+" has no status! Will check next time"
                continue 
                

            #print "Fold ",job," Status: ",status," ID:",jobs[job]["jID"]," Qstat",qstat
            if jobs[job]["jID"] in qstat:
                if 'E' not in qstat[jobs[job]["jID"]]:
                    jobs[job]["running"] = True
                else:
                    jobs[job]["running"] = False
            else:
                jobs[job]["running"] = False

            if not isJobProgressingOK(jobs[job]):
                print "Job "+job+" failed to build one or more models in getMLStatistics. It was flaged for restart"
                jobs[job]["failed"] = True
                return
            if not jobs[job]["running"]:
                    #Test if it finished properly
                    if status == "failed":
                        print "Job "+job+" failed to build all models. It was flaged for restart"
                        jobs[job]["failed"] = True
                    elif status == "finished":
                        jobs[job]["finished"] = True
                    else:
                        wait = 3*60
                        print "Job "+job+" seems to be finished ("+status+"). Waiting "+str(wait)+" sec. in case the flag is being updated"
                        time.sleep(wait)
                        st, status = commands.getstatusoutput("cat "+statusFile) 
                        if status != "finished":
                            print "Job "+job+" failed to build one or more  models. It was flaged for restart"
                            jobs[job]["failed"] = True



def restartJob(jobObj, force = False):
            job = jobObj["job"]
            print "\nJob "+job+ " is being reported as failing"
            runningJobDir = jobObj["path"]
            thisDir = os.getcwd()
            os.chdir(runningJobDir)
            if not os.path.isfile("./jID"):
                print "No job started yet for fold "+job
                
            else:       
                print "Killing eventually running job"+job+"..."
                os.system("cat jID | xargs qdel ")
                status, oldjID = commands.getstatusoutput("head -n 1 jID")
                # Save results
                print "  Backing up Job "+str(job)+"..."
                os.system("mkdir Bkup_"+oldjID)
                os.system("mv model_* Bkup_"+oldjID)
                os.system("mv jID Bkup_"+oldjID)
                os.system("mv results.* Bkup_"+oldjID)
                os.system("mv run.sh.* Bkup_"+oldjID)
                os.system("mv status Bkup_"+oldjID)
            print "  Starting Job "+str(job)+"..."
            jobFile = os.path.join(runningJobDir,"run.sh")
            status, out = commands.getstatusoutput("qsub -cwd -q batch.q " + jobFile)
            if status:
                print "  ERROR on Job "+str(job)+" (will be skipped)"
                print out
                os.chdir(thisDir)
                return
            # Your job 955801 ("template_run.sh") has been submitted
            jID = out.strip().split(" ")[2]
            print "    jID: ",jID
            os.system('echo "'+jID+'" > '+os.path.join(runningJobDir,"jID"))
            os.system('echo "Restarted at '+str(time.asctime())+'" >> '+os.path.join(runningJobDir,"restarts.log"))
            jobObj["running"] = True
            jobObj["finished"] = False
            jobObj["failed"] = False
            jobObj["jID"] = jID
            os.chdir(thisDir)
            return jobObj


def isJobProgressingOK(job):
        runningJobDir = job["path"] 
        resFile = os.path.join(runningJobDir, "results.txt")
        if not os.path.isfile(resFile):
            return True
        n = 50
        statistics = None
        while n > 0:
            try:
                statFile_h = open(resFile )
                statistics = pickle.load(statFile_h)
                statFile_h.close()
                n=0
            except:
                statistics = None
                print "   No results file in job "+job["job"]+". Waiting 5 sec and trying again..."
                n -=1
                time.sleep(5)
        if not statistics:
            print "   Results file is missing in job "+job["job"]+". will check next time!"
            return True
        mlFailed = []
        for ml in statistics:
            if ml == "Consensus":
                continue
            if statistics[ml]["CA"] is None and statistics[ml]["Q2"] is None:
                mlFailed.append(ml)
        if mlFailed:
            for ml in mlFailed:
                print "MLMethod "+ml+" of job "+job["job"]+" failed"
            return False
        else:
            return True






if __name__ == "__main__":
        data = orange.ExampleTable("./dataReg.tab")
        
        #statistics = getStatistics(data, runningDir="/home/palmeida/dev/AZOrange/azorange/AZutilities/runningJobs")
        #print "Statistics = ",statistics
        #sys.exit()

        model = getModel(data, savePath = "./MLStat_reg.txt")
        print model

        data = orange.ExampleTable("./dataClass.tab")
        model = getModel(data, savePath = "./MLStat_class.txt")
        print model

