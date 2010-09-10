import statc
import random
import time
import commands
import orange
import types
import string
import os, sys
import orngTest
import AZLearnersParamsConfig
from AZutilities import miscUtilities
from AZutilities import dataUtilities
from AZutilities import evalUtilities
import AZOrangeConfig as AZOC
from math import sqrt
from math import pow
from math import floor
import threading
import thread
from copy import deepcopy
import traceback   
from glob import glob
import sgeUtilities
 
version = 11

class Appspack:
    userVars = ("qsubFile","advancedMPIoptions","np","machinefile","externalControl","useParameters", "learner", "dataSet", "runPath", "verbose",\
                "evaluateMethod", "findMin", "samplingMethod", "nFolds","useGridSearchFirst","gridSearchInnerPoints", "queueType", "nExtFolds", \
                "useStd") 
    LEAVE_ONE_OUT         = 0
    FOLD_CROSS_VALIDATION = 1
    def __init__(self, **kwds):
        #Possible user defined Vars
        self.externalControl = 0  #if = 1, the control of the finish of optimization must be from the caller by the 
                                  # isFinished() method
        self.useParameters = None
        self.useDefaultPoint = True
        self.machinefile = None         # path to the file containing the machines nodes to use OR a list of machines
                                        #      ex: ["localhost","localhost"]
        self.np = None                  # number of processors to use when using MPI
        self.usedMPI = False            # Flag indicating if last optimization done was made using MPI version of appspack
        self.advancedMPIoptions = ""    # Other MPI valid options that user might want to set.
        self.learner = None             # The Learner to be optimized
        self.dataSet = ""               # the PATH for the dataset to use for optimization
        self.runPath = "./"             # the runPath for the optimization
        self.verbose = 0                # Verbose Flag
        self.evaluateMethod = "AZutilities.evalUtilities.RMSE"    # The evaluation methos to use
        self.findMin = True             # If the max is to be found, the caller must define findMin=False in the **kwds when calling
        self.samplingMethod = self.FOLD_CROSS_VALIDATION         # Sampling Method for the evaluation method
        self.nFolds = 5                 # the number of folds if FOLD_CROSS_VALIDATION is to be used
        self.gridSearchInnerPoints   = 5     # Number of innerpoint to split each variable interval
        self.useGridSearchFirst = False  # Flag indicating if GridSearch is to be performed to define the initial point of appspack search. If false, the midrange point will be used
        self.nExtFolds = None            # The number of folds to use in a loop over CV with different seeds. To reduce the 
                                         # influence of data sampling on the generalization accuracy of each model parameter point.
        self.useStd = True               # Do not select optimized parameter unless the accuracy differenc is significant.
        # Append arguments to the __dict__ member variable 
        self.__dict__.update(kwds)


        # Non-User defined vars
        self.parameters = None          # All the optimization parameters defined by the user in useParameters. 
                                        #     If not defined, they will be the same as the origParameters
        self.origParameters = None      # All the optimization parameters present in the static AZLearnersParamsConfig.py file
        self.finishedFlag = True        # Flag indicating the termination of Optimization
        self.appspackPID = 0            # the PID of current appspack process for this object instance
        self.tunedParameters = "No parameters yet."                     # the tuned parameters after success optimization
        self.qsubFile = None            # Name of file for sge job
        self.qsubJobId = None           # qsub job id
        self.GSRes     = None           # Results of GridSearch
        self.queueType = "batch.q"      # batch.q or quick.q - controlls which type of sge queue that will run the job
        self.STDevalRes = 0             # std in evalRes originating from data sampling effects. The learner parameters are 
                                        # only changed if the improvement in accuracy is greater than STDevalRes
        self.nStdFolds = 10              # Number of folds used to assess the std. Increase when parallel
        #self.defaultPoint = None


    def isFinished(self):
        """
        Returns True if the optimization process is finished,
                False if it still running
        """
        if self.machinefile == "qsub":
            if self.qsubJobId:
                # There is a jobid and qsub is not running -> finished
                if not self.getIsQsubRunning(): 
                    self.finishedFlag = True
                    # If the job finished put back the original ssh directory
                    #self.deleteTempSSHkey()
                    self.assignTunedParameters()
                else: self.finishedFlag = False
            ##scPA if the qsubJobId does not exist, then the process did never started in first place, then report as finished
            #else: self.finishedFlag = False
            else: self.finishedFlag = True
            ##ecPA
            return self.finishedFlag
        else:
            if self.appspackPID == 0:
                return self.finishedFlag
            else:
                readsDone = 0
                while readsDone <5:
                    try:                 
                        psLine = os.popen( "ps -p "+str(self.appspackPID) ).read()
                        readsDone = 5
                    except:
                        time.sleep(0.5)
                    readsDone += 1
                if ("appspack" in psLine or "mpirun" in psLine or "qsub" in psLine) and "<defunct>" not in psLine:
                    return False
                else:
                    self.assignTunedParameters()
                    self.finishedFlag = True
                    self.appspackPID = 0
                    return True 

    def stop(self):
        if  self.isFinished():
                self.appspackPID = 0
                self.finishedFlag = True
        else:
                if self.machinefile == "qsub":
                    if self.getIsQsubRunning() and self.qsubJobId:
                        os.system("qdel "+self.qsubJobId)
                        #self.deleteTempSSHkey()
                    else:
                        print "Qsub job is not running"
                        return False
                else:
                    os.kill(self.appspackPID,9)
                    readsDone = 0
                    while readsDone <5:
                        try:
                            psLine = os.popen( "ps -p "+str(self.appspackPID) ).read()
                            readsDone = 5
                        except:
                            time.sleep(0.5)
                        readsDone += 1
                    if "appspack" not in psLine or "<defunct>" in psLine:
                            self.tunedParameters = "The optimizer was stopped by the user"
                            self.appspackPID = 0
                            self.finishedFlag = True
                    else:
                            return False
        return True

    def __call__(self,  **kwds):
        """
        Input parameters:
                learner = <learner with GetTunableParameters function>
                dataSet = <orange Examples table with class> OBS, for qsub jobs to work this has to the the full path!!
                [ runPath = <string with path for running the APPSPACK> OBS, for qsub jobs to work this has to the the full path!!]
                [ machinefile= <list of strings with the machines to use as nodes> | <full path of file compatible with mpich containing the machines list>, or string equal to qsub to execute the mpi job on sge]
                [ np= <number of processors to use in parallel>]
        Output: 
                The learner specified at input is set up with optimal parameters
                and the "optimized" attribute is set to true or false if it was not possible to optimize parameters
                If there was an error, it will be returned a string with the occurrence
        """
        
        # Append arguments to the __dict__ member variable 
        self.__dict__.update(kwds)
        for kwd in kwds:
            if kwd in self.userVars:
                setattr(self,kwd,kwds[kwd])
            else:
                if self.verbose > 0: print "Variable "+kwd+" is not to be defined by user!" 
        #Globals: self.learner and self.RunPath
        if self.learner == None:
            return "No learner specified"
        if self.verbose > 0: print self.learner
        if self.finishedFlag == False or self.appspackPID != 0:
            return "appspack already running in background for this instance"
        # attribute "optimized" set to false until sucessfull optimization
        if hasattr(self.learner, "setattr"):
            self.learner.setattr("optimized", False)
        else:
            setattr(self.learner,"optimized", False)

        # check if learner is defined in AZLearnersParamsConfig
        self.learnerType = str(self.learner).split()[0]
        if self.learnerType[0] == "<":
            self.learnerType = self.learnerType[1:]
        if self.learnerType.rfind('.') >=0:
            self.learnerType = self.learnerType[self.learnerType.rfind('.')+1:]

        if not self.useParameters:
            if not hasattr(AZLearnersParamsConfig, self.learnerType):
                return "No configuration for specified learner " + self.learnerType +" in AZLearnersParamsConfig"
            self.parameters = eval("AZLearnersParamsConfig." + self.learnerType)
            self.origParameters = deepcopy(self.parameters)
        else:
            self.parameters = self.useParameters
            if hasattr(AZLearnersParamsConfig, self.learnerType):
                self.origParameters = eval("AZLearnersParamsConfig." + self.learnerType)
            else:
                self.origParameters = deepcopy(self.parameters)


        # Validate input parameters
        if not os.path.exists(self.dataSet):
            return "Data set not found!"
        #clean some old data
        self.tunedParameters = "No parameters yet!"        
        #Get the number of Attributies in DataSet
        dataInfo = dataUtilities.getQuickDataSize(self.dataSet)
        #Check if the evaluation method is correct according to the dataset
        try:
            exec("import "+self.evaluateMethod[:self.evaluateMethod.rfind(".")]) 
            evalType = eval(self.evaluateMethod)()
        except:
            evalType = None
            print "WARNING: The method "+self.evaluateMethod+"does not seem valid!"
        if dataInfo["discreteClass"] == -1:
            print "WARNING: It was not possible to determine the type of the class variable. It is possible that the evaluation method is not the correc for this specific problem"
        elif type(evalType) != types.DictType or "type" not in evalType:
            print "WARNING: The type of the evaluation method specified ("+self.evaluateMethod+") cannot be determined."
        elif dataInfo["discreteClass"] == 1 and ((evalUtilities.CLASSIFICATION  & evalType["type"]) != evalUtilities.CLASSIFICATION):
            print "ERROR: The evaluation method specified ("+self.evaluateMethod+") cannot be used with Classification Problems"
            return None
        elif dataInfo["discreteClass"] == 0 and ((evalUtilities.REGRESSION  & evalType["type"]) != evalUtilities.REGRESSION):
            print "ERROR: The evaluation method specified ("+self.evaluateMethod+") cannot be used with Regression Problems"
            return None

        self.nAttrInDataSet = dataInfo["N_ATTR"] #len(data.domain.attributes)
        if self.samplingMethod == self.LEAVE_ONE_OUT:        #sMethod will be orngTest.leaveOneOut
           self.nExamplesInDataSet = dataInfo["N_EX"]-1 #len(data)-1
           self.nFolds = 1
        else:    ## by default: samplingMethod = self.FOLD_CROSS_VALIDATION      # sMethod will be orngTest.crossValidation with folds = self.nFolds
           self.nExamplesInDataSet = dataInfo["N_EX"] - floor(dataInfo["N_EX"]/self.nFolds)
        if not os.path.isdir(self.runPath):
            if self.verbose > 0: print "Warning! Specified run path does not exist! ./ will be used!"
            self.runPath = "./"
        if self.runPath[-1] != "/":
            self.runPath += "/"
        paramFile=file(os.path.join(self.runPath,"AZLearnersParamsConfig.py"),"w")
        paramFile.write(self.learnerType + "= " + str(self.parameters)+"\r\n")
        paramFile.close()
 
 
        # delete the intermediate results output file if exists
        #if os.path.exists(self.runPath+"intRes.txt"):
        #    os.remove(self.runPath+"intRes.txt")
  
        for part in sorted(glob(os.path.join(self.runPath,"*intRes.txt"))):
            os.remove(part)

        # Create the input file for appspack (this includes calculating default point and midrange point)
        appsInput = self.__CreateInput()
        if appsInput == None:
            print "ERROR: Cannot create input file for appspack"
            return None
        # Run appspack
        retVal = self.__RunAppspack(appsInput)
        if retVal == None:
                if sorted(glob(os.path.join(self.runPath,"*intRes.txt"))) != []:
                #if os.path.isfile(os.path.join(self.runPath,"intRes.txt")):
                    self.assignTunedParameters()
                    return self.tunedParameters
                if self.verbose > 0: print "ERROR: __RunAppspack returned None!"
                return None
        if self.externalControl==0:     # the control is internal to this class, and then the appspack has finished
                self.assignTunedParameters()
                return self.tunedParameters
        else:                           # The control is external to this class, someone will check the finished state
                if self.machinefile == "qsub":
                    return int(self.qsubJobId)
                else:
                    return self.appspackPID
        
    def assignTunedParameters(self):

        if self.learner.optimized:
            tunedParameters = self.processAppspackResults()
        else:
            tunedParameters = self.processIntResResults()

        if tunedParameters == None:
            self.tunedParameters = "Could not optimize the parameters"
            return
        
        # Configure learner with parameters returned in tunedParameters
        #Types handle here
        for key in tunedParameters[1]:
            if hasattr(self.learner, "setattr"):
                if type(eval(self.parameters[key][0])) == types.ListType:
                    self.learner.setattr(key, [eval(self.parameters[key][0])[0](tunedParameters[1][key][1:-1])])
                elif  type(eval(self.parameters[key][0])) == types.TypeType and eval(self.parameters[key][0]) == types.TypeType:
                    self.learner.setattr(key, eval(tunedParameters[1][key]))
                else:
                    if eval(self.parameters[key][0]) == types.BooleanType:
                        self.learner.setattr(key, types.BooleanType(eval(tunedParameters[1][key])))
                    elif eval(self.parameters[key][0]) == types.IntType:
                        self.learner.setattr(key, int(float(tunedParameters[1][key])))
                    else:
                        self.learner.setattr(key, eval(self.parameters[key][0])(tunedParameters[1][key]))
            else:
                if type(eval(self.parameters[key][0])) == types.ListType:
                    setattr(self.learner, key, [eval(self.parameters[key][0])[0](tunedParameters[1][key][1:-1])])
                elif  type(eval(self.parameters[key][0])) == types.TypeType and eval(self.parameters[key][0]) == types.TypeType:
                    setattr(self.learner, key, eval(tunedParameters[1][key]))
                else:
                    if eval(self.parameters[key][0]) == types.BooleanType:
                        serattr(self.learner, key, types.BooleanType(eval(tunedParameters[1][key])))
                    elif eval(self.parameters[key][0]) == types.IntType:
                        setattr(self.learner, key, int(float(tunedParameters[1][key])))
                    else:
                        setattr(self.learner, key, eval(self.parameters[key][0])(tunedParameters[1][key]))
 
        if self.verbose > 0: print "tunedParameters = ", tunedParameters
        #optimization was successfull, so set "optimized" to true
        if hasattr(self.learner, "setattr"):
            self.learner.setattr("optimized", True)
        else:
            setattr(self.learner, "optimized", True)
        self.tunedParameters = tunedParameters

    def getTunedParameters(self):
        return {"bestRes":self.tunedParameters[0],"optParam":self.tunedParameters[1],"ResIdx":self.tunedParameters[2]}


    def getEvalRes(self):
        """
        Return the accuracy of an nFold CV.
        """

        if self.verbose > 0: print "The number of folds in the CV: "+str(self.nFolds)
       
        data = dataUtilities.DataTable(self.dataSet)
        res = orngTest.crossValidation([self.learner], data, folds=self.nFolds, strat=orange.MakeRandomIndices.StratifiedIfPossible, randomGenerator = random.randint(0, 100))
        evaluateMethod = self.evaluateMethod[self.evaluateMethod.find(".")+1:len(self.evaluateMethod)]
        evalRes = eval(evaluateMethod)(res)[0]

        return evalRes

    def getSTDevalRes(self):
        """
        Assess the standard deviation of the metrics used to evaluated the generalization error. 
        Recalcualte the generalization accuracy used in the parameter optimization 50 times 
        with different seeds in the data sampling.
        """

        jobScript = """\
import orange,orngTest,random,cPickle
from AZutilities import evalUtilities
from AZutilities import dataUtilities

paramFile=open("Params.pkl","r")
(learner,nFolds,dataSet,evaluateMethod) = cPickle.load(paramFile)
paramFile.close()
data = dataUtilities.DataTable(dataSet)
res = orngTest.crossValidation([learner], data, folds=nFolds, strat=orange.MakeRandomIndices.StratifiedIfPossible, randomGenerator = random.randint(0, 100))
evalMethod = evaluateMethod[evaluateMethod.find(".")+1:len(evaluateMethod)]
print cPickle.dumps(eval(evalMethod)(res)[0])
"""
 
        extEvalList = []
        evalResList = []

        # Set std folds depending on if we're running serial or parallel.
        if self.machinefile == "qsub":
           self.nStdFolds = 50  # parallel

        # set total folds depending if we use nExtFolds or not.
        if self.nExtFolds:
            folds = self.nExtFolds * self.nStdFolds
        else:
            folds = self.nStdFolds

        # Assess the memory requirements
        memSize = dataUtilities.getApproxMemReq(self.dataSet)

        # run on the sge as a parallel array job, or plain loop if serial
        if self.machinefile == "qsub":
            evalResList = sgeUtilities.arrayJob(jobName = "CvJob", jobNumber = folds, jobParams = [self.learner,self.nFolds,self.dataSet,self.evaluateMethod], jobQueue = "batch.q", jobScript = jobScript, memSize = str(memSize)+"M")  # parallel
        else:
            for idx in range(folds):
                evalRes = self.getEvalRes()
                evalResList.append(evalRes) #serial

        # if nExtFolds is used, compute the mean of the ExtFold.
        if self.nExtFolds:
            c = 0
            while c<folds:
                for idx in range(self.nExtFolds):
                    c=c+1
                evalRes = round(statc.mean(evalResList[c-idx-1:c]),3)
                extEvalList.append(evalRes) 
        else:
            extEvalList = evalResList # if no nExtFolds then just copy
        stdEvalRes= statc.std(extEvalList)

        return round(stdEvalRes,3)
            

    def processIntResResults(self):

        # Calculate the std in evalRes originating from data sampling effects. Needed regardless of nExtFolds.
        # Don't use the std until it runs efficiently on the SGE
        #self.STDevalRes = self.getSTDevalRes()
        self.STDevalRes = 0.000

        intResTxt=[]
        # Load all parts of the intRes.txt 
        for part in sorted(glob(os.path.join(self.runPath,"*intRes.txt"))):
                file = open(part,"r")
                intResTxt = intResTxt + file.readlines()
                file.close()

        #Retrieve Domain from the first line of intRes (This remover first line from intRes itself)
        resDomain = intResTxt.pop(0).split()
        #Convert the intResTxt into splited strings
        intRes = [line.split() for line in intResTxt]
        #ADD THE STDevalRes variable before start processing the intRes file
        resDomain.insert(-1,"STDevalRes")
        for idx,line in enumerate(intRes):
            intRes[idx].insert(-1,self.STDevalRes)

        #Create the optimizationLog.txt file for user usage
        intResFile = open(os.path.join(self.runPath,"optimizationLog.txt"),"w")
        userResDomain = [attr for attr in resDomain if ("LEARNER_" in attr or attr in ["EVAL_RES","STDevalRes"])]
        intResFile.write(string.join(userResDomain, sep='\t').replace("LEARNER_","")+"\n")
        if self.verbose > 0:
            print "Running Path:",self.runPath
            print "============== Intermediate results at intRes file ("+str(len(intRes))+" lines) =============="
        for line in intRes:
            if self.verbose > 0: print line
            # When running appspack_mpi on multiple nodes (> 4), crap is sometimes written to the last line. 
            #try: intResFile.write(string.join([line.split()[idx] for idx in [resDomain.index(i) for i in userResDomain]],"\t")+"\n")
            #except: pass
            intResFile.write(string.join([str(line[idx]) for idx in [resDomain.index(i) for i in userResDomain]],"\t")+"\n")
        if self.verbose > 0: print "\n"
        intResFile.close()

        # Find the best result!
        if self.findMin:
            bestIdx = 0
            for resIdx,res in enumerate(intRes):
                if float(res[resDomain.index("EVAL_RES")]) < float(intRes[bestIdx][resDomain.index("EVAL_RES")]):
                    bestIdx = resIdx
        else:
            bestIdx = 0
            for resIdx,res in enumerate(intRes):
                if float(res[resDomain.index("EVAL_RES")]) > float(intRes[bestIdx][resDomain.index("EVAL_RES")]):
                    bestIdx = resIdx

        if self.verbose > 0: print "The best result: ",
        evalResBest = float(intRes[bestIdx][resDomain.index("EVAL_RES")])
        if self.verbose > 0: print evalResBest
        if self.verbose > 0: print "The default result: ",
        evalResDefault = float(intRes[0][resDomain.index("EVAL_RES")])
        if self.verbose > 0: print evalResDefault
        if self.verbose > 0: print "STDevalRes: ",self.STDevalRes
   
        # Did the optimization result in a significantly better model?
        if self.useStd:
            if self.findMin:
                if evalResBest < (evalResDefault - self.STDevalRes): 
                    selectOptParam = True
                else:
                    selectOptParam = False
            else:  # Classifier
                if evalResBest > (evalResDefault + self.STDevalRes):
                    selectOptParam = True
                else:
                    selectOptParam = False
        else:
            selectOptParam = True

        if self.verbose > 0: print "The optimized parameters are selected: "+str(selectOptParam)

        # get the learner parameters used on that best point or the defautl parameters
        optParameters = {}
        for param in self.parameters:
             if selectOptParam: 
                 optParameters[param] = str(intRes[bestIdx][resDomain.index("LEARNER_"+param)])
             else: # Return default parameters
                 optParameters[param] = str(intRes[0][resDomain.index("LEARNER_"+param)])

        # Return the optimal values (check if the max of function was requested, if so, return the simetrical)
        # OBS select the optimal values only if they are more accurate than STDevalRes!!!
        if self.verbose > 0: print "Best Result from opimizer: ", intRes[bestIdx][resDomain.index("EVAL_RES")]
        if self.verbose > 0: print "Best Parameters: ",optParameters

        return [float(intRes[bestIdx][resDomain.index("EVAL_RES")]), optParameters,bestIdx]         #Ex: [2.3 ,  {"a":1, "b":"pl1"} ]


    def processAppspackResults(self):
        print "=================== WARNING ======================="
        print "Deprecated... The processAppspackResults function is about to be removed from the code."
        print "Please check why this function was called and report to Pedro Almeida!"
        return None
        # Read solution from output file self.runPath+"solution.txt"
        # Compares the result from appspack with the default point and choose the best one
        if not os.path.exists(self.runPath+"solution.txt"):
            if self.verbose > 0: print "ERROR (paramOptUtilities.py): not found solution of APPSPACK on ",self.runPath+"solution.txt"
            return None
        try:
            file = open(os.path.join(self.runPath,"solution.txt"),"r")
            solution = file.readline()
            file.close()

            inputApps = open(os.path.join(self.runPath,"input.apps"))
            paramKeys = eval(inputApps.readlines()[1][1:-2])
            inputApps.close()

            if solution.split()[1] == "]":
                if self.verbose > 0: print "ERROR (paramOptUtilities.py): Solution was not in the correct format. Probably OptScriptModel.py is the wrong version"
                return None
            f = types.FloatType(solution.split()[1])
            if not self.findMin:
                f = -f
            vars = [types.FloatType(x) for x in solution.split()[solution.split().index("x=[")+1:-1]]
        except:
            if self.verbose > 0: print "ERROR (paramOptUtilities.py): Unexpected error"
            return None
        if self.verbose > 0: print "Result from APPSPACK:"
        if self.verbose > 0: print "f (fixed signal) =    ", f
        if self.verbose > 0: print "vars = ",vars,"  ", paramKeys 
        # Transform APPSPACK variables to learner parameters
        if len(vars)!=len(paramKeys):
            if self.verbose > 0: print "ERROR (paramOptUtilities.py): APPSPACK returned different parameters that the ones defined"
            return None

        intRes=[]
        # Load all parts of the intRes.txt
        for part in sorted(glob(os.path.join(self.runPath,"*intRes.txt"))):
                file = open(part,"r")
                intRes = intRes + file.readlines()
                file.close()

        #Retrieve Domain from the first line of intRes (This remover first line from intRes itself)
        resDomain = intRes.pop(0).split()
        #Create the optimizationLog.txt file for user usage
        intResFile = open(os.path.join(self.runPath,"optimizationLog.txt"),"w")
        userResDomain = [attr for attr in resDomain if ("LEARNER_" in attr or attr=="EVAL_RES")]
        intResFile.write(string.join(userResDomain, sep='\t').replace("LEARNER_","")+"\n")
        for line in intRes:
            intResFile.write(string.join([line.split()[idx] for idx in [resDomain.index(i) for i in userResDomain]],"\t")+"\n")
        intResFile.close()
        
        #Find the intRes example that contains the best point found by appspack
        minDistIdx = 0
        minDist = 999999999
        for resIdx,res in enumerate(intRes):
            dist = pow(float(res.split()[resDomain.index("EVAL_RES")]) - f , 2)
            for idx,param in enumerate(paramKeys):
                dist += pow(float(res.split()[resDomain.index("APPS_" + param)]) - vars[idx] ,2)
            dist = sqrt(dist)
            if dist < minDist: 
                minDist = dist
                minDistIdx = resIdx
        # Check if the default point was not even better!
        if self.findMin:
            if float(intRes[0].split()[resDomain.index("EVAL_RES")]) <= float(intRes[minDistIdx].split()[resDomain.index("EVAL_RES")]):
                bestIdx = 0
            else:
                bestIdx = minDistIdx
            # double check if there is not any better result!
            for resIdx,res in enumerate(intRes):
                if float(res.split()[resDomain.index("EVAL_RES")]) < float(intRes[bestIdx].split()[resDomain.index("EVAL_RES")]):
                    bestIdx = resIdx
        else:
            if float(intRes[0].split()[resDomain.index("EVAL_RES")]) >= float(intRes[minDistIdx].split()[resDomain.index("EVAL_RES")]):
                bestIdx = 0
            else:
                bestIdx = minDistIdx
            # double check if there is not any better result!
            for resIdx,res in enumerate(intRes):
                if float(res.split()[resDomain.index("EVAL_RES")]) > float(intRes[bestIdx].split()[resDomain.index("EVAL_RES")]):
                    bestIdx = resIdx
        # get the learner parameters used on that best point 
        optParameters = {}
        for param in self.parameters:
             optParameters[param] = str(intRes[bestIdx].split()[resDomain.index("LEARNER_"+param)])


        # Return the optimal values (check if the max of function was requested, if so, return the simetrical)
        if self.verbose > 0: print "Best Result from opimizer: ", intRes[bestIdx].split()[resDomain.index("EVAL_RES")]
        if self.verbose > 0: print "Best Parameters: ",optParameters
        return [float(intRes[bestIdx].split()[resDomain.index("EVAL_RES")]), optParameters, bestIdx]         #Ex: [2.3 ,  {"a":1, "b":"pl1"} ]
 
    def __CreateInput(self):
        """Createds the input file for APPSPACK and returns its full path
             Creates the inputPath in self.runPath
                Uses the learner.GetTunableParameters() function and ignores the "StepInterval" attribute
             Returns None if it was not possible to create the input file
        """
        #self.parameters are parameters with syntax: 
        #   {str ParameterName: 
        #    [types parameterType, 
        #     str ValuesRangeType <"interval"|"values">,
        #     list ValuesRange <Values|IntervalLimits> ], ...} 
        if len(self.parameters) < 1:
            if self.verbose > 0: print "ERROR: Zero parameters specified"
            return None
        try:
            N_ATTR = self.nAttrInDataSet
            N_EX = self.nExamplesInDataSet
            upperVector = []
            lowerVector = []
            defaultX = []
            # Find in advance how many parameters are being asked for optimization
            nOptParams = 0
            for key in self.parameters:
                if self.parameters[key][5]:
                    nOptParams += 1
            # Convert the parameters list 
            paramKeys = []
            for key in self.parameters:
                if self.parameters[key][5]: # parameters with len of range =1 or optimization flag false are not to be optimized and not to be in appspack imput file 
                        valuesRange=eval(self.parameters[key][2])
                        if len(valuesRange) == 0:
                                if self.verbose > 0: print "ERROR: Empty range of values(",key,")"
                                if self.verbose > 0: print self.parameters[key][2]," = ",valuesRange
                                return None
                        if len(valuesRange) == 1 and  nOptParams > 1:
                                if self.verbose > 0: print "WARNING: Parameter ",key," will not be optimized and will be set to ",valuesRange[0]
                        elif self.parameters[key][1]=="values":
                                paramKeys.append(key)
                                upperVector.append(len(valuesRange)-1)
                                lowerVector.append(0)
                                #Get the default values for the vars to be optimized        
                                defaultX.append(self.parameters[key][4])
                        else:
                                paramKeys.append(key)
                                upperVector.append(valuesRange[1])
                                lowerVector.append(valuesRange[0])
                                #compute the initial point to be the defaults
                                defaultX.append(self.parameters[key][4])

            # Test the boundaries, we do not allow unlimit boundaries
            if "DNE" in upperVector or "DNE" in lowerVector:
                if self.verbose > 0: print "ERROR: DNE not permited in range vectors"
                return None

            #evaluate the function at the default point
            if self.useDefaultPoint:
                defaultF = self.__EvaluateDefaultPoint(defaultX,paramKeys)    #Use this one when using defaults point
            #self.defaultPoint = (defaultF, defaultX, paramKeys)

            if not self.useGridSearchFirst:
                    # Calculate mid-range point from lowerVector to upperVector as initial point (initialX)
                    initialX = len(upperVector)*[0]
                    for i in range(0,len(upperVector)):
                        initialX[i] = (upperVector[i]+lowerVector[i])/2
                    initialF = self.__EvaluateFunction(initialX, paramKeys)
            else:
                    # WARNING! Do not specify a dir inside self.runPath. This would create a circular reference!
                    #          This path must ba a path accessible by all nodes running each point
                    GSPath = miscUtilities.createScratchDir(desc ="GRidSearchFiles", baseDir = AZOC.NFS_SCRATCHDIR)
                    GS_script = self.__CreateGridSearchFiles(paramKeys)
                    from AZutilities import AZGridSearch
                    GridSearch = AZGridSearch.GridSeach(varsUpperLimits = upperVector,\
                                                        varsLowerLimits = lowerVector,\
                                                        nInnerPoints    = self.gridSearchInnerPoints,\
                                                        scriptFilesPath = self.runPath,\
                                                        scriptRunPath   = GSPath,\
                                                        ScriptFileName  = GS_script)
                    self.GSRes = GridSearch()
                    #Find the best initial point analyzing the results
                    if self.findMin:
                        bestIdx = self.GSRes["results"].index(min(self.GSRes["results"]))
                    else:
                        bestIdx = self.GSRes["results"].index(max(self.GSRes["results"]))
                    initialX = self.GSRes["varValues"][bestIdx]
                    initialF = self.GSRes["results"][bestIdx]
                    miscUtilities.removeDir(GSPath)
                    #print "GSPath ", GSPath 

            #Create the shell script that will set all the envVars and then call the runScript.py
            #This file will be called each time appspack want to calculate a point.
            #  Appspack will send 2 arguments, the input file and the output file.
            shPath = os.path.join(self.runPath,"setEnvAndCallRunScript.sh")
            profilePath = os.path.join(os.environ["AZORANGEHOME"], "templateProfile")
            file = open(shPath,"w")
            # This must be a tcsh script since it will call source to a template profile which uses tcsh commands
            file.write("#!/bin/tcsh\n")
            # Source the templateprofile which was created in the install procedure
            if os.path.isfile(profilePath):
                file.write("source " + profilePath + "\n")
            else:
                file.write("#source " + profilePath + " 2>&1\n")
                print "WARNING: Running APPSPACK scripts without a source file.\n        Missing "+ profilePath
            # Now that the env is OK, call the runScript with the arguments passed by appspack
            file.write("python "+os.path.join(self.runPath,"runScript.py")+" ${1} ${2}\n")
            file.close()
            os.system("chmod +x " + shPath)

            # Build the input file for appspack
            inFile = self.runPath + "input.apps"
            file = open(inFile, 'w')
            file.write("# APPSPACK INPUT FILE FOR LEARNER "+str(self.learner.name)+"\r\n")
            #Line to indicate the order of the parameter keys in the vectors
            file.write("#"+str(paramKeys)+"\r\n")
            file.write("@ \"Linear\"\r\n")
            file.write("\"Upper\" vector " + str(len(upperVector)) + " " + str(upperVector)[1:-1].replace(",", "").replace("'","").replace("\"","") + "\r\n")
            file.write("\"Lower\" vector " + str(len(lowerVector)) + " " + str(lowerVector)[1:-1].replace(",", "").replace("'","").replace("\"","") + "\r\n")
            #file.write("\"Scaling\" vector " + str(len(scalingVector)) + " " + str(scalingVector)[1:-1].replace(",", "").replace("'","").replace("\"","") + "\r\n")
            file.write("@@\r\n@ \"Evaluator\"\r\n")
            file.write("\"Executable Name\" string \"" + os.path.join(self.runPath,"setEnvAndCallRunScript.sh") + "\"\r\n")
            file.write("\"Input Prefix\" string \"" + self.runPath + "f_input\"\r\n")
            file.write("\"Output Prefix\" string \"" + self.runPath + "f_output\"\r\n")
            file.write("@@\r\n@ \"Solver\"\r\n")
            file.write("\"Debug\" int " + str(self.verbose) + "\r\n")
            file.write("\"Initial X\" vector "  + str(len(initialX)) + " " + str(initialX)[1:-1].replace(",", "").replace("'","").replace("\"","") + "\r\n")
            file.write("\"Initial F\" double " + str(initialF) + "\r\n")
            #file.write("\"Initial F\" double " + str(initialF)[1:-1].replace(",", "").replace("'","").replace("\"","") + "\r\n")
            #file.write("\"Step Tolerance\" double 0.001\r\n")
            #file.write("\"Minimum Step\" double 1.0\r\n")
            #file.write("\"Contraction Factor\" double 1\r\n")
            file.write("\"Solution File\" string \"" + self.runPath + "solution.txt\"\r\n")
            file.write("@@\r\n")
            file.close()
        except:
            if self.verbose  == 1:
                print "ERROR: Unexpected error in __CreateInput() :", sys.exc_info()[0]
            elif self.verbose > 1:
                exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
                traceback.print_exception(exceptionType, exceptionValue, exceptionTraceback,
                              limit=2, file=sys.stdout)
            return None
        return inFile



    def __EvaluateDefaultPoint(self, vars, paramKeys):
        """Evaluates the function at default point (defaults defined in the config file).
           It uses the same learner, dataset and evaluate method as the RunAppspack
           Returns the function result (float) in default point
        """
        if not self.parameters:
            return None

        # Create the python script to be called by appspack: self.runPath+"runScript.py"
        # Use the self.learnerType to add a configuration with only default values in a new file in self.runPath
        paramsConfigFile = os.path.join(self.runPath,"AZDefaultsParamsConfig.py")
        # Add to the paramsConfigFile the defaults configuration
        DefaultPars = deepcopy(self.origParameters)
        #Change the patameters to contain only defaults
        for pars in DefaultPars:
            #DefaultPars[pars][2] = "[" + DefaultPars[pars][4]  + "]"
            DefaultPars[pars][5] = False
        file = open(paramsConfigFile,'w')
        file.write(self.learnerType + "= " + str(DefaultPars)+"\r\n")
        file.close()
        if self.__CreateAppspackScript(file=self.runPath+"calcDefaultF.py", evaluateMethod=self.evaluateMethod, findMin=self.findMin, paramsConfigFile="AZDefaultsParamsConfig", paramKeys = paramKeys) == None:
            return None

        initialXF = open(self.runPath + "defaultX.txt","w")
        initialXF.write(str(len(vars))+"\r\n")
        initialXF.write(str(vars)[1:-1].replace(", ","\r\n")+"\r\n")
        initialXF.close()

        args = ["python", self.runPath + "calcDefaultF.py",self.runPath + "defaultX.txt",self.runPath + "defaultF.txt"]
        if self.verbose > 1: print "Command:",args
        exitCode = os.spawnvpe(os.P_WAIT, 'python', args, os.environ)
        if self.verbose > 1: print "Exited code:",exitCode 
        initialFF = open(self.runPath + "defaultF.txt","r")
        res = float(initialFF.readline())
        initialFF.close()

        return float(res)


    def __CreateGridSearchFiles(self, paramKeys):
        """Creates the files to be used with GridSearch Algorithm
        """
        SHScript="setEnvAndCallGSScript.sh"
        GSScript = "GSScript.py"
        if not os.path.isdir(self.runPath):
            print "No running Dir yet: " + self.runPath
            return None
        if not self.parameters:
            if self.verbose > 0: print "ERROR: missing Parameters definition"
            return None

        if not os.path.isfile(os.path.join(self.runPath,"AZLearnersParamsConfig.py")):
            if self.verbose > 0: print "ERROR: missing AZLearnersParamsConfig.py on running dir"
            return None
        #Create the shell script that will set all the envVars and then call the GSScript
        #This file will be called each time appspack want to calculate a point.
        #  Appspack will send 2 arguments, the input file and the output file.
        shPath = os.path.join(self.runPath,SHScript)
        profilePath = os.path.join(os.environ["AZORANGEHOME"], "templateProfile")
        file = open(shPath,"w")
        # This must be a tcsh script since it will call source to a template profile which uses tcsh commands
        file.write("#!/bin/tcsh\n")
        # Source the templateprofile which was created in the install procedure
        if os.path.isfile(profilePath):
            file.write("source " + profilePath + "\n")
        else:
            file.write("#source " + profilePath + " 2>&1\n")
            print "WARNING: Running APPSPACK scripts without a source file.\n        Missing "+ profilePath
        # Now that the env is OK, call the runScript with the arguments passed by appspack
        file.write("python ${1}/GSScript.py ${2} ${3}\n")
        file.close()
        os.system("chmod +x " + shPath)


        # Create the python script to be called by GridSearch
        if self.__CreateAppspackScript(file=os.path.join(self.runPath,GSScript), evaluateMethod=self.evaluateMethod, findMin=self.findMin, paramKeys = paramKeys) == None:
            return None

        return SHScript


    def __EvaluateFunction(self, vars, paramKeys):
        """Evaluates the function for the variables given in vars. 
           It uses the same learner, dataset and evaluate method as the RunAppspack
           Returns the function result (float) in vars point
        """
        if not self.parameters:
            if self.verbose > 0: print "ERROR: missing Parameters definition"
            return None

        if not os.path.isfile(os.path.join(self.runPath,"AZLearnersParamsConfig.py")):
            if self.verbose > 0: print "ERROR: missing AZLearnersParamsConfig.py on running dir"
            return None
        # Create the python script to be called by appspack: self.runPath+"runScript.py"
        #file = open(paramsConfigFile,'w')
        #file.write(self.learnerType + "= " + str(self.parameters)+"\r\n")
        #file.close()
        if self.__CreateAppspackScript(file=self.runPath+"calcInitialF.py", evaluateMethod=self.evaluateMethod, findMin=self.findMin, paramKeys = paramKeys) == None:
            return None
 
        initialXF = open(self.runPath + "initialX.txt","w")
        initialXF.write(str(len(vars))+"\r\n")
        initialXF.write(str(vars)[1:-1].replace(", ","\r\n")+"\r\n")
        initialXF.close()

        args = ["python", self.runPath + "calcInitialF.py",self.runPath + "initialX.txt",self.runPath + "initialF.txt"]
        if self.verbose > 1: print "Command:",args
        exitCode = os.spawnvpe(os.P_WAIT, 'python', args, os.environ)
        if self.verbose > 1: print "Exited code:",exitCode
        initialFF = open(self.runPath + "initialF.txt","r")
        res = float(initialFF.readline())
        initialFF.close()
       
        return float(res)

    def __RunAppspack(self, inputFile):
        """Runs the APPAPACK using the inputFile
           Returns a list with the minimum found and a list with the optimal values 
           for each parameter returned by the output APPSPACK file:
                [f , {"ParameterName1" : ParameterValue1, "ParameterName2": ParameterValue2, ...} ]
           Returns None if an error occurred while running APPSPACK
        """
        if not self.parameters:
            return None
        # Create the python script to be called by appspack: self.runPath+"runScript.py"
        if self.__CreateAppspackScript(file = self.runPath+"runScript.py", evaluateMethod = self.evaluateMethod, findMin = self.findMin) == None:
            return None
        # Call APPSPACK with input file = self.runPath+"input.apps" and wait it to finish  
        self.usedMPI = False    
        if self.advancedMPIoptions or (self.machinefile != None and (type(self.machinefile) in [types.ListType,types.StringType, types.IntType])):
            if type(self.machinefile) == types.StringType:
                if self.machinefile == "qsub":
                    machinefile = "qsub"
                elif os.path.isfile(self.machinefile):   
                    machinefile = ["-machinefile" , self.machinefile]
                else: machinefile = []
            elif type(self.machinefile) == types.ListType: # it was passed a list of nodes so we have to build the machinefile
                machinefile = ["-machinefile" , os.path.join(self.runPath,"appsMachines")]
                appsMachineFile = open(os.path.join(self.runPath,"appsMachines") , "w")
                for machine in self.machinefile:
                    appsMachineFile.write(machine+"\n")
                appsMachineFile.close()
            elif type(self.machinefile) == types.IntType: # it was passed an integer. if it is 0 then use as much cores available in local machine
                status,out = commands.getstatusoutput("cat /proc/cpuinfo | grep processor")
                self.np = len(out.split("\n"))
                machinefile = ["-machinefile" , os.path.join(self.runPath,"appsMachines")]
                appsMachineFile = open(os.path.join(self.runPath,"appsMachines") , "w")
                if self.np < 1 or status != 0:
                    self.np = 1
                appsMachineFile.write("localhost:"+str(self.np)+"\n")
                appsMachineFile.close()
            else:
                machinefile = []
                
            if not self.np or self.np < 1:
                np = []
            else:
                np = ["-np" , str(self.np)]
            if self.verbose > 0:
                otherOptions = ["-v"]
            else:
                otherOptions = []

            # set any other advanced options user might had set Ex:   "-allcpus"
            if self.advancedMPIoptions:
                otherOptions += self.advancedMPIoptions.split(" ")
            if machinefile == "qsub":
                self.__CreateQsubScript()
                args = None
                appspackExec = "qsub"
            else:
                appspackExec =  "mpirun"
                args = [appspackExec] + otherOptions + np + machinefile + [os.path.join(os.environ["AZORANGEHOME"],"orangeDependencies/bin/appspack_mpi"), self.runPath + "input.apps"]
            if self.np <= 1 and self.machinefile == 0:
                self.usedMPI = False
            else:
                self.usedMPI = True
            if self.verbose > 0: print "########### Using MPI ###########"
        else:
            self.usedMPI = False
            if self.verbose > 0: print "########### Using SERIAL ###########"
            appspackExec = "appspack_serial"
            args = [appspackExec, self.runPath + "input.apps"]

        if self.verbose > 0: print "Command: ", args
        if self.externalControl==0:
                self.finishedFlag = False
                if self.machinefile == "qsub":
                    exitCode = self.submitQsub()
                    if exitCode != 0:
                        #self.deleteTempSSHkey()
                        return None
                    else: self.waitForQsub()
                    self.finishedFlag = True
                    #self.deleteTempSSHkey()
                else:
                    if self.verbose > 1: print "Command:",args
                    exitCode = os.spawnvpe(os.P_WAIT, appspackExec, args, os.environ)
                    if self.verbose > 1: print "Exited code:",exitCode
                    self.finishedFlag = True
                    if exitCode != 0:
                        return None
        else:
                if self.machinefile == "qsub":
                    exitCode = self.submitQsub()
                    ##scPA   Check if the Qsub command did not fail. If it did fail, report the job as terminated
                    if exitCode != 0:
                        #self.deleteTempSSHkey()
                        self.qsubJobId = None
                        self.finishedFlag = True
                        return None
                    #ecPA
                    time.sleep(5)  # Assure that the qsub job had time to start
                    isRunning = self.getIsQsubRunning()
                    if not isRunning:
                        #self.deleteTempSSHkey()
                        return None
                    else:
                        self.finishedFlag = False
                else:
                    if not self.isFinished():
                        return None
                    if self.verbose > 1: print "Command:",args
                    self.appspackPID = os.spawnvpe(os.P_NOWAIT, appspackExec , args, os.environ)
                    if self.verbose > 1: print "Running PID:",self.appspackPID
                    if self.appspackPID == 0:
                        return None
                    else:
                        self.finishedFlag = False
        return True

    def submitQsub(self):
        """
        Submit with qsub and get the jobID
        """
        # Assure that the user can log on to all required nodes without passwd
        #self.createTempSSHkey()

        #This is to be some where else...
        #miscUtilities.autoValidateRSAKey(machine,user)

        # Assess the memory requirements
        memSize = dataUtilities.getApproxMemReq(self.dataSet)

        presentDir = os.getcwd()
        os.chdir(self.runPath)
        if self.verbose > 1: print "Command:  qsub -l mf="+str(memSize)+"M -q "+self.queueType+" "+self.qsubFile
        exitCode, qsub = commands.getstatusoutput("qsub -l mf="+str(memSize)+"M -q "+self.queueType+" "+self.qsubFile)
        os.chdir(presentDir)
        if self.verbose > 1: print "Exit Status: ",exitCode,"\n",qsub
        try: self.qsubJobId = string.split(qsub)[2]
        except: exitCode = 1
        return exitCode

    def createTempSSHkey(self):
        """
        Deprecated - Run without passwd between nodes
        """
        # Assure that authorized_keys, known_hosts and id_dsa.pub has 644 permisions
        os.umask(022)
        #scPA   User home dir may be other than /honme/<user>
        #user = os.environ["USER"]
        homeDir = os.environ["HOME"]  #"/home/"+user
        #ecPA
        self.sshDir = homeDir+"/.ssh"
        self.sshSaveDir = homeDir+"/.sshSave"+str(time.time())
        os.system("/bin/mv "+self.sshDir+" "+self.sshSaveDir)
        exitCode, stdOut = commands.getstatusoutput('yes | ssh-keygen -t dsa -f '+self.sshDir+'/id_dsa -q -N ""')
        os.system("cat "+self.sshDir+"/id_dsa.pub > "+self.sshDir+"/authorized_keys")
        

    def deleteTempSSHkey(self):
        """
        Deprecated - Put back original ssh dir if the temp dir is still there
        """
        if os.path.exists(self.sshSaveDir):
            os.system("/bin/rm -r "+self.sshDir)
            os.system("/bin/mv "+self.sshSaveDir+" "+self.sshDir)


    def waitForQsub(self):

        isRunning = self.getIsQsubRunning()
        while isRunning:
            time.sleep(2)
            isRunning = self.getIsQsubRunning()


    def getIsQsubRunning(self):

        exitCode, qstat = commands.getstatusoutput("qstat")
        if self.qsubJobId: 
            idx = string.find(qstat, self.qsubJobId)
            if idx != -1: isRunning = True
            else: isRunning = False
        else: isRunning = False
        return isRunning


    def __CreateQsubScript(self):
        """
        Creates a script to execute appspack_mpi with qsub.
        Returns list of arguments to the os.spawnvpn call.
        """

        full_hostfile = os.path.join(self.runPath, "full_hostfile")
        full_machinefile = os.path.join(self.runPath, "full_machinefile")
        appsFile = os.path.join(self.runPath, "input.apps")
        execFile = os.path.join(os.environ["AZORANGEHOME"],"orangeDependencies/bin/appspack_mpi")
        self.qsubFile = os.path.join(self.runPath, "runQsub.sh")

        fid = open(self.qsubFile, "w")
        writeStr = """#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N appspack_mpi
#$ -pe pe_mpich """+str(self.np)+"""
#$ -v MPICH_HOME=/opt/az/anlg/mpich/1.2.7p1

exec 2>&1
env | sort
echo Running on host: `hostname`
echo Fix machine file
perl -pe 's/ /:/g' """+full_hostfile+""" > """+full_machinefile+"""
echo "Start mpirun"
$MPICH_HOME/bin/mpirun -machinefile """+full_machinefile+""" -np """+str(self.np)+""" """+execFile+""" """+appsFile+"""

echo "end mpirun"
        """
        fid.write(writeStr)
        fid.close()
 
        # Assure that the qsub script can be executed
        os.system("chmod a+x "+self.qsubFile)


    def __CreateAppspackScript(self, file, evaluateMethod, findMin, paramsConfigFile = "AZLearnersParamsConfig", paramKeys = "None"):
        """
        Creates the script to be called by appspack on each evaluation of the function in a specific point
        The point being evaluated is passed by appspack on the input file of parameter 1 and the output 
        (the image of the input variables) must be placed in the output file specified by appspack on
        second parameter
        It uses the globals:
                self.runPath
                self.learner
                self.learnerType
                self.dataSet
        """
        
        if os.path.exists(file):
            os.remove(file)

        fullLearnerClass = str(self.learner.__class__).split("'")
        if len(fullLearnerClass)>1:
            fullLearnerClass = fullLearnerClass[1]
        else:
            fullLearnerClass = fullLearnerClass[0]

        if (fullLearnerClass.find(".")==-1):
            fullLearnerClass = "orange." + fullLearnerClass
        if findMin:
            sign = ""
        else:
            sign = "-"

        if self.samplingMethod == self.LEAVE_ONE_OUT:
           sMethod = "orngTest.leaveOneOut([learner], dataSet)"
        elif self.nExtFolds:    ## by default: samplingMethod = self.FOLD_CROSS_VALIDATION
           sMethod = "orngTest.crossValidation([learner], dataSet, folds=" + str(self.nFolds) + ", strat=orange.MakeRandomIndices.StratifiedIfPossible, randomGenerator = MyRandom)"
        else:
           sMethod = "orngTest.crossValidation([learner], dataSet, folds=" + str(self.nFolds) + ", strat=orange.MakeRandomIndices.StratifiedIfPossible)"
        
        #Define the scriptVars to use when filling the script model file 
        scriptVars={}
        scriptVars["evalMethod"]=evaluateMethod[:evaluateMethod.rfind(".")]
        scriptVars["fullLearner"]=fullLearnerClass[:fullLearnerClass.rfind(".")]
        scriptVars["FullLearnerClass"]=fullLearnerClass
        scriptVars["learnerType"]=self.learnerType
        scriptVars["dataset"]=self.dataSet
        scriptVars["sMethod"]=sMethod
        scriptVars["evalMethodFunc"]=evaluateMethod
        scriptVars["runPath"]=self.runPath
        scriptVars["resSign"]=sign
        scriptVars["verbose"]=self.verbose
        scriptVars["paramsConfigFile"]=paramsConfigFile
        scriptVars["paramKeys"]=paramKeys
        scriptVars["nFolds"]=self.nFolds
        scriptVars["nExtFolds"]=self.nExtFolds
        scriptVars["machinefile"]=self.machinefile
        
        #open the script model file which will be filled with the scriptVars{}
        modelFile=open(os.path.dirname(__file__)+"/OptScriptModel.py")
        strFile = ""
        for line in modelFile.readlines(): strFile+=line
        modelFile.close()
        strFile=strFile % scriptVars
        
        #Write the model file already filled into the  file on running path
        try:
            pyFile = open(file,"w")   
            pyFile.write(strFile)
            pyFile.close()    
        except:
            return None
        return True


def optimize(learner, learnerName, trainDataFile, responseType, verbose = 0, queueType = "batch.q", runPath = None, nExtFolds = None, nFolds = 5, useGrid = False):
    """
    Optimize with default values for the learner (defined in AZLearnersParmsConfig)
    Run optimization in parallel.
    Possible values:
    responseType: Classification or other
    learnerNames: RFLearner, CvSVMLearner, CvANNLearner, PLSLearner
    queueType: 'batch.q' or 'quick.q' (jobs start immediatly but are terminated after 30 min) or 'NoSGE'
    runPath: If directory not provided, will run in NFS_SCRATCHDIR
    """

    # Create an interface for setting optimizer parameters
    pars = AZLearnersParamsConfig.API(learnerName)

    # Create a directory for running the appspack (if not defined it will use the present working directory)
    if not runPath:
        runPath = miscUtilities.createScratchDir(desc ="optQsubTest", baseDir = AZOC.NFS_SCRATCHDIR)

    # Response type
    if responseType == "Classification":
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False
    else:
        evalM = "AZutilities.evalUtilities.RMSE"
        fMin = True

    # Distributed computing available?
    if queueType == "NoSGE":
        np = None
        machineFile = None
    else:
        np = 8
        machineFile = "qsub" 


    # Calculate the optimal parameters. This can take a long period of time!
    optimizer = Appspack()
    tunedPars = optimizer(learner=learner,\
                    dataSet=trainDataFile,\
                    evaluateMethod = evalM,\
                    useParameters = pars.getParametersDict(),\
                    findMin=fMin,\
                    runPath = runPath,\
                    verbose = verbose,\
                    nExtFolds = nExtFolds, \
                    nFolds = nFolds, \
                    useGridSearchFirst = useGrid,\
                    gridSearchInnerPoints = 3,\
                    np = np,\
                    machinefile = machineFile,\
                    queueType = queueType)

    #if verbose > 0:
    print "====================== optimization Done ==========================="
    print "Learner optimized flag = ", learner.optimized
    print "Tuned parameters = ", tunedPars[1]
    print "Best optimization result = ", tunedPars[0]
    print "check the file optimizationLog.txt in "+runPath+" to see the intermediate results of optimizer!"

    return learner, learner.optimized


def optimizeSelectedParam(learner, learnerName ,trainDataFile, paramList, responseType, grid = False, useGrid = False, verbose = 0, queueType = "batch.q", runPath = None, nExtFolds = None, nFolds = 5):
    """
    Optimize the parameters in paramList. Run optimization in parallel.
    Possible values:
    grid: Execute on the SGE or not
    responseType: Classification or other
    learnerNames: RFLearner, CvSVMLearner, CvANNLearner, PLSLearner
    queueType: batch.q or quick.q (jobs start immediatly but are terminated after 30 min)
    runPath: If directory not provided, will run in NFS_SCRATCHDIR
    """

    optimizer = Appspack()

    # Create an interface for setting optimizer parameters
    pars = AZLearnersParamsConfig.API(learnerName)

    # Set all parameters to not be optimized
    pars.setOptimizeAllParameters(False)

    # Set the parameters in parameterList to be optimized
    for parameter in paramList:
        pars.setParameter(parameter,"optimize",True)

    # Create a directory for running appspack (if not defined it will use the present working directory)
    if not runPath:
        runPath = miscUtilities.createScratchDir(desc ="optQsubTest", baseDir = AZOC.NFS_SCRATCHDIR)

    if responseType == "Classification":
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False
    else:
        evalM = "AZutilities.evalUtilities.RMSE"
        fMin = True

    if grid:
        machinefile = "qsub"
        np = 8
    else:
        machinefile = None
        np = None

    # Calculate the optimal parameters. This can take a long period of time!
    tunedPars = optimizer(learner=learner,\
                    dataSet=trainDataFile,\
                    evaluateMethod = evalM,\
                    useParameters = pars.getParametersDict(),\
                    findMin=fMin,\
                    runPath = runPath,\
                    useGridSearchFirst = useGrid,\
                    gridSearchInnerPoints = 3,\
                    nExtFolds = nExtFolds,\
                    nFolds = nFolds,\
                    np = np,\
                    machinefile = machinefile,\
                    verbose = 0,\
                    queueType = queueType)

    if verbose > 0:
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "Results directory"
        print runPath
        print "check the file optimizationLog.txt to see the intermediate results of optimizer!"

    return learner



