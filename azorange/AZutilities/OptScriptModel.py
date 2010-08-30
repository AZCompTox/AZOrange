# File automatically created by  paramOptUtilities.py
import os, time,random,sys, types
import orange, orngTest
import %(evalMethod)s
import %(fullLearner)s
import %(paramsConfigFile)s
from AZutilities import miscUtilities
from AZutilities import dataUtilities
from AZutilities import sgeUtilities
from math import sqrt
from math import floor
import statc
from glob import glob


# Read all parts of the intRes.txt file. 
def readIntRes():
    intRes= []
    for part in sorted(glob("%(runPath)s*intRes.txt")):
            file = open(part,"r")
            intRes = intRes + file.readlines()
            file.close()
    return intRes

version = 9
verbose = %(verbose)s
inputFile = sys.argv[1]
outputFile = sys.argv[2]
learner = %(FullLearnerClass)s()

useDefaults = False
inF = open(inputFile,"r")
if "defaultX" in inputFile:
    useDefaults = True
    #These vars are not used at all for dafaul point. they will be just used to confirn the number of parameters to optimize
    #Vars are also used to create the intRes file as "asked by appspack"
    vars = [str(x).strip() for x in inF.readlines()][1:]
else:  
    vars = [types.FloatType(x) for x in inF.readlines()][1:]
inF.close()

# All Learner's parameters from config file
parameters = %(paramsConfigFile)s.%(learnerType)s
dataSet=dataUtilities.DataTable("%(dataset)s")
N_ATTR = len(dataSet.domain.attributes)
N_EX = len(dataSet) - floor(len(dataSet)/%(nFolds)s)

if dataSet.domain.classVar.varType == orange.VarTypes.Discrete:
    isClassifier = True 
else:
    isClassifier = False
    
#Parameter names to be optimized (sent directly or loaded ahead from input.apps) 
paramKeys = %(paramKeys)s

try:
    if paramKeys == None:
        if not os.path.isfile("%(runPath)sinput.apps"):
            if verbose > 0: print "ERROR: Cannot find the correspondence parameters between names and values! No input.apps file!"
            sys.exit()
        else:
            inputApps = open("%(runPath)sinput.apps")
            paramKeys = eval(inputApps.readlines()[1][1:-2])
            inputApps.close()
    if (not paramKeys) or (type(paramKeys) != types.ListType) or (len(paramKeys)!=len(vars)):
        if verbose > 0: print "ERROR: Cannot find the correspondence parameters between names and values! Wrong parameter keys!"
        sys.exit() 
    for key in paramKeys:
        if key not in parameters:
            if verbose > 0: print "ERROR: Cannot find the correspondence parameters between names and values! Keys are different!"
            sys.exit()
     
except:
    if verbose > 0: print "ERROR: Cannot find the correspondence parameters between names and values!"
    sys.exit()


#Check if the defaults sent are according to the defaults from config file
if useDefaults:
    for idx,param in enumerate(paramKeys):
        if parameters[param][4] != vars[idx]:
            if verbose > 0 and str(parameters[param][4]).replace("'","").replace('"','') != str(vars[idx]).replace("'","").replace('"',''):
                print "WARNING!!  Default parameter "+param+" = "+vars[idx]+" was not the same as the default from the config file: "+str(parameters[param][4])
            vars[idx] = parameters[param][4]
testParameters = {}
#Fill the TestParameters to be used in the Learner in this particular step (values are required by the appspack -> "vars" variable or by the defsultX.txt file which send the values in text format already 'decoded' !!)
if len(glob("%(runPath)s*intRes.txt")):
    #resFile=open("%(runPath)sintRes.txt","r")
    #resFile=miscUtilities.lockFile("%(runPath)sintRes.txt","r")
    #headerL = resFile.readline()
    headerL = readIntRes()[0]
    #resFile.close()
    keys =  [attr.strip() for attr in headerL.split("\t") if attr.strip()[0:5] != "APPS_" and attr.strip()[0:8] != "LEARNER_" and attr.strip() != "EVAL_RES"]
else:
    keys = [key for key in parameters]
for key in keys:
    valuesRange=eval(parameters[key][2])
    if key in paramKeys:
        idx = paramKeys.index(key)
    else:
        idx = None
    if parameters[key][1]=="values":
        if type(valuesRange[0]) == types.StringType:
            if useDefaults or not parameters[key][5]:
                testParameters[key] = str(parameters[key][4])
            elif len(valuesRange) == 1:
                testParameters[key] = str(valuesRange[0])
            else:
                testParameters[key] = str(valuesRange[types.IntType(round(vars[idx]))])
            if type(eval(parameters[key][0])) == types.ListType:
                if eval(parameters[key][0])[0] != types.TypeType:
		    if eval(parameters[key][0])[0] == types.BooleanType:
                        testParameters[key] = eval(parameters[key][0])[0](eval(testParameters[key]))
		    else:
                        testParameters[key] = eval(parameters[key][0])[0](testParameters[key])
            else:
                if eval(parameters[key][0]) != types.TypeType:
		    if eval(parameters[key][0]) == types.BooleanType:
                        testParameters[key] = eval(parameters[key][0])(eval(testParameters[key]))
		    else:
                        testParameters[key] = eval(parameters[key][0])(testParameters[key])
        
        else:
            if (len(valuesRange) != 1) and parameters[key][5] and (not useDefaults):
                varVal = int(round(vars[idx]))
            if type(eval(parameters[key][0])) == types.ListType:
                if useDefaults or not parameters[key][5]:
                    if eval(parameters[key][0]) == types.BooleanType:
                        testParameters[key] = eval(parameters[key][0])[0](eval(parameters[key][4]))
                    else:
                        testParameters[key] = eval(parameters[key][0])[0](parameters[key][4])
                elif len(valuesRange) == 1:
                    testParameters[key] = eval(parameters[key][0])[0](valuesRange[0])
                else:
                    testParameters[key] = eval(parameters[key][0])[0](valuesRange[varVal])
            else:
                if useDefaults or not parameters[key][5]:
                    if eval(parameters[key][0]) == types.BooleanType:
                        testParameters[key] = eval(parameters[key][0])(eval(parameters[key][4]))
                    else:
                        testParameters[key] = eval(parameters[key][0])(parameters[key][4])
                elif len(valuesRange) == 1:
                    testParameters[key] = eval(parameters[key][0])(valuesRange[0])
                else:
                    testParameters[key] = eval(parameters[key][0])(valuesRange[varVal])
    else:
        if type(eval(parameters[key][0])) == types.ListType: 
            if useDefaults or  not parameters[key][5]:
                if eval(parameters[key][0])[0] == types.BooleanType:
                    testParameters[key] = eval(parameters[key][0])[0](eval(parameters[key][4]))
                else:
                    testParameters[key] = eval(parameters[key][0])[0](parameters[key][4])
            elif len(valuesRange) == 1:
                testParameters[key] = eval(parameters[key][0])[0](valuesRange[0])
            else:  
                testParameters[key] = eval(parameters[key][0])[0](vars[idx])
        else:
            if useDefaults or  not parameters[key][5]:
                if eval(parameters[key][0]) == types.BooleanType:
                    testParameters[key] = eval(parameters[key][0])(eval(parameters[key][4]))
                else:
                    testParameters[key] = eval(parameters[key][0])(parameters[key][4])
            elif len(valuesRange) == 1:
                testParameters[key] = eval(parameters[key][0])(valuesRange[0])
            else:  
                testParameters[key] = eval(parameters[key][0])(vars[idx])

#Test if this point was already required by appspack since we are emulating the use of discrete variables in appspack
evaluated = None
actualIntRes = []
actualIntResHeader = []
if len(glob("%(runPath)s*intRes.txt")):
        #Load the intRes.txt built so far
        #file = open("%(runPath)sintRes.txt")
        #file = miscUtilities.lockFile("%(runPath)sintRes.txt") 
        #actualIntRes = file.readlines() 
        actualIntRes = readIntRes()
        #print actualIntRes
        #file.close()
        actualIntResHeader = actualIntRes[0].split()
        parValues = []
        #Built the par values in this particular case
        for key in keys:
            if type(eval(parameters[key][0])) == types.TypeType and eval(parameters[key][0]) == types.TypeType:
                parValue= "'"+str(testParameters[key]).replace("'","")+"'"
            else:
                parValue = str(testParameters[key]).replace("'","")
            parValues.append(parValue)
        #Test the par values in this case with the parameter of each evaluation in the intRes file
        for ev in actualIntRes[1:]:
#            if len(actualIntResHeader) == len(ev.split()):
                evValues = [([x for x in ev.split()])[actualIntResHeader.index(key)] for key in keys]
                if parValues == evValues:
                    evaluated = ev.split()[-1].strip()
                    break 

#If the point was already evaluated, exit and return proper value in output file
if evaluated:
        outF = open(outputFile,"w")
        outF.write("%(resSign)s" + str(evaluated) + "\r\n")
        outF.close()

        if verbose > 0:
            print "\r\nVars from inputFile (",inputFile,")= ", vars
            print "Point was already evaluated:"
            for parameter in testParameters:
                print "\t", parameter,"= ", testParameters[parameter]
            print "Solution already found = "+"%(resSign)s" + str(evaluated)
        sys.exit()

#Assign TestParameters to the Learner
#Types handle here
for key in keys:
    if hasattr(learner, "setattr"):
        if type(eval(parameters[key][0])) == types.ListType:
            learner.setattr(key, [testParameters[key]])
        elif  type(eval(parameters[key][0])) == types.TypeType and eval(parameters[key][0]) == types.TypeType:
            learner.setattr(key, eval(testParameters[key]))
        else:
            learner.setattr(key, testParameters[key])
    else:
        if type(eval(parameters[key][0])) == types.ListType:
            setattr(learner, key, [testParameters[key]])
        elif type(eval(parameters[key][0])) == types.TypeType and eval(parameters[key][0]) == types.TypeType:
            setattr(learner, key,  eval(testParameters[key]))
        else:
            setattr(learner, key, testParameters[key])

# Evaluate function at specified point
if %(nExtFolds)s:
    evalResList = []
    if verbose > 0: print "Number of external folds"
    if verbose > 0: print %(nExtFolds)s
    if "%(machinefile)s" == "qsub":
        jobScript = """\
import orange,orngTest,random,pickle,os

paramFile=open("Params.pkl","r")
(learner,nFolds,dataSet,evaluateMethod) = pickle.load(paramFile)
paramFile.close()
MyRandom = orange.RandomGenerator(1000*int(os.environ["SGE_TASK_ID"]))
res = %(sMethod)s
print pickle.dumps(evaluateMethod(res)[0])
"""
        # Assess the memory requirements
        memSize = dataUtilities.getApproxMemReq(dataSet)

        evalResList = sgeUtilities.arrayJob(jobName = "EvalJob", jobNumber = %(nExtFolds)s, jobParams = [learner,%(nFolds)s,dataSet,%(evalMethodFunc)s], jobQueue = "batch.q", jobScript = jobScript, memSize = str(memSize)+"M")
    else:
        for idx in range(%(nExtFolds)s):
            MyRandom = orange.RandomGenerator(1000*idx+1)
            res = %(sMethod)s
            evalResList.append(%(evalMethodFunc)s(res)[0])

    if isClassifier:
        evalRes = [round(statc.mean(evalResList),3)]
    else:
        evalRes = [round(statc.mean(evalResList),2)]
    if verbose > 0: print evalRes
else:
    res = %(sMethod)s
    evalRes = %(evalMethodFunc)s(res)

# Save intermediate result
#if os.path.exists("%(runPath)sintRes.txt"):
if [os.path.basename(f) for f in glob("%(runPath)s"+"*intRes.txt")] != []:
    tmpNew=False
else:
    tmpNew=True
#tmp=miscUtilities.lockFile("%(runPath)sintRes.txt","a")
#tmp=open("%(runPath)sintRes.txt","a")
fname=str(time.time()) + "-" + str(random.randint(0,10000)) + "_intRes.txt"
tmp=open("%(runPath)s"+fname,"w")

# Header of results table file
if tmpNew:
    tmp.write(\
        str([x for x in keys])[1:-1].replace(", ","\t").replace("'","")+"\t"+               # parameters sent to the Learner \
        str(["APPS_"+x for x in keys])[1:-1].replace(", ","\t").replace("'","") +"\t"+      # parameters required by the appspack\
        str(["LEARNER_"+x for x in keys])[1:-1].replace(", ","\t").replace("'","") +"\t"+   # Effective Learner parameters\
        "EVAL_RES"+"\r\n")


#Write the parameters used at this specific point
for key in keys:
    if type(eval(parameters[key][0])) == types.TypeType and eval(parameters[key][0]) == types.TypeType:
        tmp.write("'"+str(testParameters[key]).replace("'","").replace(" ","")+"'\t")
    else:
        tmp.write(str(testParameters[key]).replace("'","").replace(" ","")+"\t")

#Write the appspack vars
for key in keys:
    if key in paramKeys:
        tmp.write(str(vars[paramKeys.index(key)]).replace(" ","")+"\t")
    else:
        tmp.write("NA\t")

#Write the effective learners paramaters
for key in keys:
    if hasattr(learner, key):
        tmp.write(str(getattr(learner,key)).replace(" ","")+"\t")
    else:
        tmp.write("NA\t")


#write result
tmp.write(str(evalRes[0])+"\r\n")
tmp.close()


outF = open(outputFile,"w")
outF.write(str(%(resSign)sevalRes[0]) + "\r\n")
outF.close()

if verbose > 0:
    print "\r\nVars from inputFile (",inputFile,")= ", vars
    print "Evaluated for:"
    for parameter in testParameters:
        print "\t", parameter,"= ", eval("learner." + parameter)
    print "Solution = %(resSign)s" + str(evalRes[0])


