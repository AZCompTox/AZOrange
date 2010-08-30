import os
from AZutilities import miscUtilities
import time

class GridSeach:
    def __init__(self, **kwds):
        #Possible user defined Vars
        self.varsUpperLimits = []
        self.varsLowerLimits = []
        self.nInnerPoints    = 5
        self.scriptFilesPath = None
        self.ScriptFileName  = None
        self.scriptRunPath = None
        self.verbose = 1
        self.runningPaths = []
        self.AllPaths = []
        # Append arguments to the __dict__ member variable 
        self.__dict__.update(kwds)


        # Non-User defined vars

    def __call__(self,  **kwds):
            # Append arguments to the __dict__ member variable 
            self.__dict__.update(kwds)
            self.runningPaths = []
            self.AllPaths = []
            initialX = len(self.varsUpperLimits)*[0]
            #number of inner points to divide each variable
            if self.nInnerPoints < 1:
                self.nInnerPoints = 1
            if (len(self.varsUpperLimits) != len(self.varsLowerLimits)) or (len(self.varsUpperLimits)==0):
                print "Vaector Variable limits not correcly defined:"
                print "Upper Limit: ", self.varsUpperLimits
                print "Lower Limit: ", self.varsLowerLimits
                return None
            for idx,upVal in enumerate(self.varsUpperLimits):
                if self.varsLowerLimits[idx] > upVal:
                    print "Variable limits not correcly defined."               
                    print "Upper Limit: ", self.varsUpperLimits
                    print "Lower Limit: ", self.varsLowerLimits
                    return None
            #Create all the ranges:  [ [Attr1Lo,Attr1.1, Attr1.2,...]  ,  []  ...  [AttrNUp, AttrN.1Up,...,AttrN.NUp] ]
            allRanges = [miscUtilities.Range(x[0],x[1],(x[1]-x[0])/(self.nInnerPoints+1.0)) for x in zip(self.varsLowerLimits,self.varsUpperLimits)]
            allPoints = [[0]*len(initialX)]*((len(allRanges[0])-2)**len(initialX))

            #product(Npoints,nAttrs)  
            allPointsIdx = self.__product(range(len(allRanges[0]))[1:-1],repeat = len(initialX))
            for (idx,pointIdx) in enumerate(allPointsIdx):
                 allPoints[idx] = [allRanges[attr][Pidx] for (attr,Pidx) in enumerate(pointIdx)]

            #If there are no points to be evaluated, evaluate at least the Upper limit point so that we have something to return.
            if allPoints == []:
                allPoints = [initialX]
                print "WARNING: AZGridSearch: It was not possible to compute the gridSearch points. Evaluating the Upper Limit point only:"
                print "Upper Limit: ", self.varsUpperLimits
                print "Lower Limit: ", self.varsLowerLimits

            #Launch system call to eavaluate each point. No Wait for it to finish. We will collect the results later
            for idx,point in enumerate(allPoints):
                self.__LaunchPointCalc(str(idx), point)

            allRes = [None]*len(allPoints)
            #Pool each 1 second to check if all "threads" already finished
            while len(self.runningPaths) > 0:
                for path in self.runningPaths:
                    idx = int(os.path.split(path)[-1].split("_")[-1])
                    resPath = os.path.join(path,"resPoint.txt")
                    if os.path.isfile(resPath):
                        File = open(resPath,"r")
                        allRes[idx] = eval(File.readline().strip())
                        File.close()
                        self.runningPaths.pop(self.runningPaths.index(path))
                time.sleep(1)
              
            return {"varValues":allPoints, "results":allRes, "nPoints":len(allPoints), "nFailedPoints":allRes.count(None)} 

    def __product(self, *args, **kwds):
            # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
            # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
            pools = map(tuple, args) * kwds.get('repeat', 1)
            result = [[]]
            for pool in pools:
                result = [x+[y] for x in result for y in pool]
            for prod in result:
                yield tuple(prod)

    def __LaunchPointCalc(self, suffix, vars):
        #Create a unique running dir for this specific point
        runPath = os.path.join(self.scriptRunPath, "runPoint_"+suffix)
        os.system("rm -rf " + runPath) 
        os.system("mkdir -p " + runPath)
        #Copy all necessary files to the point running dir
        os.system("cp -rf " + os.path.join(self.scriptFilesPath,"*") + " " +runPath)
        self.AllPaths.append(runPath)
        self.runningPaths.append(runPath)
        
        #Create the input file for this specific point
        initialXF = open(os.path.join(runPath , "varsPoint.txt"),"w")
        initialXF.write(str(len(vars))+"\r\n")
        initialXF.write(str(vars)[1:-1].replace(", ","\r\n")+"\r\n")
        initialXF.close()
        
        #Launch the thread. Don't wait for it to finish. We will check this dir periodically looking for 
        # the rersults file.
        args = [os.path.join(runPath, self.ScriptFileName),\
                os.path.realpath(runPath),\
                os.path.join(runPath,"varsPoint.txt"),\
                os.path.join(runPath, "resPoint.txt")]
        if self.verbose > 1: print "Command:",args
        exitCode = os.spawnvpe(os.P_NOWAIT, os.path.join(runPath, self.ScriptFileName), args, os.environ)
        if self.verbose > 1: print "Exited code:",exitCode

        return True


