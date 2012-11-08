from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import cPickle,commands,os
from glob import glob
from shutil import copy


def arrayJob(jobName = "AZOarray",jobNumber =1 ,jobParams = [], jobParamFile = "Params.pkl", jobQueue = "quick.q", jobScript = "", memSize = "150M", environSource = os.path.join(os.environ["AZORANGEHOME"],"templateProfile.bash")):   

        runPath = miscUtilities.createScratchDir(desc ="optQsub"+jobName, baseDir = AZOC.NFS_SCRATCHDIR)
        cwd = os.getcwd()
        os.chdir(runPath)

        paramFile = open(jobParamFile,"w")
        cPickle.dump(jobParams,paramFile)
        paramFile.close()
 
        jobFile = open(jobName + ".py","w")
        jobFile.write(jobScript)
        jobFile.close()

        if environSource:
            sgeRunFile = os.path.join(runPath, "sgeRun") 
            envSetFile = open(sgeRunFile,"w")
            envSetFile.write("#!/bin/sh\n\n")    
            envSetFile.write("source " + environSource+"\n")    
            envSetFile.write("python  " + os.path.join(runPath, str(jobName) + ".py")+"\n")    
            envSetFile.close()
            os.system("chmod a+x " + sgeRunFile)

            cmd = "echo " + sgeRunFile + \
              " | qsub -cwd -q " + str(jobQueue) + \
              " -p -800 -t 1-" + str(jobNumber) + \
              " -N " + str(jobName) + \
              " -S /bin/sh -sync yes" + \
              AZOC.SGE_QSUB_ARCH_OPTION_CURRENT + \
              " -l mf=" + str(memSize) # specify shell /bin/sh so not to get warning: no access to tty in output file.
        else:
            cmd = "echo python " + os.path.join(runPath, str(jobName) + ".py") + \
              " | qsub -cwd -V -q " + str(jobQueue) + \
              " -p -800 -t 1-" + str(jobNumber) + \
              " -N " + str(jobName) + \
              " -S /bin/sh -sync yes" + \
              AZOC.SGE_QSUB_ARCH_OPTION_CURRENT + \
              " -l mf=" + str(memSize) # specify shell /bin/sh so not to get warning: no access to tty in output file.


        (status, output) = commands.getstatusoutput(cmd)

        # Check exit status of all our jobs
        print "runDir: ",runPath
        if status != 0:
            print jobName + " failed! Code = " + str(status)
            print output
            raise ValueError
        for line in output.split("\n"):
          if not "exit code 0" in line:
            if not "Your job-array" in line:
                print jobName + " failed! " + line
                raise ValueError

        # Check if error files exist that are not empty.
        #   Pedro Almeida: We cannot relie on this since some warnings are writen to the error output...
        #for part in sorted(glob(os.path.join(runPath,jobName+".e*"))):
        #    if os.path.getsize(part) != 0:
        #        print jobName + " failed! file " + str(part)
        #        raise ValueError

        # Build result list from pickle objects
        resList = []
        for part in sorted(glob(os.path.join(runPath,"RES_out*.pkl"))):
            file = open(part,"r")
            resList.append(cPickle.load(file))
            file.close()
        
        os.chdir(cwd)
        miscUtilities.removeDir(runPath)
        return resList

class JobError(IOError):
        """JobError"""

class Job:

    def __init__(self,name = "AZOrangeJob", queue = "batch.q", priority = None, range = None, wd = "temporary", env = "full", shell = "/bin/sh", sync = "yes", resources = None, body = "", params = None , action = None, pyFile = None , hold = None):

        qsub = "#!" + shell                     + "\n\n"
        qsub = qsub + "#$ -N " + name           + "\n"
        self.name = name
        qsub = qsub + "#$ -q " + queue          + "\n"
        self.queue = queue
        self.action = action
        if priority:
            qsub = qsub + "#$ -p " + str(priority)      + "\n"
            self.priority = priority
        if wd == "temporary":
            qsub = qsub + "#$ -cwd"             + "\n"
            self.wd = miscUtilities.createScratchDir(desc ="optQsub"+name, baseDir = AZOC.NFS_SCRATCHDIR)
        else:
            qsub = qsub + "#$ -wd " + wd        + "\n"
            self.wd = wd
        if env == "full":
            qsub = qsub + "#$ -V"        + "\n"
        if sync == "yes":
            qsub = qsub + "#$ -sync yes"        + "\n"
            self.sync = True
        qsub = qsub + "#$ -S " + shell          + "\n"
        self.shell = shell   
        self.range = range
        qsub = qsub + "#$ " + AZOC.SGE_QSUB_ARCH_OPTION_CURRENT + "\n"
        if range:
            qsub = qsub + "#$ -t 1-" + str(range)      + "\n"
        if resources:
            qsub = qsub + "#$ -l " + resources  + "\n"
            self.resources = resources
        if hold:
            qsub = qsub + "#$ -hold_jid " + str(hold)   + "\n"
        if body:
            qsub = qsub + "\n" + body + "\n"
        else:
            qsub = qsub + "\npython " + self.name + ".py\n"
        self.qsubFile = os.path.join(self.wd,name) + ".sh"
        self.qsubScript = qsub
        if pyFile:
            self.textFile(os.path.join(self.wd,name) + ".py",pyFile)
        self.textFile(self.qsubFile,qsub)
        if params:
            paramFile = open(os.path.join(self.wd,name) + ".params","w")
            cPickle.dump(params,paramFile)
            paramFile.close()

    def copyFile(self,file):
	copy(file,self.wd)

    def textFile(self,filename,text):
        file = open(os.path.join(self.wd,filename),"w")
        file.write(text)
        file.close()

    def run(self):
        cwd = os.getcwd()
        os.chdir(self.wd)
        (self.status, self.output) = commands.getstatusoutput("qsub " + self.qsubFile)

        if self.status != 0:
            print self.name + " failed! Code = " + str(self.status)
            print self.output
            raise JobError

        #   Pedro Almeida: We cannot relie on this since some warnings are writen to the error output...
        #for part in sorted(glob(os.path.join(self.wd,self.name+".e*"))):
        #    if os.path.getsize(part) != 0:
        #        print self.name + " failed! file " + str(part)
        #        raise JobError

        if self.range:
            for line in self.output.split("\n"):
              if not "exit code 0" in line:
                if not "Your job-array" in line:
                    print self.name + " failed! " + line
                    raise JobError

        if self.action == "picklebuild":  # build return array from pickle objects in output file
                resList = []
                for part in sorted(glob(os.path.join(self.wd,"RES_out*.pkl"))):
                    file = open(part,"r")
                    resList.append(cPickle.load(file))
                    file.close()
                return resList

    def __call__(self):
	return self.run()        

    def info(self):
	print "Name: " + self.name
        print "Dir: " + self.wd 
        print "Shell: " + self.shell

        cwd = os.getcwd()
        os.chdir(cwd)
#        miscUtilities.removeDir(runPath)

