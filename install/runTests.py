import ConfigParser
import os, sys
import time
import types
import commands
import string
import re 

class tester:
    def __init__(self,configFile):
	if not os.path.isfile(configFile):
	    print "Missing the configuration file. "+configFile+" not found.\n You may use the setupTemplate.ini file to create your setup.ini file."
	    self.testsOK = False
	    return
	self.config = ConfigParser.ConfigParser()
        self.config.read(configFile)
        self.currentDir = os.getcwd()
	self.testsOK = True
	self.testsSkipped = False 
        self.localVars = {}
	self.readConfig()
	if not self.testsOK:
	    return
	self.addLog("Log for AZOrange Tests",True)


    def __getFromConfig(self,section,option):
        #This method automatically pars any localVars present in the data got from config file
        strToParse = self.config.get(section,option)
        for var in self.localVars:
            strToParse = strToParse.replace("%"+var+"%",self.localVars[var])
        return strToParse
        
    def readConfig(self):
        section = "LocalSetupVars"
        if section in self.config.sections() and self.config.options(section):
             for option in self.config.options(section):    #Read all localVars defined if any
                self.localVars[option.upper()] = self.config.get(section,option).strip()

        section = "Tests"
        if not self.__validateOptions(section,["logfile","detailslogfile"]):
            return
        self.logFile = os.path.realpath(self.__getFromConfig(section,"logfile"))
	if not self.logFile:
	    self.logFile = os.path.join(self.currentDir,"tests.log")
        self.detailsLogFile = os.path.realpath(self.__getFromConfig(section,"detailslogFile"))
        if not self.detailsLogFile:
            self.detailsLogFile = os.path.join(self.currentDir,"testsDetails.log")
        self.testScript = self.__getFromConfig(section,"testscript")
        if not self.testScript:
            self.testScript = "runDailyTests.sh"

	section = "Paths"
	if not self.__validateOptions(section,["installdir"]):
            return
	self.trunkDir = os.path.realpath("../")#os.path.realpath(self.__getFromConfig(section,"trunkdir"))
        self.installDir = os.path.realpath(self.__getFromConfig(section,"installdir"))

        section = "Installation"
	if not self.__validateOptions(section,["installtype"]):
            return
        if not self.__validateOptions(section,["installtype"]):
            return
        self.installType = self.__getFromConfig(section,"installtype")
	if not "sys" in self.installType.lower():
	    self.installDir = self.trunkDir
	if not os.path.isdir(os.path.join(self.installDir,"azorange")) or not os.path.isdir(os.path.join(self.installDir,"tests")):
            print "ERROR: AZOrange not found on " + self.installDir 
            self.testsOK = False
            return
	else:
	    self.testsDir=os.path.join(self.installDir,"tests")
 
        section = "FeedBack"
        if not self.__validateOptions(section,["supportemails"]):
            return
        self.supportEmails = self.__getFromConfig(section,"supportemails")

    def __validateOptions(self,section,allOptions):
	if section not in self.config.sections():
	    print "Missing section "+section+" in setup file!"
            self.testsOK = False
            return False
	if sum([x in self.config.options(section) for x in allOptions])!=len(allOptions):
            print "Some options of section \""+section+"\" are missing:"+str([x for x in allOptions if x not in self.config.options(section)])
            self.testsOK = False
            return False
	return True

    def test(self):
        self.addLog("*Testing AZOrange")
        #run the defined AfterInstall tests
        os.chdir(self.testsDir)
        sys.path.append(self.testsDir)
	if not self.testScript:
            self.addLog("#No test Script defined in setup.ini ")
            self.testsOK = False
            return
        else:
            self.addLog("#Running the tests in Test suite")
        if not os.path.isfile("./"+self.testScript):
            self.addLog("#ERROR: " + self.testScript + " not found in the tests directory.")
            self.testsOK = False
            return
        status = commands.getstatusoutput("./" + self.testScript)
        
        lines = status[1].split("\n")
        lines.pop(0)
        testsLogFileText = ""
        testsLogDir = []
        try:
          for line in lines:
            if "LOGFILE:" in line and line[-1]!=":": 
	        testsLogFile = open(line.split(":")[1])
                testsLogDir.append(os.path.split(line.split(":")[1])[0])
                testsLogFileText += testsLogFile.read()
                testsLogFile.close()
        except:
          testsLogFileText = "Tests log file not found!"
        if testsLogFileText == "":
            testsLogFileText = "No LOGFILE returned from test script " + self.testScript
        for logDir in testsLogDir:
            if os.path.isdir(logDir):
                statusRM = commands.getstatusoutput("rm -rf "+logDir)
                if statusRM[0] != 0:
                     self.addLog("#WARNING: Could nor remove dir %s : %s" % (logDir,statusRM[1]))    
            else:
                self.addLog("#WARNING: Dir '%s' not found" % logDir)

        self.addLog((status[0], status[1] + "\n" + testsLogFileText))
        logStr = ""
        for line in lines:
            if "[34m" in line:
                #line = filter(lambda x: x in string.printable, line)
                #line = re.sub(r'[0-9]* %.',r'',re.sub(r'\[[0-9]*m',r'',line))
                if "SKIPPED" in line:
                    self.testsSkipped = True    
            #if "== Report ==" in line:
            #    break
                line = filter(lambda x: x in string.printable, line)
                line = re.sub(r'[0-9]* %.',r'',re.sub(r'\[[0-9]*m',r'',line))
                logStr+=line.strip()+"\n"
        self.addLog("#"+logStr)
	
    def printConfiguration(self,sendOnlyToDetails = False):
        logStr = ""
        for section in self.config.sections():
            logStr+="["+section+"]\n"
            for option in self.config.options(section):
                logStr+=" "+option+"="+self.__getFromConfig(section, option)+"\n"
        self.addLog("#Setup Configuration:")
        if sendOnlyToDetails:
            self.addLog((0,logStr))
        else:
            self.addLog("#"+logStr)
	
    def addLog(self,status,create=False):
	    if create:
                log = open(self.logFile,"w")
                detailsLog = open(self.detailsLogFile,"w")
		logStr = status + " (" +time.asctime() + ")\n" +"="*60
		log.write(logStr+"\n")
		detailsLog.write(logStr+"\n")
		print logStr
		log.close()
		detailsLog.close()
		return
	    #status can be the output of commands.getstatusoutput: (0, "something")
	    #If status is a string, it can be a new log section if started by "*" or a
	    #   simply line if started by "#" 
	    log = open(self.logFile,"a")
	    detailsLog = open(self.detailsLogFile,"a")
	    if type(status) in types.StringTypes:
		if status[0]=="#":
                    logStr = status[1:]
		elif status[0]=="*":
	            logStr = "-"*60 + "\n"+ status[1:] + " (" + time.asctime() + ")\n"+"-"*60
		else:
		    logStr ="-"*60 + "\n"+ status + " (" + time.asctime() + ")\n"+"-"*60
		log.write(logStr+"\n")
		detailsLog.write(logStr+"\n")
		print logStr

	    else:
		if status[0]==0:
	           if not status[1]:
		       log.close()
		       detailsLog.close()
		       return
		   else:
		       logStr = "Detail Code: OK" + "_"+str(time.time())
		       detailsLogStr = logStr+"\n\t\t" + status[1].replace("\n","\n\t\t")
		else:
		   self.testsOK = False
		   logStr = "ERROR: Detail Code: " + str(status[0]) + "_"+str(time.time())
		   detailsLogStr=logStr+"\n\t\t" + status[1].replace("\n","\n\t\t")
		log.write(logStr+"\n")
		detailsLog.write(detailsLogStr+"\n")
		print logStr
	    log.close()
	    detailsLog.close()
		
    def emailLog(self):
            self.addLog("*Support FeedBack")
            if not self.supportEmails:
                self.addLog("#No emails specified for feedBack")
                return
            allLogs=open(os.path.join(self.currentDir,"allLogs.log"),"w")
            log=open(self.logFile,"r")
            allLogs.write(log.read())
            allLogs.write("\n\n"+"*"*80+"\n\n\n")
            log.close()
            log=open(self.detailsLogFile,"r")
            allLogs.write(log.read())
            times = commands.getoutput("python " + os.path.join(test.currentDir,"chekLogTestTimes.py") + " " + os.path.join(test.currentDir,"testDetails.log"))
            allLogs.write("="*40)
            allLogs.write(times)
            log.close()
            allLogs.close()
            if not self.testsOK:
                subject = "AZOrange Tests report: FAILED"
            elif self.testsSkipped:
                subject = "AZOrange Tests report: SKIPPED!"
            else:
                subject = "AZOrange Tests report: Success!"
            subject += " (" + self.testScript.replace("run","").replace(".sh","")+")" 
            for email in self.supportEmails.split(","):
                status = commands.getstatusoutput('mail -s "%s" %s < %s' % (subject, email.strip(), os.path.join(self.currentDir,"allLogs.log")) )
                if status[0] != 0:
                    self.addLog(status)
                    self.addLog("#WARNING: Could not send a report to %s.\nPlease contact support." % email.strip())
                else:
                    self.addLog("#Report sent to "+email.strip())

 
if __name__ == "__main__":
    test = tester("./setup.ini")
    if not test.testsOK:
	print "Nothing Done!"
	sys.exit(1)

    #Sends the configuration used to the details only!
    test.printConfiguration(True)
    startTestsTime = time.time()
    #===== Test procedures=====
    test.test()
    #==========================
    test.addLog("*Finished")
    test.addLog("#The process spent " + str(int((time.time() - startTestsTime)/60)) + " minutes")

    if not test.testsOK:
        test.addLog("#ERRORS while running the tests!!")
    else:
        test.addLog("#Tests were all succeeded!!")

    #get back to the initial directory
    os.chdir(test.currentDir)
    #Send the log if required in setup file
    test.emailLog()
