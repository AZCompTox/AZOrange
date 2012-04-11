import ConfigParser
import commands
from optparse import OptionParser
import os
import string
import sys
import time
import types

SHELL_TYPE_BASH = "bash"
SHELL_TYPE_TCSH = "tcsh"

class Installer(object):

    def __init__(self, commandLine):
        self.AZOrangeInstallDir = None
        self.commandLineArguments = commandLine
        self.configFile = None
        self.currentDir = os.getcwd()
        self.localVars = {}
        self.prepareOnly = False
        self.successInstall = True
        self.verbosedLogging = False


    def main(self):
        options, args = self.__parse(self.commandLineArguments)
        if not os.path.isfile(options.file):
            print "Missing the configuration file. "+options.file+" not found.\n You may use the setupTemplate.ini file to create your setup.ini file."
            self.successInstall = False
            return
        self.configFile = options.file
        self.prepareOnly = options.prepare
        self.verbosedLogging = options.verbose
        self.readConfig()
        self.addLog("Log for AZOrange installation",True)

        os.system("clear")
        
        if not self.successInstall:
            print "Nothing Done!"
            sys.exit(1)
    
        #Sends the configuration used to the details only!
        self.printConfiguration(True)
        startInstallTime = time.time()
            
        # Check for some important requirements
        self.addLog("*Requirements and System Info")
        st,out = commands.getstatusoutput('python -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()"; uname -a; lsb_release -a')
        self.addLog('#'+out) 
        if not os.path.isfile("/bin/tcsh"):
            self.addLog("#/bin/tcsh is missing. tcsh is required for azorange to work properly.")
            self.successInstall = False
        #=====install procedures=====
        if self.successInstall:
            self.prepareInstallDirs()
    
        #Checkout any 3rd party software to the proper locations
        if self.successInstall:
            self.checkoutFTM()
        if self.successInstall:
            self.checkoutStructClust()
        if self.successInstall:
            self.checkoutFMINER()
        if self.successInstall:
            self.checkOutOpenAZO()
        if self.successInstall:
            self.checkOutCinfony()
        if self.successInstall:
            self.checkOutRdkit()
        if self.successInstall:
            self.checkOutCDK()
        if self.successInstall:
            self.checkOutOrange() 
        if self.successInstall:
            self.checkOutAppspack()
        if self.successInstall:
            self.checkOutOpencv()
        if self.successInstall:
            self.checkOutBoost()
        if self.successInstall:
            self.checkOutPLearn()
        if self.successInstall:
            self.checkOutOasa()
    
        if self.successInstall: 
            if self.prepareOnly:
                self.addLog("*Everything is now unpacked and in place ready for install!")
            else:
                self.install()
    
        if not self.successInstall:
            self.addLog("*ERROR: Installation aborted!")
    
        if self.successInstall and not self.prepareOnly:
            self.createProfileExample(SHELL_TYPE_TCSH)
            self.createProfileExample(SHELL_TYPE_BASH)
            self.InstallCacheCleaner()
            self.runAfterInstallScripts()
    
        #==========================
        self.addLog("*Finished")
        self.addLog("#The process spent " + str(int((time.time() - startInstallTime)/60)) + " minutes")
    
        #get back to the initial directory
        os.chdir(self.currentDir)
    
        if not self.successInstall:
            self.addLog("#ERRORS during installation!!")
            self.emailLog()
            sys.exit(1)
        elif not self.prepareOnly :
            startAPPTemplate = os.path.join(self.trunkDir,'install/startAZOrange')
            startAPPTarget = os.path.join(self.AZOrangeInstallDir,'startAZOrange')
            AppLaucherTemplate = os.path.join(self.trunkDir,'install/AZOrange.desktop')
            AppLaucherTarget = None
            if "HOME" in os.environ:
                if os.path.isdir(os.path.join(os.getenv("HOME"),'Desktop')):
                    AppLaucherTarget = os.path.join(os.getenv("HOME"),'Desktop')
            #Create the Application GUI Starter script
            cmd = 'sed "s|AZO_INSTALL_DIR|'+self.AZOrangeInstallDir+'|g" ' + startAPPTemplate + ' > '+startAPPTarget
            self.__logAndExecute(cmd)
            self.__logAndExecute("chmod a+x " +startAPPTarget)
            #Create a Launcher in the Desktop
            thisOS = commands.getstatusoutput("uname -o")
            if "gnu/linux" in thisOS[1].lower() and AppLaucherTarget:
                cmd='sed "s|AZO_INSTALL_DIR|'+self.AZOrangeInstallDir+'|g" '+ AppLaucherTemplate + ' > ' + os.path.join(AppLaucherTarget,'AZOrange.desktop')
                self.__logAndExecute(cmd)
                self.__logAndExecute("chmod a+x " + os.path.join(AppLaucherTarget,'AZOrange.desktop'))
            elif AppLaucherTarget:
                self.__logAndExecute("ln -s " + startAPPTarget + " " + os.path.join(AppLaucherTarget,'AZOrange'))
                self.__logAndExecute("chmod a+x " + os.path.join(AppLaucherTarget,'AZOrange'))
            #associate Orange canvas files with AZOrange Application: Not critial if fails
            try:
                home = os.environ["HOME"]
                mimeStr = """
[Default Applications]
application/xml=AZOrange.desktop

[Added Associations]
application/xml=AZOrange.desktop;
                """
                appDir =  os.path.join(home,".local/share/applications")
                if not os.path.isfile(os.path.join(appDir,"mimeapps.list")):
                    os.system("mkdir -p "+appDir)
                    fileh=open(os.path.join(appDir,"mimeapps.list"),"w")
                    fileh.write(mimeStr)
                    fileh.close()
                else:
                    fileh=open(os.path.join(appDir,"mimeapps.list"),"r")
                    lines = fileh.readlines()
                    fileh.close()
                    fileh=open(os.path.join(appDir,"mimeapps.list"),"w")
                    for line in lines:
                        if "application/xml" in line:
                            fileh.write("application/xml=AZOrange.desktop\n")
                        else:
                            fileh.write(line)
                    fileh.close()
                cmd='sed "s|AZO_INSTALL_DIR|'+self.AZOrangeInstallDir+'|g" '+ AppLaucherTemplate + ' > ' + os.path.join(appDir,'AZOrange.desktop')
                self.__logAndExecute(cmd)
                self.__logAndExecute("chmod a+x " + os.path.join(appDir,'AZOrange.desktop'))
            except:
                print "Could not associate .ows files to AZOrange"
            self.addLog("#Installation done successfully!")
        else:
            self.addLog("#Preparation done successfully!")
    
        #Send the log if required in setup file
        self.emailLog()
    
        #Start the tests if required
        if self.runTests:
            os.system(os.path.join(self.currentDir,"runTests"))
            status,arch = commands.getstatusoutput('uname -i')
            if  '_64' not in arch:
                print "AZOrange tests were done using 64 bits architecture."
                print "In 32 bits is expected that some tests fail although does not necessarily means that there is any problem on running AZOrange."

    # Needed for modules to be loaded.
    def module(self,command, *arguments):
        commands = os.popen('%s/bin/modulecmd python %s %s' % (os.environ['MODULESHOME'], command, string.join(arguments))).read()
        exec commands


    def __logAndExecute(self, command):
        if self.verbosedLogging:
            self.addLog("#install - About to execute (in " + os.getcwd() + "): " + command)
        status, output = commands.getstatusoutput(command)
        self.addLog((status, "#install - Output from command: " + str(output)))
        return status, output


    def __update_EnvVars(self, envFile):
        """ This function will add to self.EnvVars the variables in the 'envFile'.
            For the variables that already exists, they will prepend the values in the file  to the values already 
            defined, if they do not exist already.
        """
        envVars = {}
        if os.path.isfile(envFile):
            file = open(envFile,"r")
            lines = file.readlines()
            file.close()
            for line in lines:
                envVars[line.split("=")[0].strip()] = line.split("=")[1].strip()
        for varName in envVars:
            values = [var.strip() for var in envVars[varName].split(":")]
            if varName not in self.EnvVars:
                self.EnvVars[varName] = values
            else:
                for value in values:
                    if value not in self.EnvVars[varName]:
                        self.EnvVars[varName].insert(0,value)


    def __getFromConfig(self,section,option):
        #This method automatically pars any localVars present in the data got from config file
        strToParse = self.config.get(section,option)
        for var in self.localVars:
            strToParse = strToParse.replace("%"+var+"%",self.localVars[var])
        return strToParse


    def readConfig(self):
        self.config = ConfigParser.ConfigParser()
        self.config.read(self.configFile)

        section = "LocalSetupVars"
        self.localVars={}
        if section in self.config.sections() and self.config.options(section):
             for option in self.config.options(section):    #Read all localVars defined if any
                self.localVars[option.upper()] = self.config.get(section,option).strip()

        section = "Installation"
        if not self.__validateOptions(section,["logfile","templateprofilefile","detailslogfile"]):
            return
        self.logFile = os.path.realpath(self.__getFromConfig(section,"logfile"))
        if not self.logFile:
            self.logFile = os.path.join(self.currentDir,"install.log")
        self.detailsLogFile = os.path.realpath(self.__getFromConfig(section,"detailslogFile"))
        if not self.detailsLogFile:
            self.detailsLogFile = os.path.join(self.currentDir,"details.log")
        self.TemplateProfileFile = os.path.realpath(self.__getFromConfig(section,"templateprofilefile"))
        if not self.TemplateProfileFile:
            self.TemplateProfileFile = os.path.join(self.currentDir,"templateProfile")
        self.TemplateProfileFileBash = os.path.realpath(self.__getFromConfig(section,"templateprofilefilebash"))
        if not self.TemplateProfileFileBash:
            self.TemplateProfileFileBash = os.path.join(self.currentDir,"templateProfileBash")

        section = "Paths"
        if not self.__validateOptions(section,["builddir","installdir"]):
            return
        self.trunkDir = os.path.realpath("../")#os.path.realpath(self.__getFromConfig(section,"trunkdir"))
        self.DepSrcDir = os.path.join(self.trunkDir,"orangeDependencies","src")
        self.buildDir = os.path.realpath(self.__getFromConfig(section,"builddir"))
        self.installDir = os.path.realpath(self.__getFromConfig(section,"installdir"))

        section = "Repository"
        if not self.__validateOptions(section,["repository","repointeraction"]):
            return
        self.repo = self.__getFromConfig(section,"repository")
        self.repoInter= self.__getFromConfig(section,"repointeraction")
        if "no" in self.repoInter.lower():
            self.repoInter = "no"
        elif "yes" in self.repoInter.lower():
            self.repoInter = "yes"
        else:
            self.repoInter = None
        if not self.repoInter:
            print "Invalid repo interaction: ",self.repoInter
            print "   Use 'yes'  or 'no' "
            self.successInstall = False
            return

        section = "Tests"
        if not self.__validateOptions(section,["runtestsafterinstall"]):
            return
        self.runTests = "yes" in self.__getFromConfig(section,"runtestsafterinstall").lower()

        section = "Installation"
        if not self.__validateOptions(section,["installtype","openinstallation","cleaninstallcache","precleanbuilddir","precleaninstalldir"]):
            return
        self.installType = self.__getFromConfig(section,"installtype")
        self.openInstallation = self.__getFromConfig(section,"openinstallation")
        if "false" in self.openInstallation.lower():
            self.openInstallation = False
        else:
            self.openInstallation = True
        if "sys" in self.installType.lower():
            self.installType = "system"
            self.AZOrangeInstallDir = self.installDir
        else:
            self.installType = "developer"
            self.AZOrangeInstallDir = self.trunkDir
            if self.repoInter == "export":
                print "For developer installation you cannot use 'export' interaction with repository"
                self.successInstall = False
                return
        self.cleanInstallCache= self.__getFromConfig(section,"cleaninstallcache")
        if "yes" in self.cleanInstallCache.lower():
            self.cleanInstallCache = True
        else:
            self.cleanInstallCache = False
        self.preCleanTrunkDir = False #"yes" in self.__getFromConfig(section,"precleantrunkdir").lower()
        self.preCleanBuildDir = "yes" in self.__getFromConfig(section,"precleanbuilddir").lower()
        self.preCleanInstallDir = "yes" in self.__getFromConfig(section,"precleaninstalldir").lower()

        section = "FeedBack"
        if not self.__validateOptions(section,["supportemails"]):
            return
        self.supportEmails = self.__getFromConfig(section,"supportemails")

        section = "Advanced"
        if not self.__validateOptions(section,["platformtype"]):
            return
        self.platformType = self.__getFromConfig(section,"platformtype")
        #Set the InstallTimeModules to be the union of the modules and InstallTimeModules. 
        #only the modules var will be used to compose the templateProfile
        #only the InstallTimeModules var will be sent to the setup script
        self.modules = []
        if "modules" in self.config.options(section):
            self.modules = self.__getFromConfig(section,"modules")
            if self.modules and self.modules.lower()!="none" and self.modules.lower()!="no":
                self.modules = [module.strip() for module in self.modules.split(",")]
            else:
                self.modules = []
        self.InstallTimeModules = []
        if "installtimemodules" in self.config.options(section):
            self.InstallTimeModules = self.__getFromConfig(section,"installtimemodules")
            if self.InstallTimeModules and self.InstallTimeModules.lower()!="none" and self.InstallTimeModules.lower()!="no":
                self.InstallTimeModules = [module.strip() for module in self.InstallTimeModules.split(",")]
            else:   
                self.InstallTimeModules = []
        #Compose the InstallTimeModules and the modules vars properly
        self.InstallTimeModules = self.InstallTimeModules + self.modules
        if len(self.InstallTimeModules)>0:
            self.InstallTimeModules = str(self.InstallTimeModules)
        else:
            self.InstallTimeModules = "None"
        if len(self.modules)>0:
            self.modules = str(self.modules)                
        else:
            self.modules = "None"

        self.RunTimeModules = "None"
        if "runtimemodules" in self.config.options(section):
            self.RunTimeModules = self.__getFromConfig(section,"runtimemodules")
            if self.RunTimeModules and self.RunTimeModules.lower()!="none" and self.RunTimeModules.lower()!="no":
                self.RunTimeModules = str([module.strip() for module in self.RunTimeModules.split(",")])
            else:
                self.RunTimeModules = "None"

        self.preInstallModules = None
        if "preinstallmodules" in self.config.options(section):
            self.preInstallModules = self.__getFromConfig(section,"preinstallmodules")
            if self.preInstallModules and self.preInstallModules.lower()!="none" and self.preInstallModules.lower()!="no":
                self.preInstallModules = [module.strip() for module in self.preInstallModules.split(",")]
            else:
                self.preInstallModules = None
            if self.preInstallModules:
                for moduleX in self.preInstallModules:
                    self.module("load", moduleX)

        self.GCCmoduleForAppspackMPI = None
        if "gccmoduleforappspackmpi" in self.config.options(section):
            self.GCCmoduleForAppspackMPI = self.__getFromConfig(section,"gccmoduleforappspackmpi").strip()
            if not self.GCCmoduleForAppspackMPI or self.GCCmoduleForAppspackMPI.lower()=="none" or self.GCCmoduleForAppspackMPI.lower=="no":
                self.GCCmoduleForAppspackMPI = None
        self.sources = None
        if "sources" in self.config.options(section):
            self.sources = self.__getFromConfig(section,"sources")
            if self.sources and self.sources.lower()!="none" and self.sources.lower()!="no":
                self.sources = [source.strip() for source in self.sources.split(",")]
            else:
                self.sources = None

        section = "EnvVars"
        self.EnvVars={}
        if section in self.config.sections() and self.config.options(section):
            for option in self.config.options(section):    #Read all envVard defined in the EnvVars section if any
                configValues = self.__getFromConfig(section,option).split(":")   #For each envVar get a list of elements
                newValues = []
                if not configValues[0]:         #If the envVar is empty, set it yo ""
                    os.environ[option.upper()] = ""
                    self.EnvVars[option.upper()] = []
                    continue 
                for value in configValues:    
                    if "$" in value[0]:       #Load other envVars to this one if defined with $
                        newValues += os.environ[value[1:]].split(":")
                    else:
                        newValues.append(value)
                #Compose the string to assign to the envVar
                strValues = "" 
                for idx,value in enumerate(newValues):
                    if idx: strValues += ":"
                    strValues += value
                    #Assure that if something was added to PYTHONPATH, that is also in the sys.path
                    if option.upper()=="PYTHONPATH":   
                        sys.path.insert(0,value)
                os.environ[option.upper()] = strValues
                self.EnvVars[option.upper()] = newValues 

        section = "Dependencies"
        self.dependencies = {}
        if section in self.config.sections() and self.config.options(section):
            for option in self.config.options(section):    #Read all dependencies if any
                self.dependencies[option] = self.__getFromConfig(section,option)


    def __validateOptions(self,section,allOptions):
        if section not in self.config.sections():
            print "Missing section "+section+" in setup file!"
            self.successInstall = False
            return False
        if sum([x in self.config.options(section) for x in allOptions])!=len(allOptions):
            print "Some options of section "+section+" are missing!"
            self.successInstall = False
            return False
        return True


    def prepareInstallDirs(self):        
        #trunkDir must already exist if repoInter is "" or "update"
        #if repoInter is checkout, clean if exists the trunkDir
        #Always clean the installDir
        #if installType is developer,the installDir will not be created and 
        #  BuildDir is not deleted
        #if installType is system, we MUST create installDir
        self.addLog("*Preparing Directories")
        if (not self.trunkDir and self.repoInter != "export") or not self.buildDir:
            self.addLog("#ERROR: Missing the definition of some Paths in Setup file.")
            self.successInstall = False
            return

        #Prepare the Trunk Dir
        if (self.repoInter == "yes") and (not os.path.isdir(self.trunkDir)):
            self.addLog("#ERROR: TrunkDir must exist for the chosen RepositoryInteraction.")
            self.successInstall = False
            return
        elif self.repoInter == "checkout":
            if os.path.isdir(self.trunkDir) and self.preCleanTrunkDir:
                self.addLog("#Removing existing TrunkDir")
                self.__logAndExecute("chmod -R 777 " + self.trunkDir)
                self.__logAndExecute("rm -Rf " + self.trunkDir)
            if not os.path.isdir(self.trunkDir): 
                self.addLog("#Creating TrunkDir")
                self.__logAndExecute("mkdir -p " + self.trunkDir)
        #elif self.repoInter!="export":  #is "" or "update"
        #    if not os.path.isdir(os.path.join(self.trunkDir,".svn")):
        #        self.addLog("#WARNING: Seems that the trunkDir is not a subversion trunk!")

        #Prepare BuildDir
        if os.path.isdir(self.buildDir) and (self.preCleanBuildDir):
            self.addLog("#Removing existing buildDir")
            self.__logAndExecute("chmod -R 777 " + self.buildDir)
            self.__logAndExecute("rm -Rf " + self.buildDir)
        if not os.path.isdir(self.buildDir):
            self.addLog("#Creating buildDir")
            self.__logAndExecute("mkdir -p " + self.buildDir)

        #Prepare InstallDir
        #Create the installDir if is a system installation. If is a dev install, we will 
        #  install back to trunk, so we do not need installDir
        if self.installType == "system":
            if self.trunkDir+"/" in self.installDir:
                self.addLog("#ERROR: Invalid installDir defined in setup file. The system installation needs a different install dir.")
                self.successInstall = False
                return
            if os.path.isdir(self.installDir) and self.preCleanInstallDir:
                self.addLog("#Removing existing installDir")
                self.__logAndExecute("chmod -R 777 " + self.installDir)
                self.__logAndExecute("rm -Rf " + self.installDir)
            if not os.path.isdir(self.installDir):
                self.addLog("#Creating installDir")
                self.__logAndExecute("mkdir -p " + self.installDir)
        if self.openInstallation:
            self.__logAndExecute("mkdir -p " + self.DepSrcDir)


    def checkOutOpenAZO(self):
        USE_INSTALLED = False
        if self.repoInter=="yes":
            #Update the AZO source from GITHub
            os.chdir(self.trunkDir)
            if self.openInstallation:
                self.addLog("*OpenAZOrange: Using current files.")
                self.addLog("#trunk: "+self.trunkDir)
                #self.__logAndExecute("git pull")
            else:
                if self.installType == "developer":
                        if "openazo" in self.dependencies:
                            depCfg = self.dependencies["openazo"].split(",")
                            if len(depCfg)>=3 and depCfg[2] == "*":
                                USE_INSTALLED = True
                            else:
                                USE_INSTALLED = False

                        if USE_INSTALLED: 
                           self.addLog("*Not cloning openAZO from GIT")
                           return

                        self.addLog("*Cloning from GIT")
                        self.addLog("#trunk: "+self.trunkDir)
                        self.__logAndExecute("rm -rf AZOrange")
                        self.__logAndExecute("git clone "+ self.repo)
                        self.__logAndExecute("cp -Rf AZOrange/* .")
                        self.__logAndExecute("cp -Rf AZOrange/.git " + self.trunkDir)
                        self.__logAndExecute("rm -rf AZOrange")
                else:
                        if "openazo" not in self.dependencies:
                            self.addLog((1,"Missing OpenAZO definition in Dependencies at setup.ini"))  
                        else:
                            depCfg = self.dependencies["openazo"].split(",")
                            URL = os.path.join(self.DepSrcDir,depCfg[0])
                            if len(depCfg)>=3 and depCfg[2] == "*":
                                USE_INSTALLED = True
                            else:
                                USE_INSTALLED = False

                        if not URL or USE_INSTALLED or self.repoInter == "no":
                           self.addLog("*Not unpacking openAZO")
                           return

                        self.addLog("*Extracting openAZO from " + URL)
                        self.addLog("#trunk: "+self.trunkDir)
                        os.chdir(self.DepSrcDir)
                        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"AZOrange"))
                        self.__logAndExecute("tar xfz " + os.path.split(URL)[-1])
                        self.__logAndExecute("cp -Rf AZOrange/* " + self.trunkDir)
                        self.__logAndExecute("cp -Rf AZOrange/.git " + self.trunkDir)
                        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"AZOrange"))
                #Revert any files replaced from the Git repo. We want to assure that the AZO inside files prevail
                os.chdir(self.trunkDir)
                self.addLog("#Reverting full RepoDir to SVN version...")
                self.__logAndExecute("svn revert -R ./*") 
                os.chdir(self.currentDir)


    def checkOutCDK(self):
        # Get the dependency Config
        name = "cdk"
        if name not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies[name].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading "+name)
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,name))
        self.__logAndExecute("mkdir " +  os.path.join(self.DepSrcDir,name))
        os.chdir(os.path.join(self.DepSrcDir,"cdk"))
        jarFile = os.path.split(URL)[-1].strip()
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + jarFile)
                self.__logAndExecute("wget " + URL )
        else:
                self.addLog("*Using "+name+" in SVN Repo (Not implemented)")

    def checkoutFMINER(self):
        # Get the dependency Config
        name = "fminer"
        if name not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies[name].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading "+name)
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,name))
        os.chdir(self.DepSrcDir)
        tarFile = os.path.split(URL)[-1].strip()
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("git clone " + URL + " " + os.path.join(self.DepSrcDir,name))
        else:
                self.addLog("*Using "+name+" in SVN Repo (Not implemented)")
                return
        if not os.path.isdir(os.path.join(self.DepSrcDir,name,"libbbrc")):
            self.addLog("#ERROR: Could not fet fminer source code.")
            self.successInstall = False
            return


    def checkoutStructClust(self):
        # Get the dependency Config
        name = "clustering"
        if name not in self.dependencies:
            self.addLog("Name " + str(name) + " not in dependencies")
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies[name].split(",")
            URL = depCfg[0]
            self.addLog(URL)
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading "+name)
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,name))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                tarFile = "structuralClustering.tar.gz"
                dwnldFile = os.path.split(URL)[-1].strip()
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + tarFile)
                self.__logAndExecute("rm -rf " + dwnldFile)
                self.__logAndExecute("wget " + URL )
                self.__logAndExecute("mv "+dwnldFile+" "+tarFile)
        else:
                self.addLog("*Using "+name+" in SVN Repo")
                tarFile = URL 
        UnpackCmd = "tar "
        if  tarFile[-6:] == "tar.gz":
            UnpackCmd += "xfz "
        elif tarFile[-6:] == "tar.bz2":
            UnpackCmd += "xfj "
        else:
            self.addLog("#ERROR: Not a known tar file.")
            self.successInstall = False
            return
        self.__logAndExecute(UnpackCmd + tarFile)

    def checkoutFTM(self):
        # Get the dependency Config
        name = "ftm"
        if name not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies[name].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading "+name)
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,name))
        os.chdir(self.DepSrcDir)
        tarFile = "ftm.tar.gz"
        dwnldFile = os.path.split(URL)[-1].strip()
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + tarFile)
                self.__logAndExecute("rm -rf " + dwnldFile)
                self.__logAndExecute("wget " + URL )
                self.__logAndExecute("mv "+dwnldFile+" "+tarFile)
        else:
                self.addLog("*Using "+name+" in SVN Repo (Not implemented yet)")
                return
        UnpackCmd = "tar "
        if  tarFile[-6:] == "tar.gz":
            UnpackCmd += "xfz "
        elif tarFile[-6:] == "tar.bz2":
            UnpackCmd += "xfj "
        else:
            self.addLog("#ERROR: Not a known tar file.")
            self.successInstall = False
            return
        self.__logAndExecute(UnpackCmd + tarFile)


    def checkOutRdkit(self):
        # Get the dependency Config
        name = "rdkit"
        if name not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies[name].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading "+name)
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,name))
        os.chdir(self.DepSrcDir)
        tarFile = os.path.split(URL)[-1].strip()
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + tarFile)
                self.__logAndExecute("wget " + URL )
        else:
                self.addLog("*Using "+name+" in SVN Repo (Not implemented)")
                return
        UnpackCmd = "tar "
        if  tarFile[-6:] == "tar.gz":
            UnpackCmd += "xfz "
            unpackDir = tarFile[0:tarFile.rfind(".tar")]
        elif  tarFile[-4:] == ".tgz":
            UnpackCmd += "xfz "
            unpackDir = tarFile[0:tarFile.rfind(".tgz")]
        elif tarFile[-6:] == "tar.bz2":
            UnpackCmd += "xfj "
            unpackDir = tarFile[0:tarFile.rfind(".tar")]
        else:
            self.addLog("#ERROR: Not a known tar file.")
            self.successInstall = False
            return 
        self.__logAndExecute(UnpackCmd + tarFile)
        self.__logAndExecute("mv " + unpackDir + " " + name )


    def checkOutCinfony(self):
        # Get the dependency Config
        name = "cinfony"
        if name not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies[name].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading "+name)
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,name))
        os.chdir(self.DepSrcDir)
        tarFile = "cinfony_Download.tar.gz"
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + tarFile)
                # Download the File 
                self.__logAndExecute("wget " + URL + " -O " + tarFile)
        else:
                self.addLog("*Using "+name+" in SVN Repo (Not implemented)")
        UnpackCmd = "tar "
        if  tarFile[-6:] == "tar.gz":
            UnpackCmd += "xfz "
        elif tarFile[-6:] == "tar.bz2":
            UnpackCmd += "xfj "
        else:
            self.addLog("#ERROR: Not a known tar file.")
            self.successInstall = False
            return 
        # This file has different folder name from the tar file
        os.mkdir("./"+name)
        UnpackCmd += " " + tarFile + " -C "+name
        self.__logAndExecute(UnpackCmd)
        folderName = os.listdir("./"+name)[0]
        self.__logAndExecute("mv "+os.path.join(name,folderName,"*") + " " + name)
        self.__logAndExecute("rmdir "+os.path.join(name,folderName))
        #self.__logAndExecute("mv " + tarFile[0:tarFile.rfind(".tar")] + " " + name )


    def checkOutOrange(self):
        # Get the dependency Config
        if "orange" not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies["orange"].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading/unpacking orange")
           return 

        #self.__logAndExecute("rm -rf " + os.path.join(self.trunkDir,"orange/*"))
        #This command may have some failures, but it's no problem. We just want to delete if there is something to delete!
        self.__logAndExecute("mkdir -p " + os.path.join(self.trunkDir,"orange"))
        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"orange"))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Checking out from orange bitbucket to trunk ("+URL+":"+REV+")")
                #os.chdir(os.path.join(self.trunkDir,"orange"))
                self.__logAndExecute("hg clone " + URL + " ./orange")
                os.chdir("orange")
                self.__logAndExecute("hg update " + REV)
        else:
                self.addLog("*Extracting orange from " + URL)
                self.__logAndExecute("tar xfz " + os.path.split(URL)[-1])
        os.chdir(self.DepSrcDir)

        # Apply Patch
        self.addLog("#Applying Patch...")
        os.chdir(os.path.join(self.currentDir,"Patches"))
        status,out = commands.getstatusoutput("./applyOrangePatch.sh %s" % (os.path.join(self.DepSrcDir,"orange")))
        if status != 0:
            self.addLog("#WARNING: Patch was not properly applied. It might be that it has been already applied.")
            self.addLog("#         Please check the next details.")
            self.addLog([0,out])
        else:
            self.addLog([status,out])

        os.chdir(self.currentDir)


    def checkOutAppspack(self):
        # Get the dependency Config
        if "appspack" not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies["appspack"].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading AppsPack")
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"appspack"))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Downloading AppsPack to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + os.path.split(URL)[-1])
                self.__logAndExecute("wget " + URL )
        else:   
                self.addLog("*Using APPSPACK in SVN Repo")
        self.__logAndExecute("tar xfz " + os.path.split(URL)[-1])
        self.__logAndExecute("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " appspack" )
        os.chdir(self.currentDir)


    def checkOutOpencv(self):
        # Get the dependency Config
        if "opencv" not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies["opencv"].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading opencv")
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"opencv"))

        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Downloading opencv to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + os.path.split(URL)[-1])
                self.__logAndExecute("wget " + URL )
        else:
                self.addLog("*Using opencv in SVN Repo")
        self.__logAndExecute("tar xfj " + os.path.split(URL)[-1])
        self.__logAndExecute("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " opencv" )

        # Apply Patch
        self.addLog("#Applying Patch...")
        os.chdir(os.path.join(self.currentDir,"Patches"))
        status,out = commands.getstatusoutput("./applyOpenCVPatch.sh %s" % (os.path.join(self.DepSrcDir,"opencv")))
        if status != 0:
            self.addLog("#WARNING: Patch was not properly applied. It might be that it has been already applied.")
            self.addLog("#         Please check the next details.")
            self.addLog([0,out])
        else:
            self.addLog([status,out])
        os.chdir(self.currentDir)


    def checkOutOasa(self):
        #Get the oasa dependency config
        if "oasa" not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies["oasa"].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading Oasa")
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"oasa"))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Oasa is not needed in this installation")
                return
                # To download and install oasa in the open installation, remove both 2 previous lines and uncomment next 3 lines
                #self.addLog("*Downloading Oasa to trunk ("+URL+":"+REV+")")
                #self.__logAndExecute("rm -rf " + os.path.split(URL)[-1])
                #self.__logAndExecute("wget " + URL )
        else:
                self.addLog("*Using oasa in SVN Repo")
        self.addLog("#Extracting oasa from " + URL )
        self.__logAndExecute("tar xfz " + os.path.split(URL)[-1])
        self.__logAndExecute("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " oasa" )
        os.chdir(self.currentDir)


    def checkOutBoost(self):
        # Get the dependency Config
        if "boost" not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies["boost"].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading Boost")
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"boost"))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Downloading Boost to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("rm -rf " + os.path.split(URL)[-1])
                self.__logAndExecute("wget " + URL )
        else:
                self.addLog("*Using boost in SVN Repo")
        self.__logAndExecute("tar xfz " + os.path.split(URL)[-1])
        self.__logAndExecute("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " boost" )

        os.chdir(self.currentDir)


    def checkOutPLearn(self):
        # Get the dependency Config
        if "plearn" not in self.dependencies:
            URL = None
            REV = None
            USE_INSTALLED = True
        else:
            depCfg = self.dependencies["plearn"].split(",")
            URL = depCfg[0]
            if len(depCfg)<2 or depCfg[1] == "":
                REV = "HEAD"
            else:
                REV = depCfg[1]
            if len(depCfg)>=3 and depCfg[2] == "*":
                USE_INSTALLED = True
            else:
                USE_INSTALLED = False

        if not URL or USE_INSTALLED or self.repoInter == "no":
           self.addLog("*Not downloading plearn")
           return

        self.__logAndExecute("rm -rf " + os.path.join(self.DepSrcDir,"plearn"))
        if self.openInstallation:
                self.addLog("*Checking out from plearn SVN to trunk ("+URL+":"+REV+")")
                self.__logAndExecute("mkdir -p " + os.path.join(self.DepSrcDir,"plearn"))
                os.chdir(os.path.join(self.DepSrcDir,"plearn"))
                self.__logAndExecute("svn export" + " -r " + REV + " " + URL + " ./ --force")
        else:
                self.addLog("*Extracting PLearn from " + URL)
                os.chdir(self.DepSrcDir)
                self.__logAndExecute("tar xfz " + os.path.split(URL)[-1])
        self.__logAndExecute('find '+os.path.join(self.DepSrcDir,"plearn")+' -name ".svn" | xargs rm -Rf')

        # Apply Patch
        self.addLog("*Applying Patch...")
        os.chdir(os.path.join(self.currentDir,"Patches"))
        status,out = commands.getstatusoutput("./applyPLearnPatch.sh %s" % (os.path.join(self.DepSrcDir,"plearn")))
        if status != 0:
            self.addLog("#WARNING: Patch was not properly applied. It might be that it has been already applied.")
            self.addLog("#         Please check the next details.")
            self.addLog([0,out])
        else:
            self.addLog([status,out])
        os.chdir(self.currentDir)


    def install(self):
        # Apply Patch
        if not self.openInstallation:
            self.addLog("*Applying AZInHouse Patch...")
            os.chdir(os.path.join(self.currentDir,"Patches"))
            status,out = commands.getstatusoutput("./applyAZInHousePatch.sh %s" % (self.trunkDir))
            if status != 0:
                self.addLog("#WARNING: Patch was not properly applyed. It might be that it has been already applyed.")
                self.addLog("#         Please check the next details.")
                self.addLog([0,out])
            else:
                self.addLog([status,out])

        self.addLog("*Compiling and Installing")
        os.chdir(self.buildDir)
        self.addLog("#Copying from Trunk to build dir")
        self.addLog("#"+self.trunkDir+" -> "+self.buildDir)

        #Export from the Trunk working copy
        self.__logAndExecute("cp -Rf " + os.path.join(self.trunkDir,"*") + " ./ ")
        # Remove all .svn and .git folders
        self.__logAndExecute('find ./ -name ".svn" | xargs rm -Rf')
        self.__logAndExecute('find ./ -name ".git" | xargs rm -Rf')

        #The orange core setup.py that is placed in orange/install-scripts/linux must be copied to orange
        #self.__logAndExecute("/bin/cp -f " + os.path.join(self.buildDir,"orange/install-scripts/linux/setup.py") + " " + os.path.join(self.buildDir,"orange"))

        #In buildDir> python setup.py build ...
        self.addLog("#Compiling by running setup.py")
        self.addLog("#Installation will be to made to " + self.AZOrangeInstallDir)
        envTmp = os.path.join(self.currentDir,"env.tmp")
        if os.path.isfile(envTmp):
            os.system("/bin/rm -f " + envTmp)
        cmdFormat = 'python setup.py --action build --platform %s --builddir %s --installdir %s --modulestoload "%s"' + \
                             ' --envfile %s --dependencies "%s" --appspackgccmodule "%s" --openinstall %s' + \
                             ' --logfile "%s" --detailslogfile "%s" --verbose %s'
        self.__logAndExecute( cmdFormat % (self.platformType, self.buildDir, self.AZOrangeInstallDir, self.InstallTimeModules,
                                           envTmp , self.dependencies, self.GCCmoduleForAppspackMPI, self.openInstallation,
                                           str(self.logFile), str(self.detailsLogFile), str(self.verbosedLogging)))
        #Update the local self.EnvVars with the ones passed by setup.py which will be needed at runtime
        self.__update_EnvVars(envTmp) 


    def runAfterInstallScripts(self):        
        # This runs  once the .../azorange/bin/archive.sh in order to create the AZO_NFZ_scratchDir if it does not exists
        self.addLog("#Running after-install scripts")
        self.__logAndExecute(os.path.join(self.AZOrangeInstallDir,'azorange/bin/clean.sh'))


    def createProfileExample(self, shellType = SHELL_TYPE_TCSH):
        # In addition to the place defined in setup.ini for the template profile, a file called templateProfile will always be
        # placed in:
        #    - Install Dir ($AZORANGEHOME):
            #       - system installation: place Profile Template at root of installDir
            #       - developer installation : place Profile template at root of trunkDir
        #    - within this running install dir
        # NOTE: The $AZORANGEHOME/templateProfile will be used by MPI calls. Do not remove it from there and
        #       do not rename it.

        self.addLog("*Create a Template  profile")
        #Compose the template profile content
        if shellType == SHELL_TYPE_BASH:
            strFile =  "#!/bin/bash\n"
        else:
            strFile =  "#!/bin/tcsh\n"
        localVars = {"installDir":self.AZOrangeInstallDir}
        strFile +=  "# Template profile for azorange installation at %(installDir)s\n" % localVars
        if self.sources:
            strFile += "\n# Additional sources needed\n"
            for source in self.sources:
                if shellType == SHELL_TYPE_BASH:
                    if os.path.isfile(source):
                        strFile += ". "+ source  + "\n"
                    else:
                        strFile += "#. "+ source  + "\n"
                        self.addLog("#The specified source file does not exist: " + source)
                else:
                    if os.path.isfile(source):
                        strFile += "source "+ source  + "\n"
                    else:
                        strFile += "#source "+ source  + "\n"
                        self.addLog("#The specified source file does not exist: " + source)
        #Variables
        if shellType == SHELL_TYPE_BASH:
            strFile += "export AZORANGEHOME=%(installDir)s\n" % localVars
        else:
            strFile += "setenv AZORANGEHOME %(installDir)s\n" % localVars

        #  LD_LIBRARY_PATH  space separated paths in tcsh!!
        #LD_LIBPaths = [localVars["installDir"]+"/orange", localVars["installDir"]+"/orangeDependencies/bin"]
        LD_LIBPaths = [os.path.join("$AZORANGEHOME", "orange")]
        if "LD_LIBRARY_PATH" in self.EnvVars:
            for value in self.EnvVars["LD_LIBRARY_PATH"]:
                if value not in LD_LIBPaths: LD_LIBPaths.insert(0,value)
        self.EnvVars["LD_LIBRARY_PATH"] = LD_LIBPaths
        libPathsStr = ""
        for idx,value in enumerate(self.EnvVars["LD_LIBRARY_PATH"]):
            if idx: libPathsStr += ":"
            if str(value).startswith(self.AZOrangeInstallDir):
                value = "$AZORANGEHOME" + str(value)[len(self.AZOrangeInstallDir):] 
            libPathsStr += value
        if shellType == SHELL_TYPE_BASH:
            strFile += "if [[ ! -z $LD_LIBRARY_PATH ]] ; then\n"
            strFile += "    export LD_LIBRARY_PATH=\"%s:${LD_LIBRARY_PATH}\"\n" % libPathsStr
            strFile += "else\n"
            strFile += "    export LD_LIBRARY_PATH=\"%s\"\n" % libPathsStr
            strFile += "fi\n"
        else:
            strFile += "if ( $?LD_LIBRARY_PATH ) then\n"
            strFile += "    setenv LD_LIBRARY_PATH \"%s:${LD_LIBRARY_PATH}\"\n" % libPathsStr
            strFile += "else\n"
            strFile += "    setenv LD_LIBRARY_PATH \"%s\"\n" % libPathsStr
            strFile += "endif\n"

        #  PATH
        PATHPaths = [os.path.join("$AZORANGEHOME", "orangeDependencies", "bin"),"${PATH}"]
        if "PATH" in self.EnvVars:
            for value in self.EnvVars["PATH"]:
                if value not in PATHPaths: PATHPaths.insert(0,value)
        self.EnvVars["PATH"] = PATHPaths
        #  PYTHONPATH
        pythonPaths = [".",os.path.join("$AZORANGEHOME", "orange"), os.path.join("$AZORANGEHOME", "azorange"), os.path.join("$AZORANGEHOME", "tests")]
        if "PYTHONPATH" in self.EnvVars:      
            for value in self.EnvVars["PYTHONPATH"]:
                if value not in pythonPaths: pythonPaths.insert(0,value)
        self.EnvVars["PYTHONPATH"] = pythonPaths
        pythonPathsStr = ""
        for idx,value in enumerate(self.EnvVars["PYTHONPATH"]):
            if idx: pythonPathsStr += ":"
            if str(value).startswith(self.AZOrangeInstallDir):
                value = "$AZORANGEHOME" + str(value)[len(self.AZOrangeInstallDir):] 
            pythonPathsStr += value
        if shellType == SHELL_TYPE_BASH:
            strFile += "if [[ ! -z $PYTHONPATH ]] ; then\n"
            strFile += "    export PYTHONPATH=%s:${PYTHONPATH}\n" % pythonPathsStr
            strFile += "else\n"
            strFile += "    export PYTHONPATH=%s\n" % pythonPathsStr
            strFile += "fi\n"
        else:
            strFile += "if ( $?PYTHONPATH ) then\n"
            strFile += "    setenv PYTHONPATH %s:${PYTHONPATH}\n" % pythonPathsStr
            strFile += "else\n"
            strFile += "    setenv PYTHONPATH %s\n" % pythonPathsStr
            strFile += "endif\n"

        #Update and assign ALL Variables but "AZORANGEHOME" and "PYTHONPATH"
        for envVar in [x for x in self.EnvVars if x.upper() not in ["AZORANGEHOME" , "PYTHONPATH","LD_LIBRARY_PATH"]]:
            strValues = ""
            for idx,value in enumerate(self.EnvVars[envVar]):
                if idx: strValues += ":"
                strValues += value
            if shellType == SHELL_TYPE_BASH:
                strFile += "export %s=%s\n" % (envVar.upper() ,strValues)
            else:
                strFile += "setenv %s %s\n" % (envVar.upper() ,strValues)
        #Aliases
        strFile += "\n# AZOrange canvas alias\n"
        if shellType == SHELL_TYPE_BASH:
            strFile += "alias azorange='python " + os.path.join("$AZORANGEHOME", "orange", "OrangeCanvas", "orngCanvas.pyw") + "'\n"
        else:
            strFile += "alias azorange python " + os.path.join("$AZORANGEHOME", "orange", "OrangeCanvas", "orngCanvas.pyw") + "\n"
        #Modules
        if eval(self.modules):
            strFile += "\n# AZOrange module dependencies\n"
            if shellType == SHELL_TYPE_BASH:
                strFile += ". /etc/profile.d/modules.sh\n\n"
            else:
                strFile += "source /etc/profile.d/modules.csh\n\n"
            for mname in eval(self.modules):
                strFile += "module load %s\n" % (mname)
        else:
            strFile += "\n# Using NO modules!\n"

        #RunTimeModules
        if eval(self.RunTimeModules):
            strFile += "\n# AZOrange module dependencies needed at runtime only\n"
            for mname in eval(self.RunTimeModules):
                strFile += "module load %s\n" % (mname)
        else:
            strFile += "\n# Using NO specific modules for runtime only!\n"

        #Scripts to run upon setting the envitonment or loading the respective module
        strFile += "\n# Startup scripts\n"
        strFile += os.path.join("$AZORANGEHOME", "azorange", "bin", "clean.sh") + "\n"
        #strFile += os.path.join(self.AZOrangeInstallDir, "azorange/bin/ssh_testcfg.sh")  # This will be uncommented when using local mpi for the optimizer

        #Write the template file to current dir
        try:
            if shellType == SHELL_TYPE_BASH:
                localTemplateFile = os.path.join(self.currentDir,"templateProfile.bash")
            else:
                localTemplateFile = os.path.join(self.currentDir,"templateProfile")
            pyFile = open(localTemplateFile,"w")
            pyFile.write(strFile)
            pyFile.close()
        except:
            self.addLog((1,"Failed to create profile template file"))     
            return
        # #Write the template file to the user defined place (defined in setup.ini)
        if shellType == SHELL_TYPE_BASH:
            if localTemplateFile != self.TemplateProfileFileBash:
                self.__logAndExecute("cp -p" + localTemplateFile  +" "+ self.TemplateProfileFileBash)
            self.addLog("#Profile template file created in "+self.TemplateProfileFileBash)
        else:
            if localTemplateFile != self.TemplateProfileFile:
                self.__logAndExecute("cp -p" + localTemplateFile  +" "+ self.TemplateProfileFile)
            self.addLog("#Profile template file created in "+self.TemplateProfileFile)

        #Write the template file to the install dir depending on the installType
        if self.installType == "system" or self.repoInter=="export":
            self.addLog("#Profile template file copied into installDir")
            self.__logAndExecute("cp " + localTemplateFile + " " + self.AZOrangeInstallDir)
        else:            
            self.addLog("#Profile template file copied into trunkDir")
            self.__logAndExecute("cp " + localTemplateFile + " " + self.trunkDir)


    def InstallCacheCleaner(self):
        if self.cleanInstallCache:
            #For system installation, removes trunk and build dir
            #For developer installation, removes build dir
            self.addLog("*Cleaning install Cache")
            self.addLog("#Removing BuildDir")
            self.__logAndExecute("chmod -R 777 " + self.buildDir)
            self.__logAndExecute("rm -Rf " + self.buildDir)
            #export did not create trunkDir, so, there will be no trunk to delete
            #if self.installType == "system" and self.repoInter != "export":
            #    self.addLog("#Removing TrunkDir")
            #    self.__logAndExecute("chmod -R 777 " + self.trunkDir)
            #    self.__logAndExecute("rm -Rf " + self.trunkDir)


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
        #status can be the output of comands.getstatusoutput: (0, "something")
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
               self.successInstall = False
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
        log.close()
        allLogs.close()
        if self.successInstall:
            subject = "AZOrange Installation report: Success!"
        else:
            subject = "AZOrange Installation report: FAILED"
        for email in self.supportEmails.split(","):
            status = commands.getstatusoutput('mail -s "%s" %s < %s' % (subject, email.strip(), os.path.join(self.currentDir,"allLogs.log")) ) 
            if status[0] != 0:
                self.addLog(status)
                self.addLog("#WARNING: Could not send a report to %s.\nPlease contact support." % email.strip())
            else:
                self.addLog("#Report sent to "+email.strip())


    def __parse(self, arguments):
        opt = OptionParser(usage="%prog [options]")
        opt.add_option("-f", "--file", default="./setup.ini", dest='file', help='Path to INI file.')
        opt.add_option("-p", "--prepare", default=False, action="store_true", help="Prepare the install dir extracting/getting all 3rd party code into place. It does not start the installation procedure.")
        opt.add_option("-v", "--verbose", default=False, action="store_true", help="Enable verbose logging.")
        return opt.parse_args(arguments)


if __name__ == "__main__":
    import sys
    installer = Installer(sys.argv[1:])
    installer.main()
