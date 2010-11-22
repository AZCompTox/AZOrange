import ConfigParser
import os, sys
import time
import types,string
import commands

class installer:
    def __init__(self,configFile):
	if not os.path.isfile(configFile):
	    print "Missing the configuration file. "+configFile+" not found.\n You may use the setupTemplate.ini file to create your setup.ini file."
	    self.successInstall = False
	    return
	self.config = ConfigParser.ConfigParser()
        self.config.read(configFile)
        self.currentDir = os.getcwd()
	self.successInstall = True
        self.localVars = {}
        self.AZOrangeInstallDir = None
	self.readConfig()
	if not self.successInstall:
	    return
	self.addLog("Log for AZOrange installation",True)

    # Needed for modules to be loaded.
    def module(self,command, *arguments):
        commands = os.popen('%s/bin/modulecmd python %s %s' % (os.environ['MODULESHOME'], command, string.join(arguments))).read()
        exec commands
       
    def __update_EnvVars(self, envFile):
        """ This function will add to self.EnvVars  the variabled in the 'envFile'.
            For the vars that already exists, they will prepend the values in the file  to the values already 
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
                    #Assure that if something was added to pythonpath, that is also in the sys.path
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
                self.addLog(commands.getstatusoutput("chmod -R 777 " + self.trunkDir))
                self.addLog(commands.getstatusoutput("rm -Rf " + self.trunkDir))
            if not os.path.isdir(self.trunkDir): 
                self.addLog("#Creating TrunkDir")
                self.addLog(commands.getstatusoutput("mkdir -p " + self.trunkDir))
	#elif self.repoInter!="export":  #is "" or "update"
	#    if not os.path.isdir(os.path.join(self.trunkDir,".svn")):
	#        self.addLog("#WARNING: Seems that the trunkDir is not a subversion trunk!")
 
	#Prepare BuildDir
        if os.path.isdir(self.buildDir) and (self.preCleanBuildDir):
            self.addLog("#Removing existing buildDir")
            self.addLog(commands.getstatusoutput("chmod -R 777 " + self.buildDir))
            self.addLog(commands.getstatusoutput("rm -Rf " + self.buildDir))
	if not os.path.isdir(self.buildDir):
            self.addLog("#Creating buildDir")
            self.addLog(commands.getstatusoutput("mkdir -p " + self.buildDir))

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
                self.addLog(commands.getstatusoutput("chmod -R 777 " + self.installDir))
                self.addLog(commands.getstatusoutput("rm -Rf " + self.installDir))
            if not os.path.isdir(self.installDir):
                self.addLog("#Creating installDir")
                self.addLog(commands.getstatusoutput("mkdir -p " + self.installDir))
        if self.openInstallation:
            self.addLog(commands.getstatusoutput("mkdir -p " + self.DepSrcDir))

    def checkOutOpenAZO(self):
        if self.repoInter=="yes":
            #Update the AZO source from GITHub
            os.chdir(self.trunkDir)
            if self.openInstallation:
                self.addLog("*Using current files.")
                self.addLog("#trunk: "+self.trunkDir)
                #self.addLog(commands.getstatusoutput("git pull"))
            else:
                if self.installType == "developer":
                        self.addLog("*Cloning from GIT")
                        self.addLog("#trunk: "+self.trunkDir)
                        self.addLog(commands.getstatusoutput("rm -rf AZOrange"))
                        self.addLog(commands.getstatusoutput("git clone "+self.repo))
                        self.addLog(commands.getstatusoutput("cp -Rf AZOrange/* ."))
                        self.addLog(commands.getstatusoutput("cp -Rf AZOrange/.git " + self.trunkDir))
                        self.addLog(commands.getstatusoutput("rm -rf AZOrange"))
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
                        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"AZOrange")))
                        self.addLog(commands.getstatusoutput("tar xfz " + os.path.split(URL)[-1]))
                        self.addLog(commands.getstatusoutput("cp -Rf AZOrange/* " + self.trunkDir))
                        self.addLog(commands.getstatusoutput("cp -Rf AZOrange/.git " + self.trunkDir))
                        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"AZOrange")))
                #Revert any files replaced from the Git repo. We want to assure that the AZO inside files prevail
                os.chdir(install.trunkDir)
                self.addLog("#Reverting full RepoDir to SVN version...")
                self.addLog(commands.getstatusoutput("svn revert -R ./*")) 
                os.chdir(install.currentDir)

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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,name)))
        self.addLog(commands.getstatusoutput("mkdir " +  os.path.join(self.DepSrcDir,name)))
        os.chdir(os.path.join(self.DepSrcDir,"cdk"))
        jarFile = os.path.split(URL)[-1].strip()
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.addLog(commands.getstatusoutput("rm -rf " + jarFile))
                self.addLog(commands.getstatusoutput("wget " + URL ))
        else:
                self.addLog("*Using "+name+" in SVN Repo (Not implemented yet)")


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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,name)))
        os.chdir(self.DepSrcDir)
        tarFile = os.path.split(URL)[-1].strip()
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.addLog(commands.getstatusoutput("rm -rf " + tarFile))
                self.addLog(commands.getstatusoutput("wget " + URL ))
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
        self.addLog(commands.getstatusoutput(UnpackCmd + tarFile))
        self.addLog(commands.getstatusoutput("mv " + tarFile[0:tarFile.rfind(".tar")] + " " + name ))


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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,name)))
        os.chdir(self.DepSrcDir)
        tarFile = "cinfony_Download.tar.gz"
        if self.openInstallation:
                self.addLog("*Downloading "+name+" to trunk ("+URL+":"+REV+")")
                self.addLog(commands.getstatusoutput("rm -rf " + tarFile))
                # Download the File 
                self.addLog(commands.getstatusoutput("wget " + URL + " -O " + tarFile))
        else:
                self.addLog("*Using "+name+" in SVN Repo")
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
        self.addLog(commands.getstatusoutput(UnpackCmd))
        folderName = os.listdir("./"+name)[0]
        self.addLog(commands.getstatusoutput("mv "+os.path.join(name,folderName,"*") + " " + name))
        self.addLog(commands.getstatusoutput("rmdir "+os.path.join(name,folderName)))
        #self.addLog(commands.getstatusoutput("mv " + tarFile[0:tarFile.rfind(".tar")] + " " + name ))


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

        #self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.trunkDir,"orange/*")))
        #This command may have some failures, but it's no problem. We just want to delete if there is something to delete!
        self.addLog(commands.getstatusoutput("mkdir -p " + os.path.join(self.trunkDir,"orange")))
        if self.installType == "developer":
            commands.getstatusoutput('find ' + os.path.join(self.trunkDir,"orange") + '| grep -v "\.svn" | xargs rm -f')
        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"orange")))
        if self.openInstallation:
                self.addLog("*Checking out from orange SVN to trunk ("+URL+":"+REV+")")
                #os.chdir(os.path.join(self.trunkDir,"orange"))
                self.addLog(commands.getstatusoutput("mkdir -p " + os.path.join(self.DepSrcDir,"orange")))
                os.chdir(os.path.join(self.DepSrcDir,"orange"))
                self.addLog(commands.getstatusoutput("svn export --force" + " -r " + REV + " " + os.path.join(URL,"add-ons") + " ./add-ons"))
                self.addLog(commands.getstatusoutput("svn export --force" + " -r " + REV + " " + os.path.join(URL,"install-scripts/linux") + " ./install-scripts/linux"))
                self.addLog(commands.getstatusoutput("svn export --force" + " -r " + REV + " " + os.path.join(URL,"source") + " ./source"))
                self.addLog(commands.getstatusoutput("svn export --force" + " -r " + REV + " " + os.path.join(URL,"testing") + " ./testing"))
                self.addLog(commands.getstatusoutput("svn export --force" + " -r " + REV + " " + os.path.join(URL,"orange") + " ./ --force"))
        else:
                self.addLog("*Extracting orange from " + URL)
                os.chdir(self.DepSrcDir)
                self.addLog(commands.getstatusoutput("tar xfz " + os.path.split(URL)[-1]))
        os.chdir(self.DepSrcDir)
        self.addLog(commands.getstatusoutput('find '+os.path.join(self.DepSrcDir,"orange")+' -name ".svn" | xargs rm -Rf'))
        if self.installType == "developer":
            self.addLog(commands.getstatusoutput("cp -rf orange/* " + os.path.join(self.trunkDir,"orange/") ))
        else:
            self.addLog(commands.getstatusoutput("yes n | cp -R -i orange/* " + os.path.join(self.trunkDir,"orange/") ))

        # Apply Patch
        self.addLog("#Applying Patch...")
        os.chdir(os.path.join(self.currentDir,"Patches"))
        status,out = commands.getstatusoutput("./applyOrangePatch.sh %s" % (os.path.join(self.trunkDir,"orange")))
        if status != 0:
            self.addLog("#WARNING: Patch was not properly applied. It might be that it has been already applied.")
            self.addLog("#         Please check the next details.")
            self.addLog([0,out])
        else:
            self.addLog([status,out])
        #Revert to the Orange GIT files for the developer version since they were removed previously and wre copied with -rf flaga
        if self.installType == "developer":
            os.chdir(os.path.join(self.trunkDir,"orange"))
            self.addLog("#Reverting orange to Git version...")
            self.addLog(commands.getstatusoutput("git checkout -f ./"))

        if not self.openInstallation:
            os.chdir(os.path.join(self.trunkDir,"orange"))
            self.addLog("#Reverting orange to SVN version...")
            self.addLog(commands.getstatusoutput("svn revert -R ./*"))

	os.chdir(install.currentDir)


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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"appspack")))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Downloading AppsPack to trunk ("+URL+":"+REV+")")
                self.addLog(commands.getstatusoutput("rm -rf " + os.path.split(URL)[-1]))
                self.addLog(commands.getstatusoutput("wget " + URL ))
        else:   
                self.addLog("*Using APPSPACK in SVN Repo")
        self.addLog(commands.getstatusoutput("tar xfz " + os.path.split(URL)[-1]))
        self.addLog(commands.getstatusoutput("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " appspack" ))
        os.chdir(install.currentDir)


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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"opencv")))

        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Downloading opencv to trunk ("+URL+":"+REV+")")
                self.addLog(commands.getstatusoutput("rm -rf " + os.path.split(URL)[-1]))
                self.addLog(commands.getstatusoutput("wget " + URL ))
        else:
                self.addLog("*Using opencv in SVN Repo")
        self.addLog(commands.getstatusoutput("tar xfj " + os.path.split(URL)[-1]))
        self.addLog(commands.getstatusoutput("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " opencv" ))

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
        os.chdir(install.currentDir)

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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"oasa")))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Oasa is not needed in this installation")
                return
                # To download and install oasa in the open installation, remove both 2 previous lines and uncomment next 3 lines
                #self.addLog("*Downloading Oasa to trunk ("+URL+":"+REV+")")
                #self.addLog(commands.getstatusoutput("rm -rf " + os.path.split(URL)[-1]))
                #self.addLog(commands.getstatusoutput("wget " + URL ))
        else:
                self.addLog("*Using oasa in SVN Repo")
        self.addLog("#Extracting oasa from " + URL )
        self.addLog(commands.getstatusoutput("tar xfz " + os.path.split(URL)[-1]))
        self.addLog(commands.getstatusoutput("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " oasa" ))
        os.chdir(install.currentDir)


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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"boost")))
        os.chdir(self.DepSrcDir)
        if self.openInstallation:
                self.addLog("*Downloading Boost to trunk ("+URL+":"+REV+")")
                self.addLog(commands.getstatusoutput("rm -rf " + os.path.split(URL)[-1]))
                self.addLog(commands.getstatusoutput("wget " + URL ))
        else:
                self.addLog("*Using boost in SVN Repo")
        self.addLog(commands.getstatusoutput("tar xfz " + os.path.split(URL)[-1]))
        self.addLog(commands.getstatusoutput("mv " + os.path.split(URL)[-1][0:os.path.split(URL)[-1].rfind(".tar")] + " boost" ))

        os.chdir(install.currentDir)
 

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

        self.addLog(commands.getstatusoutput("rm -rf " + os.path.join(self.DepSrcDir,"plearn")))
        if self.openInstallation:
                self.addLog("*Checking out from plearn SVN to trunk ("+URL+":"+REV+")")
                self.addLog(commands.getstatusoutput("mkdir -p " + os.path.join(self.DepSrcDir,"plearn")))
                os.chdir(os.path.join(self.DepSrcDir,"plearn"))
                self.addLog(commands.getstatusoutput("svn export" + " -r " + REV + " " + URL + " ./ --force"))
        else:
                self.addLog("*Extracting PLearn from " + URL)
                os.chdir(self.DepSrcDir)
                self.addLog(commands.getstatusoutput("tar xfz " + os.path.split(URL)[-1]))
        self.addLog(commands.getstatusoutput('find '+os.path.join(self.DepSrcDir,"plearn")+' -name ".svn" | xargs rm -Rf'))

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
        os.chdir(install.currentDir)



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
        self.addLog(commands.getstatusoutput("cp -Rf " + os.path.join(self.trunkDir,"*") + " ./ "))
        # Remove all .svn and .git folders
        self.addLog(commands.getstatusoutput('find ./ -name ".svn" | xargs rm -Rf'))
        self.addLog(commands.getstatusoutput('find ./ -name ".git" | xargs rm -Rf'))
        
        #The orange core setup.py that is placed in orange/install-scripts/linux must be copied to orange
        self.addLog(commands.getstatusoutput("/bin/cp -f " + os.path.join(self.buildDir,"orange/install-scripts/linux/setup.py") + " " + os.path.join(self.buildDir,"orange")))

        #In buildDir> python setup.py build ...
        self.addLog("#Compiling by running setup.py")
	self.addLog("#Installation will be to made to " + self.AZOrangeInstallDir)
        envTmp = os.path.join(self.currentDir,"env.tmp")
        if os.path.isfile(envTmp):
            os.system("/bin/rm -f " + envTmp)
        self.addLog(commands.getstatusoutput('python setup.py build %s %s %s "%s" %s "%s" "%s" %s' % (self.platformType,self.buildDir,self.AZOrangeInstallDir,self.InstallTimeModules , envTmp , self.dependencies,self.GCCmoduleForAppspackMPI,self.openInstallation)))
        #Update the local self.EnvVars with the ones passed by setup.py which will be needed at runtime
        self.__update_EnvVars(envTmp) 

    def runAfterInstallScripts(self):	
        # This runs  once the .../azorange/bin/archive.sh in order to create the AZO_NFZ_scratchDir if it does not exists
        self.addLog("#Running after-install scripts")
        self.addLog(commands.getstatusoutput(os.path.join(self.AZOrangeInstallDir,'azorange/bin/archive.sh')))

    def createProfileExample(self):
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
        localVars = {"installDir":self.AZOrangeInstallDir}
        strFile =  "# Template profile for azorange installation at %(installDir)s\n" % localVars
        if self.sources:
            strFile += "\n# Additional sources needed\n"
            for source in self.sources:
                if os.path.isfile(source):
                    strFile += "source "+ source  + "\n"
                else:
                    strFile += "#source "+ source  + "\n"
                    self.addLog("#The specified source file does not exist: " + source)
        #Variables
        strFile += "setenv AZORANGEHOME %(installDir)s\n" % localVars

        #  LD_LIBRARY_PATH  space separated paths in tcsh!!
        #LD_LIBPaths = [localVars["installDir"]+"/orange", localVars["installDir"]+"/orangeDependencies/bin"]
        LD_LIBPaths = [localVars["installDir"]+"/orange"]
        if "LD_LIBRARY_PATH" in self.EnvVars:
            for value in self.EnvVars["LD_LIBRARY_PATH"]:
                if value not in LD_LIBPaths: LD_LIBPaths.insert(0,value)
        self.EnvVars["LD_LIBRARY_PATH"] = LD_LIBPaths
        libPathsStr = ""
        for idx,value in enumerate(self.EnvVars["LD_LIBRARY_PATH"]):
            if idx: libPathsStr += ":"
            libPathsStr += value
        strFile += "if ( $?LD_LIBRARY_PATH ) then\n"
        strFile += "    setenv LD_LIBRARY_PATH \"%s:${LD_LIBRARY_PATH}\"\n" % libPathsStr
        strFile += "else\n"
        strFile += "    setenv LD_LIBRARY_PATH \"%s\"\n" % libPathsStr
        strFile += "endif\n"

        #  PATH
        PATHPaths = [localVars["installDir"]+"/orangeDependencies/bin","${PATH}"]
        if "PATH" in self.EnvVars:
            for value in self.EnvVars["PATH"]:
                if value not in PATHPaths: PATHPaths.insert(0,value)
        self.EnvVars["PATH"] = PATHPaths
        #  PYTHONPATH
        pythonPaths = [".",localVars["installDir"]+"/orange", localVars["installDir"]+"/azorange", localVars["installDir"]+"/tests"]
        if "PYTHONPATH" in self.EnvVars:      
            for value in self.EnvVars["PYTHONPATH"]:
                if value not in pythonPaths: pythonPaths.insert(0,value)
        self.EnvVars["PYTHONPATH"] = pythonPaths
        pythonPathsStr = ""
        for idx,value in enumerate(self.EnvVars["PYTHONPATH"]):
            if idx: pythonPathsStr += ":"
            pythonPathsStr += value
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
            strFile += "setenv %s %s\n" % (envVar.upper() ,strValues)
        #Aliases
        strFile += "\n# AZOrange canvas alias\n"
        strFile += "alias azorange python %(installDir)s/orange/OrangeCanvas/orngCanvas.pyw\n" % localVars
        #Modules
        if eval(self.modules):
            strFile += "\n# AZOrange module dependencies\n"
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
        strFile += os.path.join(self.AZOrangeInstallDir, "azorange/bin/archive.sh\n")
        #strFile += os.path.join(self.AZOrangeInstallDir, "azorange/bin/ssh_testcfg.sh")  # This will be uncommented when using local mpi for the optimizer
        

        #Write the template file to current dir
        try:
            localTemplateFile = os.path.join(self.currentDir,"templateProfile")
            pyFile = open(localTemplateFile,"w")
            pyFile.write(strFile)
            pyFile.close()
        except:
	    self.addLog((1,"Failed to create profile template file"))     
            return
        # #Write the template file to the user defined place (defined in setup.ini)
        if localTemplateFile != self.TemplateProfileFile:
            self.addLog(commands.getstatusoutput("cp " + localTemplateFile  +" "+ self.TemplateProfileFile))
	self.addLog("#Profile template file created in "+self.TemplateProfileFile)
        #Write the template file to the install dir depending on the installType
	if self.installType == "system" or self.repoInter=="export":
	    self.addLog("#Profile template file copied into installDir")
	    self.addLog(commands.getstatusoutput("cp " + localTemplateFile + " " + self.AZOrangeInstallDir))
	else:	    
	    self.addLog("#Profile template file copied into trunkDir")
	    self.addLog(commands.getstatusoutput("cp " + localTemplateFile + " " + self.trunkDir))

	
    def InstallCacheCleaner(self):
	if self.cleanInstallCache:
	    #For system installation, removes trunk and build dir
	    #For developer installation, removes build dir
	    self.addLog("*Cleaning install Cache")
	    self.addLog("#Removing BuildDir")
	    self.addLog(commands.getstatusoutput("chmod -R 777 " + self.buildDir))
            self.addLog(commands.getstatusoutput("rm -Rf " + self.buildDir))
            #export did not create trunkDir, so, there will be no trunk to delete
            #if self.installType == "system" and self.repoInter != "export":
            #    self.addLog("#Removing TrunkDir")
            #    self.addLog(commands.getstatusoutput("chmod -R 777 " + self.trunkDir))
            #    self.addLog(commands.getstatusoutput("rm -Rf " + self.trunkDir))

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

 
if __name__ == "__main__":
    os.system("clear")
    PREPARE_ONLY = False
    if len(sys.argv) >= 2:
        if "prepare" in sys.argv[1].lower():
             PREPARE_ONLY = True
        else:
            print "Invalid input parameters."
            print "Usage:"
            print "    python install.py [option]"
            print "  option:"
            print "      prepare         - Prepare the install dir extracting/getting all 3rd party code into place"
            print "                        It does not start the installation procedure"
            sys.exit(1)


    install = installer("./setup.ini")
    
    if not install.successInstall:
	print "Nothing Done!"
	sys.exit(1)
    #Sends the configuration used to the details only!
    install.printConfiguration(True)
    startInstallTime = time.time()

    # Check for some important requirements
    install.addLog("*Requirements")
    if not os.path.isfile("/bin/tcsh"):
        install.addLog("#/bin/tcsh is missing. tcsh is required for azorange to work properly.")
        install.successInstall = False
    #=====install procedures=====
    if install.successInstall:
        install.prepareInstallDirs()

    #Checkout any 3rd party software to the proper locations
    if install.successInstall:
        install.checkOutOpenAZO()
    if install.successInstall:
        install.checkOutCinfony()
    if install.successInstall:
        install.checkOutRdkit()
    if install.successInstall:
        install.checkOutCDK()
    if install.successInstall:
        install.checkOutOrange() 
    if install.successInstall:
        install.checkOutAppspack()
    if install.successInstall:
        install.checkOutOpencv()
    if install.successInstall:
        install.checkOutBoost()
    if install.successInstall:
        install.checkOutPLearn()
    if install.successInstall:
        install.checkOutOasa()

       
    if install.successInstall: 
        if PREPARE_ONLY:
            install.addLog("*Everything is now unpacked and in place ready for install!")
        else:
            install.install()

    if not install.successInstall:
        install.addLog("*ERROR: Installation aborted!")

    if install.successInstall and not PREPARE_ONLY:
        install.createProfileExample()
        install.InstallCacheCleaner()
        install.runAfterInstallScripts()
	
    #==========================
    install.addLog("*Finished")
    install.addLog("#The process spent " + str(int((time.time() - startInstallTime)/60)) + " minutes")

    #get back to the initial directory
    os.chdir(install.currentDir)

    if not install.successInstall:
        install.addLog("#ERRORS during installation!!")
        install.emailLog()
        sys.exit(1)
    elif not PREPARE_ONLY :
        startAPPTemplate = os.path.join(install.trunkDir,'install/startAZOrange')
        startAPPTarget = os.path.join(install.AZOrangeInstallDir,'startAZOrange')
        AppLaucherTemplate = os.path.join(install.trunkDir,'install/AZOrange.desktop')
        AppLaucherTarget = None
        if "HOME" in os.environ:
            if os.path.isdir(os.path.join(os.getenv("HOME"),'Desktop')):
                AppLaucherTarget = os.path.join(os.getenv("HOME"),'Desktop')
        #Create the Application GUI Starter script
        cmd = 'sed "s|AZO_INSTALL_DIR|'+install.AZOrangeInstallDir+'|g" ' + startAPPTemplate + ' > '+startAPPTarget
        install.addLog(commands.getstatusoutput(cmd))
        install.addLog(commands.getstatusoutput("chmod a+x " +startAPPTarget))
        #Create a Launcher in the Desktop
        thisOS = commands.getstatusoutput("uname -o")
        if "gnu/linux" in thisOS[1].lower() and AppLaucherTarget:
            cmd='sed "s|AZO_INSTALL_DIR|'+install.AZOrangeInstallDir+'|g" '+ AppLaucherTemplate + ' > ' + os.path.join(AppLaucherTarget,'AZOrange.desktop')
            install.addLog(commands.getstatusoutput(cmd))
            install.addLog(commands.getstatusoutput("chmod a+x " + os.path.join(AppLaucherTarget,'AZOrange.desktop')))
        elif AppLaucherTarget:
            install.addLog(commands.getstatusoutput("ln -s " + startAPPTarget + " " + os.path.join(AppLaucherTarget,'AZOrange')))
            install.addLog(commands.getstatusoutput("chmod a+x " + os.path.join(AppLaucherTarget,'AZOrange')))
        install.addLog("#Installation done successfully!")
    else:
        install.addLog("#Preparation done successfully!")

    #Send the log if required in setup file
    install.emailLog()


    #Start the tests if required
    if install.runTests:
       os.system(os.path.join(install.currentDir,"runTests"))
 
