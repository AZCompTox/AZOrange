"""
Module to build install the azorange package. 
Specific versions of dependencies are needed, eg gcc.
"""
import sys
import os
import string
import commands
import types
import traceback

# Needed for modules to be loaded.
def module(command, *arguments):
    #status = commands.getstatusoutput('%s/bin/modulecmd python %s %s' % (os.environ['MODULESHOME'], command, string.join(arguments)))
    #print status
    commands = os.popen('%s/bin/modulecmd python %s %s' % (os.environ['MODULESHOME'], command, string.join(arguments))).read()
    exec commands
    

def checkStatus(status, output, message,check=True):
    if status != 0:
        print output
        print "\n\n++++++++++++++++++\n\n\n", message
        if check == True:
            sys.exit(1)
        else:
            print "Above Error ignored."


def syncCPATHandINCLUDE_DIR():
        """ This function should be called every time a change is made in CPATH or in INCLUDE_DIR"""
        CPATH = []
        INCLUDE_DIR = []
        # Get the current paths in both vars
        if "CPATH" in os.environ:
            CPATH = os.environ["CPATH"].split(":")
        if "INCLUDE_DIR" in os.environ:
            INCLUDE_DIR = os.environ["INCLUDE_DIR"].split(":")
        #Make CPATH handle ALL paths
        for path in INCLUDE_DIR:
            if path not in CPATH:
                CPATH.insert(0,path)
        #Assign the same paths to both variables (all are now present in CPATH)
        paths = ""
        #remove any empty paths
        INCLUDE_DIR=[]
        for idx,path in enumerate(CPATH):
            if path: INCLUDE_DIR.insert(0,path)
                
        for idx,path in enumerate(INCLUDE_DIR):
            if idx: paths += ":"
            paths += path
        os.environ["CPATH"] = paths 
        os.environ["INCLUDE_DIR"] = paths 



class installerClass:
    
    def __init__(self, rootDir, installDir, platform, modules):
        # rootDir is the current Build dir where this script is executed!
        # Dirs definition
        self.orangeInstallDir = installDir
        self.orangeDir = os.path.join(rootDir,"orange")
        self.azorangeDir = os.path.join(rootDir,"azorange")
        self.orngCRSDir = os.path.join(rootDir,"orangeDependencies/src/orngCRS")
        self.mpichDir = os.path.join(rootDir,"orangeDependencies/src/mpich-1.2.7p1")
        # The boost source path is the dir created when uncompressing the orangeDependencies/src/boost_1_34_1.tar.gz
        self.orngBoostDir = os.path.join(rootDir,"orangeDependencies/src/boost")
        self.APPSPackDir = os.path.join(rootDir,"orangeDependencies/src/appspack")
        self.azFannDir = os.path.join(rootDir,"orangeDependencies/src/azFann-2.0.0")
        self.opencvDir = os.path.join(rootDir,"orangeDependencies/src/opencv")
        self.oasaDir = os.path.join(rootDir,"orangeDependencies/src/oasa")
        self.plearnDir = os.path.join(rootDir,"orangeDependencies/src/plearn")
        self.R8Dir = os.path.join(rootDir,"orangeDependencies/src/R8/Src")
        self.trainingDir = os.path.join(rootDir,"azorange/trainingMethods")
        self.rootDir = rootDir
        # dir to install all other orange dependencies packages as opencv, fann, PLearn, etc...
        self.orangeDependenciesDir = os.path.join(installDir,"orangeDependencies")
        #make sure .../bin is also created
        os.system("mkdir -p " + os.path.join(self.orangeDependenciesDir,"bin"))

        #Variables
        self.buildPythonPath = self.orangeDir+":"+self.azorangeDir
        self.envFile = None     # new input argument to define additional envVars needed at runtime back to the mail install script
        self.mpiCVerOK=False    # used for MPI compatibility check
        self.modulesToLoad = modules
        self.platform = platform
        
        self.LoadModules()
        


    def __prependEnvVar(self,envName, value):
        """ This function will add to ./env.tmp a definition of the var 'envName' with the 'value'.
            If the envVar already exists in the file, it will prepend the value to the values already defined.
            If the envVar already have the value, it will NOT be doubled.
        """
        if not self.envFile:
            self.envFile = os.path.realpath("./env.tmp")
            print "WARNING!! It was not defines a file to report back new environment variables needed at run time.\n      " +\
                  self.envFile + "will be used."
        if not envName or not value or type(envName)!=types.StringType or type(value)!=types.StringType:
            print "Error appending envVar, envName and value cannot be empty, and have to be strings"
            sys.exit(1)
        envVars = {}
        envName = envName.upper()
        if os.path.isfile(self.envFile):
            file = open(self.envFile,"r")
            lines = file.readlines()
            file.close()
            for line in lines:
                envVars[line.split("=")[0].strip()] = line.split("=")[1].strip()
        if envName in envVars:
            values = [var.strip() for var in envVars[envName].split(":")]
            if value not in values:
                envVars[envName] = value + ":" + envVars[envName]
        else:
            envVars[envName] = value
        file = open(self.envFile,"w")
        for var in envVars:
            file.write(var+"="+envVars[var]+"\n")
        file.close()



    def LoadModules(self):
        # Setup up the correct env using modules for each target.
        ##scPA
        if self.modulesToLoad:
            sys.stdout.write("Loading module: ")
            #Check out what vital modules define in setup.ini  are missing
            if "LOADEDMODULES" in os.environ:
                missingModules = [m for m in self.modulesToLoad if m not in (os.environ["LOADEDMODULES"]).split(":")]
            else:
                missingModules = self.modulesToLoad
            #Load the missing modules
            for mname in missingModules:
                #print commands.getstatusoutput("gcc -dumpversion")
                sys.stdout.write("    %s " % mname)
                module("load", mname)
                #print commands.getstatusoutput("gcc -dumpversion")

            #Cheack if the needed modules were loaded correctly
            for mname in self.modulesToLoad:
                if mname not in (os.environ["LOADEDMODULES"]).split(":"):
                    print "\nWARNING!! Required Module '"+mname+"' was not able to be loaded...Installation may fail!"
        else:
            sys.stdout.write("Not using modules!")
        

    def installPythonLayer(self):
        """Once the orange/source directory is compiled, the python layer can be built using the 
           orange/setup.py script. """

        print "InstallPythonLayer"
        # Set the python path and install the python layer
        #setupFile = os.path.join(self.orangeDir,"setup.py")
        os.chdir(self.orangeDir)
        try:
            os.system("python setup.py install --orangepath="+self.orangeInstallDir)
        except:
            print "Error installing the orange in the installation directory."
            sys.exit(1)
        try:
            os.remove(os.path.join(self.orangeInstallDir,"orange","liborange.so"))
        except:
            print "Unable to remove orange/liborange.so"
        try:
            os.symlink(os.path.join(self.orangeInstallDir,"orange","orange.so"), os.path.join(self.orangeInstallDir,"orange","liborange.so"))
        except:
            print "Failed to link liborange.so"
            sys.exit(1)
        # Also copy the shared libraries that were built separately.
        #print "Copying c45.so and _orngCRS.so to installation directory."
        #try:
        #    os.system("cp c45.so %s/orange/c45.so" % self.orangeInstallDir)
        #    os.system("cp _orngCRS.so %s/orange/_orngCRS.so" % self.orangeInstallDir)
        #except:
        #    print "Error copying c45.so and/or _orngCRS.so from orange build directory to installation directory."

        print "Copying trainingMethods, AZutilities and AZOrangeConfig.py to installation directory."
        sumStatus = 0
        try:
	    ##scPA    Changed the source and destination of these files
            sumStatus += os.system("mkdir -p %s/azorange/trainingMethods" % self.orangeInstallDir)
            #Not included in the sumStatus because this command will generate a warning because being sipping directories.
            os.system("cp -f ../azorange/trainingMethods/* %s/azorange/trainingMethods/" % self.orangeInstallDir)
            #Not being used... it is still empty
            #sumStatus += os.system("cp -f ../azorange/trainingMethods/bin/* %s/azorange/trainingMethods/bin/" % self.orangeInstallDir)
            sumStatus += os.system("cp -Rf ../azorange/AZutilities %s/azorange/." % self.orangeInstallDir)
            sumStatus += os.system("cp -Rf ../azorange/statlib %s/azorange/." % self.orangeInstallDir)
            sumStatus += os.system("cp -f ../azorange/*.py %s/azorange/." % self.orangeInstallDir)

            if not self.OpenInstallation:
                # The azorange/bin may be needed for the appspackMPI for the open installation
                sumStatus += os.system("cp -Rf ../doc %s/." % self.orangeInstallDir)
                sumStatus += os.system("cp -Rf ../azorange/bin %s/azorange/." % self.orangeInstallDir)
                sumStatus += os.system("cp -Rf ../azorange/documentation %s/azorange/." % self.orangeInstallDir)
                sumStatus += os.system("cp -f ../azorange/*.txt %s/azorange/." % self.orangeInstallDir)
                sumStatus += os.system("cp -Rf  ../exampleScripts %s/." % self.orangeInstallDir)
	    #  Added the azorange to pythonpath in order to maintain the use of modules
            #  inside it accessible by the same way
            if "PYTHONPATH" in os.environ:
                os.environ["PYTHONPATH"] =  self.orangeDir + "/azorange:" +  os.environ["PYTHONPATH"]
            else:
                os.environ["PYTHONPATH"] =  self.orangeDir + "/azorange"
	    ##ecPA
            #print "copy env:"
            #os.system("ls -l ")
            #os.system("pwd")
        except:
            print "Error copying trainingMethods directory from orange build directory to installation directory."
        if sumStatus:
            print "Error copying azorange files to the install Dir."
            sys.exit(1)
        

    def test(self, installDir):
 
        print "Testing the installation!"
        self.testsDir = os.path.join(installDir,"azorange/azorange/tests")

        self.defPaths(installDir)
        os.putenv("PYTHONPATH", self.PYTHONPATH)
        os.putenv("LD_LIBRARY_PATH", self.LD_LIBRARY_PATH) 

        # Execute in the shell to get the PYTHONPATH defined in installPythonLayer
        print "Testing the ANN implementation!!!"
        os.system("python "+self.testsDir+"/AZorngANNTest.py")
        print "Testing the SVM implementation!!!"
        os.system("python "+self.testsDir+"/AZongSVMTest.py")
        print "Testing the RF implementation!!!"
        os.system("python "+self.testsDir+"/AZorngRFTest.py")
        print "Testing the PLS implementation!!!"
        os.system("python "+self.testsDir+"/AZorngPLSTest.py")


    def compileOrange(self):
        # Orange is being compiled three times because there is something wrong with how the dependencies and the entire orange make procedure is set up.
        print "CompileOrange"
        os.chdir(self.orangeDir)
        stat=1
        c=0
        stat, out = commands.getstatusoutput("python setup.py compile")
        #while stat and c<10:
        #    print "Running orange setup.py ... try " + str(c+1)
        #    stat, out = commands.getstatusoutput("python setup.py compile")
        #    c=c+1

        # checkStatus(stat, out,"Error compiling orange.") Can't check
        # error here since there is one due to the faulty make
        # procedure..
        #stat, out = commands.getstatusoutput("python setup.py compile")
        checkStatus(stat, out,"Error compiling orange.")


    def compileOrngCRS(self):
        if "orngcrs" not in self.dependencies:
            print "Not using the local orngCRS"
            return
        if self.dependencies["orngcrs"]:   # compile and install
                print "compile OrngCRS"
                os.chdir(self.orngCRSDir+"/src")
                stat,out = commands.getstatusoutput("make")
                checkStatus(stat, out,"Error compiling _orngCRS.")
                
                stat, out = commands.getstatusoutput("cp -f _orngCRS.so " + self.orangeDir + "/_orngCRS.so")
                checkStatus(stat, out,"Error copying _orngCRS.so to orange build directory.")
        else:
                print "Not reinstalled"

    def compileR8(self):
        if "r8" not in self.dependencies:
            print "Not using the local R8"
            return
        if self.dependencies["r8"]:   # compile and install
                print "compileR8"
                print "Need to tell the buildC45 script where orange resides."
                if "PYTHONPATH" in os.environ:
                    os.environ["PYTHONPATH"] =  self.orangeDir + ":" +  os.environ["PYTHONPATH"] 
                else:
                    os.environ["PYTHONPATH"] =  self.orangeDir
                if "LD_LIBRARY_PATH" in os.environ:
                    os.environ["LD_LIBRARY_PATH"] =   self.orangeDir + ":" + os.environ["LD_LIBRARY_PATH"]  
                else:
                    os.environ["LD_LIBRARY_PATH"] =   self.orangeDir

                os.chdir(self.R8Dir)
                stat, out = commands.getstatusoutput("cp -f %s/buildC45.py ." % self.orangeDir)
                checkStatus(stat, out,"Error copying build script for R8 (buildC45.py).")
                stat, out = commands.getstatusoutput("cp -f %s/ensemble.c ." % self.orangeDir)
                checkStatus(stat, out,"Error copying ensemble.c.")
                print "Build c45 start"
                stat, out = commands.getstatusoutput("python buildC45.py")
                checkStatus(stat, out,"Error compiling R8 (c45).")
                print "Build c45 end"
        else:
                print "Not reinstalled"

    def compileMPICH(self):
        if "mpich.tar.gz" not in self.dependencies:
            print "Not using the local mpich"
            return
        mpichInstallDir = os.path.join(self.orangeDependenciesDir,os.path.split(self.mpichDir)[1])
        if self.dependencies["mpich.tar.gz"]:   # compile and install
                # uncompress the .tar.gz in the <rootDir>/orangeDependencies/src dir, 
                #which will create a dir with name self.mpichDir
                os.chdir(os.path.join(self.rootDir,"orangeDependencies/src"))
                stat,out = commands.getstatusoutput("tar xfz mpich.tar.gz")
                checkStatus(stat, out,"Error extracting 'MPICH' source files")
                #Compile
                print "Compiling MPICH"
                os.chdir(self.mpichDir)
                print "Building in:   ",self.mpichDir
                print "Installing in: ",mpichInstallDir
                stat, out = commands.getstatusoutput("./configure --prefix=\"" + mpichInstallDir  + "\"")
                checkStatus(stat, out,"Error configuring mpich.")
                stat, out = commands.getstatusoutput("make")
                checkStatus(stat, out,"Error compiling mpich.")
                stat, out = commands.getstatusoutput("make install")
                checkStatus(stat, out,"Error installing mpich.")
        else:   
                print "Not reinstalled"
        # At runtime we will need PATH to include the location of new installed mpich
        self.__prependEnvVar("PATH" , os.path.join(mpichInstallDir,"bin"))
        #The version will also be tested when compiling  appspack
        os.environ["PATH"] =  os.path.join(mpichInstallDir,"bin") + ":" + os.environ["PATH"]
        sys.path.insert(0 , os.path.join(mpichInstallDir,"bin"))
 
        

    def compileAPPSPack(self):
        if "appspack" not in self.dependencies:
            print "Not using the local appspack"
            return
        if not self.dependencies["appspack"]:
            print "Not reinstalled"
            return 
        # Also copy appspack that was built separately.
        if self.platform != "workstation" and self.platform[0:3] != "GAS":
            print "No installation of APPSPACK on the cluster."
            return

        # Compile only APPSPACK with gcc 3 (system default) as we have not been able to build mpi version with gcc 4
        # Find the GCC correct module name on this system if any.
        gccModuleName = None
        if self.modulesToLoad:
            for moduleName in self.modulesToLoad:
                if "gcc" in moduleName.lower():
                    gccModuleName = moduleName
                    #Unload the GCC module
                    print "Unloading module ",gccModuleName
                    module("unload", gccModuleName) 
                    break

        if self.GCCmoduleForAppspackMPI:
            print "Loading module ",self.GCCmoduleForAppspackMPI
            module("load", self.GCCmoduleForAppspackMPI)


        print "Compiling APPSPack."
        print commands.getstatusoutput("env")[1]
        print commands.getstatusoutput("which mpicxx")[1]
        os.chdir(self.APPSPackDir)
        print "Building in:   ",self.APPSPackDir
        #minVer = "4.1.2"  # Version of mpich. This should be controlled by the gcc version selected for appspack. To be done
        minVer = "3.4.6"
        mpiC = "mpicxx"
        status = commands.getstatusoutput(mpiC + " -dumpversion")
        if status[0]!=0:
            self.mpiCVerOK=False
            version = None
        else:
            self.mpiCVerOK=True
            version = status[1] 
            minVerList = minVer.split(".")
            runningVerList = version.split(".")
            subVerN = min(len(minVerList),len(runningVerList))

            for idx in range(subVerN):
                if runningVerList[idx] < minVerList[idx]:
                    print "MPI set to false here"
                    self.mpiCVerOK = False
                    break
                elif runningVerList[idx] > minVerList[idx]:
                    break

        ################################################################
        # Report
        ################################################################
        print "#######  appspack MPI compiler #######"
        if self.mpiCVerOK :
            print "Running version '"+version+"' of "+mpiC+" (This version should be compatible!)"
        else:
            if not version:
                print "WARNING: "+mpiC+" was not found: "+str(status)
            else:
                print "WARNING: Running version '"+version+"' of "+mpiC+"\n         At least version ",minVer, " is required!"
            print "WARNING: MPI version of APPSPACK for the optimizer  will not be available!"

        if self.mpiCVerOK: 
            mpiFlags = "--with-mpi-compilers"
        else:
            mpiFlags = ""
        if self.platform[0:3] == "GAS":
            stat, out = commands.getstatusoutput("./configure CC=gcc CXX=g++ --prefix=%s --with-blas=%s --with-lapack=%s %s" % (self.orangeDependenciesDir,os.environ["ATLAS"],os.environ["ATLAS"],mpiFlags))
        else:
            stat, out = commands.getstatusoutput("./configure CC=gcc CXX=g++ --prefix=%s %s" % (self.orangeDependenciesDir,mpiFlags))
            
        checkStatus(stat, out,"Error configuring APPSPack.")

        # Make sure we are clean at the start.
        stat, out = commands.getstatusoutput("make clean")
        checkStatus(stat, out,"Error cleaning APPSPack")

        # Appspack cant be built in paralell
        stat, out = commands.getstatusoutput("make -j1 ")
        checkStatus(stat, out,"Error compiling APPSPack")
                 
        # Cp the binaries to the install bin dir
        os.system("cp "+self.APPSPackDir + "/src/appspack_serial " + os.path.join(self.orangeDependenciesDir,"bin"))
        # make it executable
        stat, out = commands.getstatusoutput("chmod a+x %s" % os.path.join(self.orangeDependenciesDir,"bin/appspack_serial"))  
        checkStatus(stat, out,"Error making APPSPack_serial  executable")

        if self.mpiCVerOK:
            os.system("cp "+self.APPSPackDir + "/src/appspack_mpi " + os.path.join(self.orangeDependenciesDir,"bin"))
            # make it executable
            stat, out = commands.getstatusoutput("chmod a+x %s" % os.path.join(self.orangeDependenciesDir,"bin/appspack_mpi"))
            checkStatus(stat, out,"Error making APPSPack_mpi  executable")

        # Compile the rest of the code with gcc 4
        # Unload the module used to compile appspack
        if self.GCCmoduleForAppspackMPI:
            print "Unloading module ",self.GCCmoduleForAppspackMPI
            module("unload", self.GCCmoduleForAppspackMPI)

        # Load the gcc module if it was unloaded before
        if gccModuleName:
            print "Reloading module ",gccModuleName
            module("load", gccModuleName)
        

    def getDirOfFile(self, fileName):
        """Return the directory in which fileName resides (as returned by the locate command)"""
        try:             
            filePath = commands.getoutput("locate "+fileName).split("\n")[-1]
        except:
            print "Unable to locate the file."
            return None
        dirPath = os.path.dirname(filePath)
        if os.path.isdir(dirPath):
            return dirPath 
        else: 
            print "WARNING: The search of '"+fileName+"' returned:",filePath
            print "         Try to update your DB first by running: sudo updatedb"
            return None


    def setEnv(self):
        print "The rootDir: ",self.rootDir
        # Add nspr directory to CPATH
        if self.platform[0:3] == "GAS":
	    try:
                pcfile=open("/usr/lib/pkgconfig/mozilla-nspr.pc")
                line=pcfile.readline()
                while(line):
                    if line.find("Cflags:") == 0:
                        nsprPath=line[10:].rstrip()
                    line=pcfile.readline()
                pcfile.close()
            except:
                nsprPath = self.getDirOfFile("prlong.h")
                if not nsprPath:
                    print "ERROR: NSPR is missing! Please install NSPR and try again."
                    sys.exit(1)

        else:
            nsprPath = self.getDirOfFile("prlong.h")
            if not nsprPath:
                print "ERROR: NSPR is missing! Please install NSPR and try again."
                sys.exit(1)

            idx = string.rfind(nsprPath, "/")
            nsprPath = nsprPath[:idx]

        print "nsprPath "+nsprPath
        if nsprPath.find("nspr4") != -1:
            # nspr4 not what our code expect.. Make fake link in .
            newnsprPath=os.path.realpath("nspr") #Find complete path
            print "Faking nspr path to %s to workaround the nspr4 problem" % newnsprPath
            try:
                if not os.path.exists(newnsprPath):
                    os.symlink(nsprPath, newnsprPath) # Create a link to the nspr4 include files
                nsprPath=newnsprPath+":"+os.path.realpath(".") # Set nsprPath to link... and include . 
            except:
                print "Error creating the symlink to nspr include dir. You need to be su to create the link."
                sys.exit(1)


        if nsprPath:
            print os.environ["CPATH"]
            print nsprPath
            if "CPATH" in os.environ:
                os.environ["CPATH"] =  nsprPath + ":"+ os.environ["CPATH"]
            else:
                os.environ["CPATH"] =  nsprPath
            syncCPATHandINCLUDE_DIR()
            print "Setting env CPATH to %s to make azorange find nspr include files" % os.environ["CPATH"]
        else:
            print "Terminating AZOrange installation because the nspr package was not found."
            sys.exit(1)

    def compileOasa(self):
        # Remove "or not self.OpenInstallation" from next test if want to install oasa in the open installation
        if ("oasa" not in self.dependencies) or self.OpenInstallation:
            print "Not using the local oasa"
            return
        oasainstallDir = os.path.join(self.orangeDependenciesDir,os.path.split(self.oasaDir)[1])
        if self.dependencies["oasa"]:   #compile and install 
                print "Compiling oasa"
                os.chdir(self.oasaDir)
                print "Building in:   ",self.oasaDir
                print "Installing in: ",oasainstallDir
                # for debug include:    --enable-debug=\"yes\"
                stat, out = commands.getstatusoutput("python setup.py build")
                checkStatus(stat, out,"Error compiling oasa.")
                stat, out = commands.getstatusoutput("python setup.py install --prefix=\"" + oasainstallDir  + "\"")
                checkStatus(stat, out,"Error installing oasa")
        else:
                print "Not reinstalled"

        status, out = commands.getstatusoutput("find \""+oasainstallDir+"\" | grep molecule.py | grep -v \".svn\"")
        checkStatus(status, out,"Error: molecule.py not found while finding site-packages of oasa!")
        sitePackagesPath = os.path.split(os.path.split(out.split("\n")[-1])[0])[0]
        if not sitePackagesPath:
            print "ERROR: molecule.py not found while finding site-packages of oasa"
            sys.exit(1)

        # At runtime we will need PYTHONPATH to include the location of new installed modules
        self.__prependEnvVar("PYTHONPATH" , sitePackagesPath)


    def compileOpenCV(self):
        if ("opencv" not in self.dependencies):
            print "Not using the local openCV"
            return
        openCVinstallDir = os.path.join(self.orangeDependenciesDir,os.path.split(self.opencvDir)[1])
        if self.dependencies["opencv"]:   #compile and install 
                print "Compiling openCV"
                os.chdir(self.opencvDir)
                print "Building in:   ",self.opencvDir
                print "Installing in: ",openCVinstallDir
                # for debug include:    --enable-debug=\"yes\"
                stat, out = commands.getstatusoutput("./configure --prefix=\"" + openCVinstallDir  + "\" --with-openmp --enable-apps=\"no\"")
                checkStatus(stat, out,"Error configuring opencv.")
                stat, out = commands.getstatusoutput("make")
                checkStatus(stat, out,"Error compiling opencv.")
                stat, out = commands.getstatusoutput("make install")
                checkStatus(stat, out,"Error installing opencv.")
        else:
                print "Not reinstalled"

        #This env var will be used for compiling AZOrange
        os.environ["PKG_CONFIG_PATH"] =  os.path.join(openCVinstallDir,"lib/pkgconfig")
        # At runtime we will need LD_LIBRARY_PATH to include the location of new installed libs
        self.__prependEnvVar("LD_LIBRARY_PATH" , os.path.join(openCVinstallDir,"lib"))
 
        #This env var will be used for compiling AZOrange
        status, out = commands.getstatusoutput("find \""+openCVinstallDir+"\" | grep cv.py | grep -v \".svn\"")
        checkStatus(status, out,"Error: cv.py not found while finding site-packages of opencv!")
        sitePackagesPath = os.path.split(os.path.split(out.split("\n")[-1])[0])[0]
        if not sitePackagesPath:
            print "ERROR: cv.py not found while finding site-packages of opencv"
            sys.exit(1)
            
        if "PYTHONPATH" not in os.environ:
            os.environ["PYTHONPATH"] =  sitePackagesPath
        else:
            os.environ["PYTHONPATH"] =  sitePackagesPath+":"+os.environ["PYTHONPATH"]
        # At runtime we will need PYTHONPATH to include the location of new installed modules
        self.__prependEnvVar("PYTHONPATH" , sitePackagesPath)


    def compileAZFann(self):
        if ("azfann-2.0.0" not in self.dependencies):
            print "Not using the local Fann"
            return
        FANNinstallDir = os.path.join(self.orangeDependenciesDir,os.path.split(self.azFannDir)[1])
        if self.dependencies["azfann-2.0.0"]:   #compile and install
                print "Compiling Fann"
                os.chdir(self.azFannDir)
                print "Building in:   ",self.azFannDir
                print "Installing in: ",FANNinstallDir
                stat, out = commands.getstatusoutput("./configure --prefix=\"" + FANNinstallDir  + "\"")
                checkStatus(stat, out,"Error configuring Fann.")
                stat, out = commands.getstatusoutput("make")
                checkStatus(stat, out,"Error compiling Fann.")
                stat, out = commands.getstatusoutput("make install")
                checkStatus(stat, out,"Error installing Fann.")

                print "compiling pyfann" # the python bindings
                os.chdir(os.path.join(self.azFannDir,"python"))
                stat, out = commands.getstatusoutput("make")
                checkStatus(stat, out,"Error compiling pyFann.")
                stat, out = commands.getstatusoutput("python setup.py install --install-lib " + os.path.join(FANNinstallDir,"lib"))
                checkStatus(stat, out,"Error installing pyfann.")
        else:
                print "Not reinstalled"
        # At runtime we will neew LD_LIBRARY_PATH to include the location of new installed libs
        self.__prependEnvVar("LD_LIBRARY_PATH" , os.path.join(FANNinstallDir,"lib"))
        self.__prependEnvVar("PYTHONPATH" , os.path.join(FANNinstallDir,"lib/pyfann"))


    def compileAZOrange(self):
        # These settings are very specific and might cause problems when eg the compiler is changed.
        print "compileAZOrange..."
        if "plearn" not in self.dependencies:
                print "Not using the local plearn"
                return
        else:
                os.environ["PLEARNDIR"] =  self.plearnDir
                if "PYTHONPATH" in os.environ:
                    os.environ["PYTHONPATH"] =  os.environ["PLEARNDIR"] + "/python_modules:" + os.environ["PYTHONPATH"]
                else:
                    os.environ["PYTHONPATH"] =  os.environ["PLEARNDIR"] + "/python_modules"
                os.environ["PATH"] =  os.environ["PLEARNDIR"] + "/scripts:"+ os.environ["PLEARNDIR"] + "/commands:"+os.environ["PATH"]
                if "CPATH" in os.environ:
                    os.environ["CPATH"] =  os.environ["PLEARNDIR"] + ":" + os.environ["CPATH"] 
                else:
                    os.environ["CPATH"] =  os.environ["PLEARNDIR"]
                syncCPATHandINCLUDE_DIR()
                # Make sure that whoever is running this does not have $HOME/.pymake or $HOME/.plearn.
                if ( os.path.isdir(os.environ['HOME']+"/.pymake") or os.path.isdir(os.environ['HOME']+"/.plearn") ):
                    print "You can't build azorange with .pymake or .plearn in your home directory."
                    sys.exit(1)

        if "boost" not in self.dependencies:
                print "Not using the local boost"
        else:
                #include for compiling the source files of boost
                os.environ["CPATH"] = self.orngBoostDir + ":" +  os.environ["CPATH"]
                syncCPATHandINCLUDE_DIR()

        print "Environment when compiling azorange:"
        print commands.getstatusoutput("env")[1]
        os.chdir(os.path.join(self.trainingDir, "src"))


        # Build and bind to python layer

        #opencv now have the python bindings 
        #print "Comping randomForest.cc"
        #stat,out = commands.getstatusoutput("g++ -c -fPIC `pkg-config --cflags opencv` -o randomForest.o randomForest.cc")
        #checkStatus(stat, out,"Error compiling randomForest")

        #print "swig -c++ -python randomForest.i"
        #stat,out = commands.getstatusoutput("swig -c++ -python randomForest.i")
        #checkStatus(stat, out,"Error binding randomForest")

        #print "Compile randomForest_wrap.cxx"
        #stat,out = commands.getstatusoutput("g++ -c -fPIC randomForest_wrap.cxx `pkg-config --cflags opencv`")
        #checkStatus(stat, out,"Error compiling randomForest_wrap.cxx")

        #print "Creating shared object _randomForest.so"
        #stat,out = commands.getstatusoutput("g++ -shared `pkg-config --libs opencv` randomForest.o randomForest_wrap.o -o _randomForest.so")
        #checkStatus(stat, out,"Error creating _randomForest.so")

        print " Compile PLS.cc and PlsAPI.cc"
        stat,out = commands.getstatusoutput("pymake -nopython PLS.cc")
        checkStatus(stat, out,"Error compiling PLS.cc")
        stat,out = commands.getstatusoutput("pymake -nopython PlsAPI.cc")
        checkStatus(stat, out,"Error compiling PlsAPI.cc")

        print "swig -c++ -python PlsAPI.i"
        stat,out = commands.getstatusoutput("swig -c++ -python PlsAPI.i")
        # In python layer we are NEVER using anything of PLearn. None of the PLearn methods are 
        # binding to python layer, and from python layer we will never use any object of PLearn.  Just the API methodas will
        # interact with PLearn, but they were compiled just before, and now we are only binding some of those methods.
        #if stat and "Nothing known about namespace 'PLearn'" not in out:
        checkStatus(stat, out,"Error binding PlsAPI")

        # SCPA ====================================== TMP DISABLING =================
        print "Compiling PlsAPI_wrap.cxx"
        stat,out = commands.getstatusoutput("pymake -so -nopython PlsAPI_wrap.cxx")
        checkStatus(stat, out,"Error compiling PlsAPI_wrap.cxx")
        # Move required binaries to the bin dir  (INSTALLING)
        sumStatus = 0
        try:
            sumStatus += os.system("/bin/mv pls.py ../.")
            sumStatus += os.system("/bin/cp -fL libPlsAPI_wrap.so ../_pls.so")
        except:
            print "Error copying azorange files after compilation."
            sys.exit(1)
        if sumStatus:
            print "Error copying azorange files after compilation. Probably the compilation was not successful."
            sys.exit(1)
        # ECPA ======================================================================= 
        

    def build(self):
        """Install azorange in rootDir. """

        #os.system("env")
        #sys.exit(1)
        if "LOADEDMODULES" in os.environ:
            print "Loaded modules: %s"%os.environ["LOADEDMODULES"]
        else:
            print "No Modules loaded!"
        print "Orange is being built in "+self.rootDir
        
        # Setup the correct environment.
        self.setEnv()


        # Compile the Orange C++ layer
        self.compileOrange()
    
        # R8 and copy the C4.5.so to the orange dir.
        self.compileR8()
        
        # APPSPack.  moved to last
        
        # Compile and copy _orngCRS.so to the orange dir.
        self.compileOrngCRS()

        # Compile and install the opencv.
        self.compileOpenCV()

        # Compile and install the oasa.
        self.compileOasa()

        # Compile and install the azFann.
        self.compileAZFann()

        # Compile and install the AZOrange C++ layer
        self.compileAZOrange()

        # APPSPack.
        if self.platform == "workstation" or self.platform[0:3] == "GAS":
            self.compileMPICH()
            self.compileAPPSPack()




        # Set the python path and install the python layer
        self.installPythonLayer()


        # Copy the test directory
        os.system("cp -r "+os.path.join(self.rootDir, "tests")+" "+self.orangeInstallDir)
        #Expand profiling Data
        if not self.OpenInstallation: 
            #Dir of dataSets used by profiling tools
            dataDir = os.path.realpath(os.path.join(self.orangeInstallDir,'tests','profiling','dataSets'))
            # Zipped file containing the full set of datasets for profiling
            ZippedData=os.path.realpath(os.path.join(dataDir,'profilingDataSuite.tar.gz'))
            #Dir where the original datasets are going to be unzipped
            fullTrunkDir=os.path.join(dataDir,'fullTrunk')
            os.system("rm -rf " + fullTrunkDir)
            os.system("mkdir -p " + fullTrunkDir)
            os.system("tar xfz " + ZippedData + " -C " + fullTrunkDir)



if __name__ == "__main__":
    """Installation procedure for AZOrange. 
       python setup.py build workstation/GAS/cluster build_dir install_dir
       """
    ##scPA
    # Modules to load specified in the setup.ini

    if len(sys.argv) < 10:
        print "Missing input arguments. Expected 8 but got "+str(len(sys.argv))
        sys.exit(1)
    modulesToLoad = eval(sys.argv[5])
    if sys.argv[6].lower() == "none":
        envFile = None
    else:
        envFile = sys.argv[6]

    # Only compile and install if the dependency is present in the dict passed and is set to True.
    # If it is set to False, do not compile or install, but keep the envVars pointing to the supposed installation place
    # If not present at all, the dependency is not to be used, it is supposed that it is already installed and loaded
    dependencies = eval(sys.argv[7]) 

    # GCC module to be used when compiling appspack MPI
    if sys.argv[8].lower() == "none":
        GCCmoduleForAppspackMPI = None
    else:   
        GCCmoduleForAppspackMPI = sys.argv[8]
    if sys.argv[9].lower() == "true":
        OpenInstallation = True
    else:
        OpenInstallation = False
    ##ecPA

    try:
        action = sys.argv[1]
    except: 
        print "No action was specified for the setup script."
        print "Building in current directory."
        action = "build"

    gotDir=False
    ##scPA
    syncCPATHandINCLUDE_DIR()
    ##ecPA
    if action == "build":
        try:
            platform = sys.argv[2]
            if platform == "workstation": 
                print "Building workstation version."
            elif platform[0:3] == "GAS": 
                print "Building %s version." % platform
            elif platform == "cluster": 
                print "Building cluster version."
            else:
                print "ERROR"
                print "The second argument must be 'workstation', 'cluster' or 'GAS-py2.[45]'" 
                print "AZOrange is not being built."
                sys.exit(1)
            try: 
                buildDir=sys.argv[3]
                installDir=sys.argv[4]
                if os.path.exists(installDir) and os.path.exists(buildDir): gotDir = True
            except: 
                print "No directory is specified for the installation of AZOrange."
                print "AZOrange could be run from current directory."
                gotDir = False
            if gotDir:
                installer = installerClass(buildDir, installDir, platform, modulesToLoad)
                installer.envFile = envFile
                installer.dependencies = dependencies
                installer.GCCmoduleForAppspackMPI = GCCmoduleForAppspackMPI
                installer.OpenInstallation = OpenInstallation
                installer.build() 
            else:
                print "AZOrange is not installed because the build or install directory does not exist."
                sys.exit(1)
        except: 
            traceback.print_exc() 
            print "AZOrange is not built."
            sys.exit(1)
    else:
        print "The first argument bust be build."
        print "Only build argument is supported and it also installs."
        sys.exit(1)


