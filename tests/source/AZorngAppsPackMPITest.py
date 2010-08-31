from AZutilities import dataUtilities
import unittest
import os,sys
import time

import orange
from trainingMethods import AZorngPLS
#from trainingMethods import AZorngSVM
from trainingMethods import AZorngCvSVM
from trainingMethods import AZorngANN
from trainingMethods import AZorngRF
from trainingMethods import AZorngCvANN

from AZutilities import evalUtilities
from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import AZLearnersParamsConfig
import commands
from AZutilities import paramOptUtilities


class optimizerTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        """Creates the training and testing data set attributes. """
        # The original templateProfile created at installation is needed in order to call the runScript of optimizer from 
        #appspack with proper environment.
        self.assert_(os.path.isfile(os.path.join(os.environ["AZORANGEHOME"], "templateProfile")),"Missing file: "+os.path.join(os.environ["AZORANGEHOME"], "templateProfile"))
        self.contTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")
        self.contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")
        self.discTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")
        self.discTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")

        # Read in the data
        self.discTrain = dataUtilities.DataTable(self.discTrainDataPath)
        self.discTest = dataUtilities.DataTable(self.discTestDataPath)
        self.contTrain = dataUtilities.DataTable(self.contTrainDataPath)
        self.contTest = dataUtilities.DataTable(self.contTestDataPath)


    def testRF_MPI(self):
        """
        Tests changing the default range of the optimizer.
        Use MPI versio0n of appspack
        """
        # Classification accuracy:
        ExpectedCA = [0.903] #opencv1.1: 0.90480000000000005

        optimizer = paramOptUtilities.Appspack()

        learner = AZorngRF.RFLearner()
        learnerName = "RFLearner"

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)

        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["nActVars"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        # Set the NumThreads
        pars.setParameter("NumThreads","optimize",False)
        # Change the default
        pars.setParameter("NumThreads","default","1")

        trainFile=self.discTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_RF_MPI")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        useDefaultPoint = False,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useStd = False,\
                        verbose = 0,\
                        #advancedMPIoptions = "-all-local -allcpus")  # to use this the 
                        # file "<MPICHDIR>/share/machines.LINUX must be properly configured"
                        # Alternatively, we can set machinefile=0 to us also all available cores
                        machinefile =0)

        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"
        print "Number of cores used: ",optimizer.np

        verbTunedPars = optimizer.getTunedParameters()

        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)
        
        #Check if the number of processors used are all the core available
        status,out = commands.getstatusoutput("cat /proc/cpuinfo | grep processor")
        self.assertEqual(optimizer.np, len(out.split("\n")))

        # Check if the MPI version was used
        self.assertEqual(optimizer.usedMPI, True)

        # Check the number of optimized parameters
        self.assert_(len(verbTunedPars["optParam"]) in [8,9,10])

        # Check the accuracy
        self.assert_(round(verbTunedPars["bestRes"],3) in [round(x,3) for x in ExpectedCA],"Got:" + str(verbTunedPars["bestRes"]))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=3) #  Must be > 2

        miscUtilities.removeDir(runPath)


    def testCvANN_MPI(self):
        """             
        Tests changing the default range of the optimizer.
        Use MPI versio0n of appspack
        """             
        # Classification accuracy:
        ExpectedCA = [0.846,0.864761904762] #Orange 1: [0.895238095238] #opencv1.1: [0.86619047619, 0.87619047619000001]
                        
        optimizer = paramOptUtilities.Appspack()

        learner = AZorngCvANN.CvANNLearner()
        learnerName = "CvANNLearner"
        
        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)
        
        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["maxIter","nHidden"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        trainFile=self.discTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_CvANN_MPI")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        runPath = runPath,\
                        verbose = 0,\
                        useStd = False,\
                        #advancedMPIoptions = "-all-local -allcpus")  # to use this the 
                        # file "<MPICHDIR>/share/machines.LINUX must be properly configured"
                        np = 4,\
                        machinefile =["localhost:2","localhost:2"])

        verbTunedPars = optimizer.getTunedParameters()

        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "Best result index from intRes file:",verbTunedPars["ResIdx"]
        print "Optimizer runPath:",runPath
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)

        # Check if the MPI version was used
        self.assertEqual(optimizer.usedMPI, True)


        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 14)

        # Check the accuracy
        self.assert_(round(verbTunedPars["bestRes"],3) in [round(x,3) for x in ExpectedCA],"Got:" + str(verbTunedPars["bestRes"]))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=5) # Must be > 2

        miscUtilities.removeDir(runPath)




    def testANN_MPI(self):
        """
        Tests changing the default range of the optimizer.
        Use MPI versio0n of appspack
        """
        # Classification accuracy:
        ExpectedCA = [0.846]


        optimizer = paramOptUtilities.Appspack()

        learner = AZorngANN.ANNLearner()
        learnerName = "ANNLearner"

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)

        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["nEpochs","nHidden"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        trainFile=self.discTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_ANN_MPI")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        runPath = runPath,\
                        verbose = 0,\
                        useStd = False,\
                        #advancedMPIoptions = "-all-local -allcpus")  # to use this the 
                        # file "<MPICHDIR>/share/machines.LINUX must be properly configured"
                        np = 4,\
                        machinefile =["localhost:2","localhost:2"])

        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)

        # Check if the MPI version was used
        self.assertEqual(optimizer.usedMPI, True)

        verbTunedPars = optimizer.getTunedParameters()

        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 6)

        # Check the accuracy
        self.assert_(round(verbTunedPars["bestRes"],3) in [round(x,3) for x in ExpectedCA],"Got:" + str(verbTunedPars["bestRes"]))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=14) # Must be > 2

        miscUtilities.removeDir(runPath)

    def testCvSVM_MPI(self):
        """
        Tests changing the default range of the optimizer.
        Use MPI versio0n of appspack
        """
        # Classification accuracy:
        ExpectedCA = [0.827] #New at orange2.0
                                           
        optimizer = paramOptUtilities.Appspack()

        learner = AZorngCvSVM.CvSVMLearner()
        learnerName = "CvSVMLearner"

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)

        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["C", "gamma"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        # Change the range
        pars.setParameter("C","range",miscUtilities.power2Range(-5,2,1))
        pars.setParameter("priors","default",{"POS":2, "NEG":4})

        trainFile=self.discTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_CvSVM_MPI")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        runPath = runPath,\
                        verbose = 0,\
                        useStd = False,\
                        #advancedMPIoptions = "-all-local -allcpus")  # to use this the 
                        # file "<MPICHDIR>/share/machines.LINUX must be properly configured"
                        np = 4,\
                        machinefile = os.path.realpath(os.path.join(os.environ["AZORANGEHOME"], "tests/source/APPS_machines")))

        verbTunedPars = optimizer.getTunedParameters()

        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "Best result index from intRes file:",verbTunedPars["ResIdx"]
        print "Optimizer runPath:",runPath
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)

        # Check if the MPI version was used
        self.assertEqual(optimizer.usedMPI, True)

        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 12)

        # Check the accuracy
        self.assert_(round(verbTunedPars["bestRes"],3) in [round(x,3) for x in ExpectedCA],"Got:"+str(verbTunedPars["bestRes"]))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=12) # (orig: 14)  Must be > 2


        #Check Priors
        self.assertEqual( tunedPars[1]["priors"],"None")
        learner.priors = {'NEG':4,'POS':2}
        classifier = learner(self.discTest)
        classifier.write(os.path.join(runPath,"CvSVMModel"))
        file = open(os.path.join(runPath,"CvSVMModel/model.svm"),"r")
        lines = file.readlines()
        file.close()
        priors = [round(x,2) for x in eval((lines[18].strip()).replace("data:",""))]
        self.assertEqual(len(priors),2)
        self.assertEqual(priors[self.discTest.domain.classVar.values.index("POS")],2.0*float(tunedPars[1]["C"]))
        self.assertEqual(priors[self.discTest.domain.classVar.values.index("NEG")],4.0*float(tunedPars[1]["C"]))
        miscUtilities.removeDir(runPath)


    def testSVM_MPI(self):
        """
        Tests changing the default range of the optimizer.
        Use MPI versio0n of appspack
        """
        # Classification accuracy:
        ExpectedCA = [0.847, 0.838] #New at orange2.0

        optimizer = paramOptUtilities.Appspack()

        learner = AZorngCvSVM.CvSVMLearner()
        learnerName = "CvSVMLearner"
   
        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)
        
        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["C", "gamma"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        # Change the range
        pars.setParameter("C","range",miscUtilities.power2Range(-5,2,1))

        trainFile=self.discTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_SVM_MPI")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        runPath = runPath,\
                        verbose = 0,\
                        useStd = False,\
                        #advancedMPIoptions = "-all-local -allcpus")  # to use this the 
                        # file "<MPICHDIR>/share/machines.LINUX must be properly configured"
                        np = 4,\
                        machinefile = os.path.realpath(os.path.join(os.environ["AZORANGEHOME"], "tests/source/APPS_machines")))

        verbTunedPars = optimizer.getTunedParameters()

        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)

        # Check if the MPI version was used
        self.assertEqual(optimizer.usedMPI, True)

        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 12)

        # Check the accuracy
        self.assert_(round(verbTunedPars["bestRes"],3) in  [round(x,3) for x in ExpectedCA],"Got:"+str(verbTunedPars["bestRes"]))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=12) # (orig: 14)  Must be > 2

        miscUtilities.removeDir(runPath)


    def testSVM_MPI_2(self):
        ###################################################################
        #       Test other way of setting appspack
        ###################################################################
         # Classification accuracy:
        ExpectedCA = 0.847 #Orange 1:   0.837619047619

        optimizer = paramOptUtilities.Appspack()

        learner = AZorngCvSVM.CvSVMLearner()
        learnerName = "CvSVMLearner"

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)

        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["C", "gamma"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        # Change the range
        pars.setParameter("C","range",miscUtilities.power2Range(-5,2,1))

        trainFile=self.discTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_SVM_MPI_2")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        #[<minimum of objective function found>, <optimized parameters>]
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useParameters = pars.getParametersDict(),\
                        verbose = 0,
                        useStd = False,\
                        advancedMPIoptions = None,
                        np = 4,
                        machinefile = ["localhost:2","localhost:2"])
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        self.assertEqual(learner.optimized,True)
        # Check if the MPI version was used
        self.assertEqual(optimizer.usedMPI, True)
        self.assertEqual(round(tunedPars[0],3),round(ExpectedCA,3))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=12) # (orig 14) Must be > 2
        #print runPath
        miscUtilities.removeDir(runPath)


    def testSVM_MPI_3(self):
        ###################################################################
        #       Test other way of setting appspack
        ###################################################################
        # Classification accuracy:
        ExpectedCA = 0.847 #orange1: 0.837619047619

        optimizer = paramOptUtilities.Appspack()

        learner = AZorngCvSVM.CvSVMLearner()
        learnerName = "CvSVMLearner"

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)

        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["C", "gamma"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        # Change the range
        pars.setParameter("C","range",miscUtilities.power2Range(-5,2,1))

        trainFile=self.discTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_SVM_MPI_3")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        #[<minimum of objective function found>, <optimized parameters>]
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useParameters = pars.getParametersDict(),\
                        verbose = 0,\
                        useStd = False,\
                        advancedMPIoptions = "-v -np 4",\
                        machinefile = ["localhost:2","localhost:2"])
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        # Check if the MPI version was used
        self.assertEqual(learner.optimized,True)
        # Check if the MPI version was used
        self.assertEqual(optimizer.usedMPI, True)
        self.assertEqual(round(tunedPars[0],3),round(ExpectedCA,3))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=12) # (orig 14) Must be > 2
        #print runPath
        miscUtilities.removeDir(runPath)


    def test_PLS_MPI(self):
        """Test of optimizer with advanced configuration
        """
        # Classification accuracy:
        ExpectedCA = [0.852, 0.865]
        ExpectedCAwithTest = [0.837142857143, 0.865238095238, 0.85619047619, 0.875] #New at orange2.0 

        #Create  the appspack instance
        opt=paramOptUtilities.Appspack()
        #Learner to be optimized
        learner=AZorngPLS.PLSLearner()
        #dataset to use in the parameters optimization (Discrete class in this example)
        dataSet=self.discTrainDataPath
        # Define the objective function. This requires:
        #    defining the extreme to find (the min or max): findMin=True or findMin=False
        fMin=False
        #    defining the method for evaluation (must be a method that accepts as input an orngTest.ExperimentResults): 
        #       evaluateMethod="AZutilities.evalUtilities.CA"
        evalM="AZutilities.evalUtilities.CA"

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_PLS_MPI")

        # Load the optimization parameters from the default configuration (AZLearnersParamsConfig.py)
        parameters = AZLearnersParamsConfig.API("PLSLearner")
        parameters.setParameter("method","default",'pls1')

        # change the optimization parameters
        parameters.setParameter("method","default",'pls1')      #   make the method fixed (do not optimize) to be pls1
        parameters.setParameter("method","optimize",False)         
        parameters.setParameter("method","rangeType","values")  #   assure that the keyword for the values range type is 
                                                                #set correctly for values instead of interval

        parameters.setParameter("k","range",[1 , 3 , 5 , 6 , 10])      #   make the method fixed (do not optimize) to be pls1
        parameters.setParameter("k","optimize",True)
        parameters.setParameter("k","rangeType","values")  #   assure that the keyword for the values range type is
                                                           #set correctly for values instead of interval

        # Run the appspack which will configure the input learner and aditionaly return 
        #[<minimum of objective function found>, <optimized parameters>]
        tunedPars = opt(learner=learner,\
                        dataSet=dataSet,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useStd = False,\
                        useParameters = parameters.getParametersDict(),\
                                                        #  The 'useParameters' is mandatory, even placing a file with the new configurations in the
                                                        # running directory, that we pass to the optimizer the correct parameters to use.
                                                        #  The parameters placed on the running directory are for appspack usage, and the 
                                                        # optimizer needs to know what parameters appspack will use, otherwise, it will 
                                                        # load the default ones 
                        verbose = 0,
                        np = 6,
                        machinefile = ["localhost:2","localhost:2"])
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        
        self.assertEqual(learner.optimized,True)
        # Check if the MPI version was used
        self.assertEqual(opt.usedMPI, True)
        self.assert_(round(tunedPars[0],3) in [round(x,3) for x in ExpectedCAwithTest],"Got:" + str(tunedPars[0]))

        #The learner is now with its optimized parameters already set, so we can now make a classifier out of it
        classifier = learner(self.discTrain)
        CA = evalUtilities.getClassificationAccuracy(self.discTest,classifier)
        self.assert_(round(CA,3) in [round(x,3) for x in ExpectedCA])
        resData1 = dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt"))
        self.assert_(len(resData1)>=4) # (orig: 5)  Must be > 2
        #print runPath
        miscUtilities.removeDir(runPath)

    def test_PLS_MPI_2(self):
        ###################################################################
        #       Test other way of setting appspack
        ###################################################################i
        # Classification accuracy:
        ExpectedCA = [0.851851851852, 0.865]
        ExpectedCAwithTest = [0.865238095238, 0.884285714286, 0.85619047619, 0.837] #New at orange2.0

        #Create  the appspack instance
        opt=paramOptUtilities.Appspack()
        #Learner to be optimized
        learner=AZorngPLS.PLSLearner()
        #dataset to use in the parameters optimization (Discrete class in this example)
        dataSet=self.discTrainDataPath
        # Define the objective function. This requires:
        #    defining the extreme to find (the min or max): findMin=True or findMin=False
        fMin=False
        #    defining the method for evaluation (must be a method that accepts as input an orngTest.ExperimentResults): 
        #       evaluateMethod="AZutilities.evalUtilities.CA"
        evalM="AZutilities.evalUtilities.CA"

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_PLS_MPI_2")

        # Load the optimization parameters from the default configuration (AZLearnersParamsConfig.py)
        parameters = AZLearnersParamsConfig.API("PLSLearner")
        parameters.setParameter("method","default",'pls1')

        # change the optimization parameters
        parameters.setParameter("method","default",'pls1')      #   make the method fixed (do not optimize) to be pls1
        parameters.setParameter("method","optimize",False)
        parameters.setParameter("method","rangeType","values")  #   assure that the keyword for the values range type is 
                                                                #set correctly for values instead of interval

        parameters.setParameter("k","range",[1 , 3 , 5 , 6 , 10])      #   make the method fixed (do not optimize) to be pls1
        parameters.setParameter("k","optimize",True)
        parameters.setParameter("k","rangeType","values")  #   assure that the keyword for the values range type is
                                                           #set correctly for values instead of interval
        #[<minimum of objective function found>, <optimized parameters>]
        tunedPars = opt(learner=learner,\
                        dataSet=dataSet,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useParameters = parameters.getParametersDict(),\
                        verbose = 0,
                        useStd = False,\
                        advancedMPIoptions = None,
                        np = 4,
                        machinefile = ["localhost:2","localhost:2"])
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"


        self.assertEqual(learner.optimized,True)
        # Check if the MPI version was used
        self.assertEqual(opt.usedMPI, True)
        self.assert_(round(tunedPars[0],3) in [round(x,3) for x in ExpectedCAwithTest],"Got:" + str(tunedPars[0]))

        #The learner is now with its optimized parameters already set, so we can now make a classifier out of it
        classifier = learner(self.discTrain)
        CA = evalUtilities.getClassificationAccuracy(self.discTest,classifier)
        self.assert_(round(CA,3)  in [round(x,3) for x in ExpectedCA])
        resData2 = dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt"))
        self.assert_(len(resData2)>=4) # (orig 5) Must be > 2
        #print runPath
        miscUtilities.removeDir(runPath)
 


if __name__ == "__main__":
        ##scPA   Test first if the service is available from this machine. If not, report it and make no tests
        if not os.path.isfile(os.path.join(os.environ["AZORANGEHOME"],"orangeDependencies/bin/appspack_mpi")):
            print "SKIP (No appspack_mpi installed)"
            sys.exit(1)

        #Test SSH connection to localhost
#        import AutoConfig
#        report = AutoConfig.config()()
#        if report != "":
#            print "AutoConfig report:\n"+report
#            print "SKIP (Unable to configure properly interfaceless SSH)"     
#            sys.exit(1)
         
        ##ecPA
        suite = unittest.TestLoader().loadTestsFromTestCase(optimizerTest)
        unittest.TextTestRunner(verbosity=2).run(suite)


