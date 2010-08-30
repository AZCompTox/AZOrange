from AZutilities import dataUtilities
import unittest
import os
import time

import orange
from trainingMethods import AZorngPLS
#from trainingMethods import AZorngSVM
#from trainingMethods import AZorngANN
from trainingMethods import AZorngRF
from trainingMethods import AZorngCvANN
from trainingMethods import AZorngCvSVM

from AZutilities import evalUtilities
from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import AZLearnersParamsConfig

from AZutilities import paramOptUtilities


class optimizerTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        """Creates the training and testing data set attributes. """
        # The original templateProfile created at installation is needed in order to call the runScript of optimizer from 
        #appspack with proper environment.
        self.assert_(os.path.isfile(os.path.join(os.environ["AZORANGEHOME"], "templateProfile")),"Missing file: "+os.path.join(os.environ["AZORANGEHOME"], "templateProfile"))

        self.contTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")
        self.contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")
        self.discTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyTrain.tab")
        self.discTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyTest.tab")

        # Read in the data
        self.discTrain = dataUtilities.DataTable(self.discTrainDataPath)
        self.discTest = dataUtilities.DataTable(self.discTestDataPath)
        self.contTrain = dataUtilities.DataTable(self.contTrainDataPath)
        self.contTest = dataUtilities.DataTable(self.contTestDataPath)

        

    def test_GridSearch(self):
        """
        Test GridSearch Module
        """
        #Create  the appspack instance
        opt=paramOptUtilities.Appspack()
        #Learner to be optimized
        learner=AZorngRF.RFLearner()
        #dataset to use in the parameters optimization (Discrete class in this example)
        dataSet=self.discTestDataPath
        # Define the objective function. This requires:
        #    defining the extreme to find (the min or max): findMin=True or findMin=False
        fMin=False
        #    defining the method for evaluation (must be a method that accepts as input an orngTest.ExperimentResults): 
        #       evaluateMethod="AZutilities.evalUtilities.CA"
        evalM="AZutilities.evalUtilities.CA"

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="RFTest")

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API("RFLearner")
        # Set the parameters in parameterList to be optimized
        pars.setParameter("NumThreads","optimize",False)
        # Change the default
        pars.setParameter("NumThreads","default","1")


        # Run the appspack which will configure the input learner and aditionaly return
        #[<minimum of objective function found>, <optimized parameters>]
        tunedPars = opt(learner=learner,\
                        dataSet=dataSet,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useGridSearchFirst = True,\
                        gridSearchInnerPoints = 3,\
                        useDefaultPoint = False,\
                        useStd = False,\
                        useParameters = pars.getParametersDict(),\
                        verbose = 0)
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"
        print "CheckSum:", round(sum(opt.GSRes["results"]),2)
        print "Number of results: ",len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))
        print "Running Path:",runPath
        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)

        # Check the accuracy
        self.assertEqual(round(tunedPars[0],2), round(0.766700000000000,2)) #0.853333333333,4

        #Check if the number of results remain equal
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=5)

        #Check that all points were evaluated
        self.assert_(opt.GSRes["nFailedPoints"]==0)
        self.assert_(opt.GSRes["nPoints"]==3)
        #CheckSum to assure results are the same
        self.assert_(round(sum(opt.GSRes["results"]),2) == -2.09)  # -2.44

        #Check if the best result was not the one with numThreads different of 1 since that way we can get 
        #different results among runs
        self.assertEqual(int(tunedPars[1]["NumThreads"]),1)

        miscUtilities.removeDir(runPath)
 



    def no_test_SVMClassificationRange(self): #Disabled since we do not have the AZSVMLearner anymore
        """
        SVM - Tests changing the default range of the optimizer.
        """
        optimizer = paramOptUtilities.Appspack()

        learner = AZorngSVM.AZSVMLearner()
        learnerName = "AZSVMLearner"
   
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
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        useStd = False,\
                        runPath = runPath,\
                        verbose = 0)
        self.assertEqual(optimizer.usedMPI,False)

        verbTunedPars = optimizer.getTunedParameters()

        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        
        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)

        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 12)

        # Check the accuracy
        self.assertEqual(round(verbTunedPars["bestRes"],2), round(0.837619047619,2))

        miscUtilities.removeDir(runPath)

    def no_test_ANNregressionDefault(self): #Disabled since we do not have the AZZLearner anymore
        """
        ANN - Tests changing the default range of the optimizer.
        """
        optimizer = paramOptUtilities.Appspack()

        learner = AZorngANN.ANNLearner()
        learnerName = "ANNLearner"
   
        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API(learnerName)
        
        # Set all parameters to not be optimized
        pars.setOptimizeAllParameters(False)

        parameterList = ["nHidden"]
        # Set the parameters in parameterList to be optimized
        for parameter in parameterList:
            pars.setParameter(parameter,"optimize",True)

        # Change the ANN minimization algorithm, still not optimized
        pars.setParameter("optAlg","default","RProp")

        trainFile=self.contTrainDataPath

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest")

        evalM = "AZutilities.evalUtilities.RMSE"
        fMin = True

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        runPath = runPath,\
                        useStd = False,\
                        verbose = 0)
        self.assertEqual(optimizer.usedMPI,False)
        verbTunedPars = optimizer.getTunedParameters()

        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"
        

        # Check that the learner was optimized
        self.assertEqual(learner.optimized,True)

        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 6)

        # Check the accuracy
        self.assertEqual(round(verbTunedPars["bestRes"],2), round(0.9142000000,2))

        miscUtilities.removeDir(runPath)

    def test_PLSAdvanced_Usage(self):
        """PLS - Test of optimizer with advanced configuration
        """
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
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest")

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
                        verbose = 0)
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"
        self.assertEqual(opt.usedMPI,False)
        self.assertEqual(learner.optimized,True)
        self.assertEqual(round(tunedPars[0],2),round(0.86523799999999995,2)) #0.87571399999999,6

        #The learner is now with its optimized parameters already set, so we can now make a classifier out of it
        classifier = learner(self.discTrain)
        CA = evalUtilities.getClassificationAccuracy(self.discTest,classifier)
        self.assertEqual(round(CA,2),round(0.851851851852,2))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=5) # Must be > 2
        miscUtilities.removeDir(runPath)

    def test_PLS_Classification(self):
        """PLS - Test of optimizer with discrete class data
        """
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
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest")

        # Run the appspack which will configure the input learner and aditionaly return 
        #[<minimum of objective function found>, <optimized parameters>]
        tunedPars = opt(learner=learner,\
                        dataSet=dataSet,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useStd = False,\
                        verbose = 0)
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"
        self.assertEqual(opt.usedMPI,False)
        self.assertEqual(learner.optimized,True)
        self.assertEqual(round(tunedPars[0],2),round(0.87523799999999996,2))


        #The learner is now with its optimized parameters already set, so we can now make a classifier out of it
        classifier = learner(self.discTrain)
        CA = evalUtilities.getClassificationAccuracy(self.discTest,classifier)
        self.assertEqual(round(CA,2),round(0.851851851852,2))

        miscUtilities.removeDir(runPath)

    def test_PLS_Regression(self):
        """PLS - Test of optimizer with continuous class data        
        """        
        #Create  the appspack instance
        opt=paramOptUtilities.Appspack()
        #Learner to be optimized
        learner=AZorngPLS.PLSLearner()
        #dataset to use in the parameters optimization 
        dataSet=self.contTrainDataPath
        # Define the objective function. This requires:
        #    defining the extreme to find (the min or max): findMin=True or findMin=False
        fMin=True
        #    defining the method for evaluation (must be a method that accepts as input an orngTest.ExperimentResults): 
        #       evaluateMethod="AZutilities.evalUtilities.R2"
        evalM="AZutilities.evalUtilities.RMSE"
        
        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest")
        
        # Run the appspack which will configure the input learner and aditionaly return 
        #[<minimum of objective function found>, <optimized parameters>]        
        tunedPars = opt(learner=learner,\
                        dataSet=dataSet,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useStd = False,\
                        verbose = 0)
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"
        self.assertEqual(opt.usedMPI,False)
        self.assertEqual(learner.optimized,True)
        self.assertEqual(round(tunedPars[0],2),round(0.858060000000,2))
        #The learner is now with its optimized parameters already set, so we can now make a classifier out of it
        classifier = learner(self.contTrain)
        RMSE = evalUtilities.getRMSE(self.contTest,classifier)
        self.assertEqual(round(RMSE,2),round(0.656979500000,2))

        miscUtilities.removeDir(runPath)

    def test_RFClassification(self):
        """RF - Test of optimizer with discrete class data
        """
        #Create  the appspack instance
        opt=paramOptUtilities.Appspack()
        #Learner to be optimized
        learner=AZorngRF.RFLearner()
        #dataset to use in the parameters optimization (Discrete class in this example)
        dataSet=self.discTrainDataPath
        # Define the objective function. This requires:
        #    defining the extreme to find (the min or max): findMin=True or findMin=False
        fMin=False
        #    defining the method for evaluation (must be a method that accepts as input an orngTest.ExperimentResults): 
        #       evaluateMethod="AZutilities.evalUtilities.CA"
        evalM="AZutilities.evalUtilities.CA"

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="RFTest")

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API("RFLearner")
        # Set the parameters in parameterList to be optimized
        pars.setParameter("NumThreads","optimize",False)
        # Change the default
        pars.setParameter("NumThreads","default","1")


        # Run the appspack which will configure the input learner and aditionaly return 
        #[<minimum of objective function found>, <optimized parameters>]
        tunedPars = opt(learner=learner,\
                        dataSet=dataSet,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useDefaultPoint = False,\
                        useStd = False,\
                        useParameters = pars.getParametersDict(),\
                        verbose = 0)
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"
        print "Number of optimization steps: ",len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))
        print "Number of Threads used: ",learner.NumThreads
        #The learner is now with its optimized parameters already set, so we can now make a classifier out of it
        learner.NumThreads = 1 
        classifier = learner(self.discTrain)
        CA = evalUtilities.getClassificationAccuracy(self.discTest,classifier)
        print "CA of optimized Learner: ",CA 

        self.assertEqual(opt.usedMPI,False)

        self.assertEqual(learner.optimized,True)
        self.assertEqual(round(tunedPars[0],2),round(0.90333333333300003,2)) #  0.92380952381000003      opencv1.1: 0.92380952381000003 (0.93333333333299995)

        self.assertEqual(round(CA,2),round(0.81481499999999996,2)) #opencv 1.1: 0.77777799999999997
        #Check if the best result was not the one with numThreads different of 1 since that way we can get 
        #different results among runs
        self.assertEqual(int(tunedPars[1]["NumThreads"]),1)

        miscUtilities.removeDir(runPath)


    def test_RFRegression(self):
        """RF - Test of optimizer with continuous class data        
        """
        #Create  the appspack instance
        opt=paramOptUtilities.Appspack()
        #Learner to be optimized
        learner=AZorngRF.RFLearner()
        #dataset to use in the parameters optimization (Discrete class in this example)
        dataSet=self.contTrainDataPath
        # Define the objective function. This requires:
        #    defining the extreme to find (the min or max): findMin=True or findMin=False
        fMin=True
        #    defining the method for evaluation (must be a method that accepts as input an orngTest.ExperimentResults): 
        #       evaluateMethod="AZutilities.evalUtilities.R2"
        evalM="AZutilities.evalUtilities.RMSE"

        # Create an interface for setting optimizer parameters
        pars = AZLearnersParamsConfig.API("RFLearner")
        # Set the parameters in parameterList to be optimized
        pars.setParameter("NumThreads","optimize",False)
        # Change the default
        pars.setParameter("NumThreads","default","1")

        # Create a directory for running the appspack (if not defined it will use the present working directory)
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest")

        # Run the appspack which will configure the input learner and aditionaly return 
        #[<minimum of objective function found>, <optimized parameters>]        
        tunedPars = opt(learner=learner,\
                        dataSet=dataSet,\
                        evaluateMethod = evalM,\
                        findMin=fMin,\
                        runPath = runPath,\
                        useStd = False,\
                        useParameters = pars.getParametersDict(),\
                        verbose = 0)
        print "Returned: ", tunedPars
        print "====================== optimization Done ==========================="
        print "Learner optimized flag = ", learner.optimized
        print "Tuned parameters = ", tunedPars[1]
        print "Best optimization result = ", tunedPars[0]
        print "check the file intRes.txt to see the intermediate results of optimizer!"

        self.assertEqual(opt.usedMPI,False)
        self.assertEqual(learner.optimized,True)
        self.assertEqual(round(tunedPars[0],2),round(0.72,2))

        #The learner is now with its optimized parameters already set, so we can now make a classifier out of it
        classifier = learner(self.contTrain)
        RMSE = evalUtilities.getRMSE(self.contTest,classifier)
        self.assertEqual(round(RMSE,2),round(0.3387,2)) #opencv1.1: 0.33829999999999999

        #Check if the best result was not the one with numThreads different of 1 since that way we can get 
        #different results among runs
        self.assertEqual(int(tunedPars[1]["NumThreads"]),1)

        miscUtilities.removeDir(runPath)


    def testCvANN(self):
        """             
        Tests changing the default range of the optimizer.
        """             
        # Classification accuracy:
        ExpectedCA = [0.846]  #0.864761904762  #0.895238095238 #opencv1.1: 0.87619047619
                        
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
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_CvANN")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        useStd = False,\
                        runPath = runPath,\
                        verbose = 0)

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

        # Check if the MPI version was not used
        self.assertEqual(optimizer.usedMPI, False)

        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 14)

        # Check the accuracy
        nOptPoints = len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))
        self.assert_(nOptPoints > 5,"N. of optimization points:" + str(nOptPoints)) # Must be > 2
        self.assert_(round(verbTunedPars["bestRes"],3) in [round(x,3) for x in ExpectedCA],"Actual result:"+str(verbTunedPars["bestRes"]))


        miscUtilities.removeDir(runPath)

    def testCvSVM(self):
        """
        Tests changing the default range of the optimizer.
        """
        # Classification accuracy:
        ExpectedCA = [0.827142857143] # New at orange2.0

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
        runPath = miscUtilities.createScratchDir(desc="ParamOptTest_CvSVM")
        evalM = "AZutilities.evalUtilities.CA"
        fMin = False

        # Calculate the optimal parameters. This can take a long period of time!
        tunedPars = optimizer(learner=learner,\
                        dataSet=trainFile,\
                        evaluateMethod = evalM,\
                        useParameters = pars.getParametersDict(),\
                        findMin=fMin,\
                        useStd = False,\
                        runPath = runPath,\
                        verbose = 0)

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

        # Check if the MPI version was not used
        self.assertEqual(optimizer.usedMPI, False)

        # Check the number of optimized parameters
        self.assertEqual(len(verbTunedPars["optParam"]), 12)

        # Check the accuracy
        self.assert_(round(verbTunedPars["bestRes"],2) in [round(x,2) for x in ExpectedCA],"Got:"+str(verbTunedPars["bestRes"]))
        self.assert_(len(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt")))>=12) # (orig: 14)  Must be > 2


        #Check Priors
        self.assertEqual(dataUtilities.DataTable(os.path.join(runPath,"optimizationLog.txt"))[1]["priors"].value,"{'NEG':4,'POS':2}")
        self.assertEqual( tunedPars[1]["priors"],"None") 

        #Set the priors since it could be choosing the first row as the best, which would be the default values, without the priors
        learner.priors = {"POS":2, "NEG":4}

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



if __name__ == "__main__":
        suite = unittest.TestLoader().loadTestsFromTestCase(optimizerTest)
        unittest.TextTestRunner(verbosity=2).run(suite)


