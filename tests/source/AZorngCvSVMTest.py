from AZutilities import dataUtilities
import unittest
import os
import time

import orange
from trainingMethods import AZorngCvSVM
from AZutilities import miscUtilities
from AZutilities import evalUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngTest

class SVMClassifierTest(AZorngTestUtil.AZorngTestUtil):
    """Class to test the model IO of the svm and svm Jakulin models."""

    def setUp(self):
        """Creates the training and testing data set attributes. """
        self.dataPathD = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris.tab")
        self.dataPathC = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Test.tab")


        # Read in the data
        self.inDataD = dataUtilities.DataTable(self.dataPathD)
        self.inDataC = dataUtilities.DataTable(self.dataPathC)


        # Full path to saved svm model
        global scratchdir
        self.modelPath = os.path.join(scratchdir,"model.svm")

        """Other datasets..."""
        contDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Imp_Test.tab")
        SVMregDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")
        contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Imp_Train.tab")
        dataNoMetaTrainPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train.tab")
        missingTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train_missing.tab")

        #These 2 datasets are equal apart from the meta atribute
        dataNoMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_SmallTest.tab")
        dataWMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_W_metas_SmallTest.tab")

        # Read in the data
        #IrisData
        trainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisTrain.tab")
        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisTest.tab")
        self.train_data =  dataUtilities.DataTable(trainDataPath)
        self.test_data =  dataUtilities.DataTable(testDataPath)

        missingInData = dataUtilities.DataTable(missingTestDataPath)
        contTrainData = dataUtilities.DataTable(contTrainDataPath)
        self.regTrainData = dataUtilities.DataTable(SVMregDataPath)
        contData = dataUtilities.DataTable(contDataPath)
        self.NoMetaTrain = dataUtilities.DataTable(dataNoMetaTrainPath)
        self.NoMetaTest = dataUtilities.DataTable(dataNoMetaTestPath)
        self.WMetaTest = dataUtilities.DataTable(dataWMetaTestPath)

        self.missingTrain = missingInData
        self.missingTest = missingInData
        self.contTrain = contTrainData
        self.contTest = contData
        
        #Data for domain fix handling
        
        badVarTypePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarType.tab")
        badVarNamePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarName.tab")
        badVarOrderPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarOrder.tab")
        badVarCountPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarCount.tab")
        # Read in the data
        self.noBadDataTrain = self.NoMetaTrain
        self.noBadDataTest = self.NoMetaTest
        self.badVarTypeData = dataUtilities.DataTable(badVarTypePath)
        self.badVarNameData = dataUtilities.DataTable(badVarNamePath)
        self.badVarOrderData = dataUtilities.DataTable(badVarOrderPath)
        self.badVarCountData = dataUtilities.DataTable(badVarCountPath)   #One less example 

    def test__PredictionVarImportance(self):
        """Test the variable importance for single predictions"""
        def TopVarImportanceTest(data, expectNone = False):
            resA = []
            resB = []
            learner = AZorngCvSVM.CvSVMLearner(gamma=1.0,svm_type=103,
                          C=1, coef0=0, degree=3, epsR=0.001, kernel_type=2,
                          nu=0.5, p=0.1, probability=0, shrinking=1)
            CvSVM = learner(data)

            for ex in data:
                resA.append(CvSVM.getTopImportantVars(ex,1))
                
            scratchdir = miscUtilities.createScratchDir(desc="TopVarImportanceTest")
            modelPath = os.path.join(scratchdir,"CvSVNModel")
            CvSVM.write(modelPath)
            LoadedCvSVM = AZorngCvSVM.CvSVMread(modelPath)
            miscUtilities.removeDir(scratchdir) 

            for ex in data:
                resB.append(LoadedCvSVM.getTopImportantVars(ex,1))
            if expectNone:
               return resA == resB == [None]*len(data)
            else:
                return resA == resB and None not in resA and resA.count(resA[0]) != len(resA)


        self.assertEqual(TopVarImportanceTest(self.inDataD, True), True, "Failed testing with uncompatible data (class with 3 values)")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisCont.tab"))
        self.assertEqual(TopVarImportanceTest(data), True,"Failed testing with Regression data.")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab"))
        self.assertEqual(TopVarImportanceTest(data), True, "Failed testing with Binary Classifier")





    def test_DFV(self):
        """ Test the Decision Function Value Return"""
        CvSVM = AZorngCvSVM.CvSVMLearner(self.NoMetaTrain)
        #Testsing with return of DFV
        RDFV = True
        for ex in self.NoMetaTest:
            predictedClass = CvSVM(ex)
            a = CvSVM(ex,returnDFV = RDFV)
            b = CvSVM(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = CvSVM(ex,resultType = orange.GetBoth,returnDFV = RDFV)
            #All must return tuples    
            self.assert_(type(a)==type(b)==type(c)==tuple)
            # Second element of the tupple must be the DFV
            self.assert_(type(a[1])==type(b[1])==type(c[1])==float)
            self.assert_(a[1]==b[1]==c[1])
            # check if if the class can be always predicted based on the DFV
            # Positive values will correspond always to the fisrt element of the class variable
            #   and negative values to the second element of the class variabel
            if a[1] > 0:
                guessedClass = ex.domain.classVar[0]
            else:
                guessedClass = ex.domain.classVar[1]
            self.assertEqual(predictedClass,guessedClass)
            #asking for GetValue
            self.assert_(type(a[0])==orange.Value)
            #asking for GetProbabilities
            self.assert_(type(b[0])==orange.DiscDistribution)
            #asking for GetBoth...
            self.assert_(type(c[0])==tuple)
            # ... where first element is the orange value...
            self.assert_(type(c[0][0])==orange.Value)
            # ... and second element is the distribution (so called probabilities)
            self.assert_(type(c[0][1])==orange.DiscDistribution)
            # Although, SVM does always return fake probabilities
            self.assertEqual(CvSVM.isRealProb(),False)
        expectedExtremes = {'max': 1.0, 'min': -1.11053}
        self.assertEqual([round(x,5) for x in CvSVM.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(CvSVM.nPredictions,4*len(self.NoMetaTest))

        #Testsing without return of DFV
        RDFV = False
        for ex in self.NoMetaTest:
            a = CvSVM(ex,returnDFV = RDFV)
            b = CvSVM(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = CvSVM(ex,resultType = orange.GetBoth,returnDFV = RDFV)
            #asking for GetValue
            self.assert_(type(a)==orange.Value)
            #asking for GetProbabilities
            self.assert_(type(b)==orange.DiscDistribution)
            #asking for GetBoth...
            self.assert_(type(c)==tuple)
            # ... where first element is the orange value...
            self.assert_(type(c[0])==orange.Value)
            # ... and second element is the distribution (so called probabilities)
            self.assert_(type(c[1])==orange.DiscDistribution)
            # Although, SVM does always return fake probabilities
            self.assertEqual(CvSVM.isRealProb(),False)

        self.assertEqual([round(x,5) for x in CvSVM.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(CvSVM.nPredictions,(4+3)*len(self.NoMetaTest))




    def test_SVMreg(self):
        CvSVMmodel = AZorngCvSVM.CvSVMLearner(self.regTrainData)
        predList = [5.482803, 4.889269,5.188474,5.528782,6.224637,4.679743,2.062022,7.878900,5.603292,5.905775,] # Ver. 0.3 
        for idx, ex in enumerate(self.regTrainData[0:10]):
            self.assertEqual(round(CvSVMmodel(ex),4), round(predList[idx], 4))


    def test_Priors(self):
        """Test to assure that priors are set correcly."""
        # Create a CvSVM model
        CvSVMlearner = AZorngCvSVM.CvSVMLearner(C=3, priors = {"Iris-versicolor":2, "Iris-virginica":4, "Iris-setosa":6})
        CvSVMmodel = CvSVMlearner(self.inDataD)

        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.inDataD, CvSVMmodel)

        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"modelPriors.CvSVM")
        CvSVMmodel.write(modelPath)

        # Read in the model
        newCvSVMmodel = AZorngCvSVM.CvSVMread(modelPath)

        # Calculate classification accuracy 
        savedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, CvSVMmodel)

        # Test that the accuracy of the two classifiers is the exact same
        self.assertEqual(Acc, savedAcc)

        #Check the priors saved in the model
        file = open(os.path.join(modelPath,"model.svm"),"r")
        lines = file.readlines()
        file.close()

        priors = [round(x,2) for x in eval((lines[18].strip()).replace("data:",""))]
        self.assertEqual(len(priors),3)
        self.assertEqual(priors[self.inDataD.domain.classVar.values.index("Iris-setosa")],18.0)
        self.assertEqual(priors[self.inDataD.domain.classVar.values.index("Iris-versicolor")],6.0)
        self.assertEqual(priors[self.inDataD.domain.classVar.values.index("Iris-virginica")],12.0)

        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


    def test_SVMD(self):

        # Train a svm
        svm = AZorngCvSVM.CvSVMLearner(self.inDataD,scaleData = False,gamma=4,C = 1,nu=0.5,p=0.1,eps=0.001, coef0=0, degree=3)
        trainedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, svm)

        self.assertEqual(round(trainedAcc,7),round(0.986666666667,7))
        # Save model 
        rc = svm.write(self.modelPath)
        self.assertEqual(rc,True)
        # Load the saved model
        loadedsvm = AZorngCvSVM.CvSVMread(self.modelPath)
        loadedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, loadedsvm)
        # Assure equal accuracy
        self.assertEqual(trainedAcc, loadedAcc)



        svmLearner = AZorngCvSVM.CvSVMLearner(scaleData = False)

        svmLearner.name = "CvSVMLearner"
        svmLearner.eps = 0.001
        svmLearner.p = 0.0
        svmLearner.nu = 0.6
        svmLearner.kernel_type = 2
        svmLearner.svm_type = 101
        svmLearner.gamma=0.0033
        svmLearner.C = 47
        svmLearner.scaleData = True
        svmLearner.scaleClass = False
        

        Res = orngTest.crossValidation([svmLearner], self.inDataD, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible)
        CA = evalUtilities.CA(Res)[0]
        self.assertEqual(round(CA,2),round(0.96666666666666667,2)) # Before in AZSVM: 0.95999999999999996

        newSVM=svmLearner(self.inDataD)
        trainedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, newSVM)
        # Save model 
        rc = newSVM.write(self.modelPath)
        self.assertEqual(rc,True)
        # Load the saved model
        loadedsvm = AZorngCvSVM.CvSVMread(self.modelPath)
        loadedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, loadedsvm)
        # Assure equal accuracy
        self.assertEqual(round(trainedAcc,7), round(0.96666669999999999,7)) #Before in AZSVM: 0.953333300000
        self.assertEqual(round(trainedAcc,1), round(loadedAcc,1))
 

    def test_SVM_Priors_D(self):
        """Test SVM with priors """
        # Train a svm
        svm = AZorngCvSVM.CvSVMLearner(self.inDataD, priors = {"Iris-setosa":0.2,"Iris-versicolor":0.3,"Iris-virginica":0.5})
        trainedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, svm)

        self.assertEqual(round(trainedAcc,7),round(0.73333329999999997,7))
        # Save model 
        rc = svm.write(self.modelPath)
        self.assertEqual(rc,True)
        # Load the saved model
        loadedsvm = AZorngCvSVM.CvSVMread(self.modelPath)
        loadedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, loadedsvm)
        # Assure equal accuracy
        self.assertEqual(trainedAcc, loadedAcc)



        svmLearner = AZorngCvSVM.CvSVMLearner(scaleData = False,priors = {"Iris-setosa":0.2,"Iris-versicolor":0.3,"Iris-virginica":0.5})

        svmLearner.name = "CvSVMLearner"
        svmLearner.shrinking = 1
        svmLearner.eps = 0.001
        svmLearner.p = 0.0
        svmLearner.nu = 0.6
        svmLearner.kernel_type = 2
        svmLearner.svm_type = 103
        svmLearner.gamma=0.0033
        svmLearner.C = 47
        svmLearner.probability = 1
        svmLearner.scaleData = True
        svmLearner.scaleClass = False
        #svmLearner.for_nomogram=1


        Res = orngTest.crossValidation([svmLearner], self.inDataD, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible)
        CA = evalUtilities.CA(Res)[0]
        self.assertEqual(round(CA,2),round(0.940000000,2))  # orange1.0: 0.93333333333333335])

        svmLearner.priors = None
        Res = orngTest.crossValidation([svmLearner], self.inDataD, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible)
        CA = evalUtilities.CA(Res)[0]
        self.assertEqual(round(CA,2),round(0.94666666666666666,2)) 


        newSVM=svmLearner(self.inDataD)
        trainedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, newSVM)
        # Save model 
        rc = newSVM.write(self.modelPath)
        self.assertEqual(rc,True)
        # Load the saved model
        loadedsvm = AZorngCvSVM.CvSVMread(self.modelPath)
        loadedAcc = evalUtilities.getClassificationAccuracy(self.inDataD, loadedsvm)
        # Assure equal accuracy
        self.assertEqual(round(trainedAcc,7), round(0.95999999999999996,7)) #Before in AZSVM: 0.953333300000
        self.assertEqual(round(trainedAcc,1), round(loadedAcc,1))

    def test_SVMC(self):

        # Train a svm
        svmL = AZorngCvSVM.CvSVMLearner(scaleData = False,svm_type = 103,gamma=0.01, C = 1,nu=0.5,p=1,eps=0.001, coef0=0, degree=3)
        svm = svmL(self.inDataC)
        trainedAcc = evalUtilities.getRMSE(self.inDataC, svm)

        self.assertEqual(round(trainedAcc,7),round(2.8525863999999999,7))# ver 0.3
        
        # Save model 
        rc = svm.write(self.modelPath)
        self.assertEqual(rc,True)
        # Load the saved model
        loadedsvm = AZorngCvSVM.CvSVMread(self.modelPath)
        loadedAcc = evalUtilities.getRMSE(self.inDataC, loadedsvm)
        # Assure equal accuracy
        self.assertEqual(trainedAcc, loadedAcc)



        svmLearner = AZorngCvSVM.CvSVMLearner(scaleData = False)

        svmLearner.name = "CvSVMLearner"
        svmLearner.eps = 0.001
        svmLearner.p = 1
        svmLearner.nu = 0.6
        svmLearner.kernel_type = 2
        svmLearner.svm_type = 103
        svmLearner.gamma=0.0033
        svmLearner.C = 47
        svmLearner.scaleData = True
        svmLearner.scaleClass = False


        Res = orngTest.crossValidation([svmLearner], self.inDataC, folds=5, strat=orange.MakeRandomIndices.StratifiedIfPossible)
        RMSE = evalUtilities.RMSE(Res)[0]
        self.assertEqual(round(RMSE,2),round(2.96,2)) #Ver 0.3


        newSVM = svmLearner(self.inDataC)
        trainedAcc = evalUtilities.getRMSE(self.inDataC, newSVM)
        # Save model 
        rc = newSVM.write(self.modelPath)
        self.assertEqual(rc,True)
        # Load the saved model
        loadedsvm = AZorngCvSVM.CvSVMread(self.modelPath)
        loadedAcc = evalUtilities.getRMSE(self.inDataC, loadedsvm)
        # Assure equal accuracy
        self.assertEqual(round(trainedAcc,4), round(2.8289,4)) #Ver 0.3
        self.assertEqual(round(trainedAcc,4), round(loadedAcc,4))




    def test_TwoWays(self):
        """Test two ways svm creation
        Test that an svm created in one or two steps give the same results
        """

        # One step svm creation
        svm = AZorngCvSVM.CvSVMLearner(self.train_data)

        # Calculate classification accuracy for the classifier trained in one step
        oneStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, svm)

        # Two step svm creation
        learner = AZorngCvSVM.CvSVMLearner()
        svm = learner(self.train_data)
        
        # Calculate classification accuracy for the classifier trained in two steps
        twoStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, svm) 

        # Test that the accuracy of the classifiers created in different ways is the exact same
        self.assertEqual(oneStepAcc, twoStepAcc)

    def test_PredictionWithDiffVarType(self):
        """Test prediction with diff. VarType
        Test the prediction of examples with different varType
        """
        expectedAcc = 0.666666666667 
        # Create a svm model
        svm = AZorngCvSVM.CvSVMLearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc2 = evalUtilities.getClassificationAccuracy(self.noBadDataTest[3:],svm)
        Acc1 = evalUtilities.getClassificationAccuracy(self.badVarTypeData[3:],svm)
        self.assertEqual(round(Acc1,7),round(expectedAcc,7),"The Accuracy is not the expected. Got: "+str(Acc1))
        self.assertEqual(round(Acc2,7),round(expectedAcc,7),"The Accuracy is not the expected. Got: "+str(Acc2))    
        self.assert_(('Fixed Types of variables' in svm.examplesFixedLog) and (svm.examplesFixedLog['Fixed Types of variables']==27), "No report of fixing in classifier class")
        self.assert_(('Vars needing type fix' in svm.examplesFixedLog) and (svm.examplesFixedLog['Vars needing type fix']['[Br]([C])']=="EnumVariable to FloatVariable"), "No report of fixing in classifier class")


    def test_PredictionWithDiffVarOrder(self):
        """Test Prediction with diff. VarOrder
        Test the prediction  examples with different varOrder
        """
        expectedAcc = 0.7 # 0.59999999999999998 #0.7 # Ver 0.3
        # Create a svm model
        svm = AZorngCvSVM.CvSVMLearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,svm)
        Acc2 = evalUtilities.getClassificationAccuracy(self.badVarOrderData,svm)

        self.assertEqual(round(Acc1,9),round(expectedAcc,9),"The Accuracy is not the expected. Got: "+str(Acc1)) #Ver 0.3
        self.assertEqual(round(Acc2,9),round(expectedAcc,9),"The Accuracy is not the expected. Got: "+str(Acc2))
        #we do not report order fix anymore! 
        #self.assertEqual(str(svm.examplesFixedLog),"{'Fixed Order of variables': 27}", "No report of fixing in classifier class")

    def test_PredictionWithIncompatibleDomain(self):
        """Test prediction with uncompatible domain
        Test the non-prediction of examples with an incompatible domain  
        """
        expectedAcc1 = 0.7   #Ver 0.3
        # Create a svm model
        svm = AZorngCvSVM.CvSVMLearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,svm)

        self.assertEqual(round(Acc1,9),round(expectedAcc1,9),"The Accuracy is not the expected. Got: "+str(Acc1))
        self.assertEqual(svm(self.badVarTypeData[0]),'NEG',"This example could still be predicted. Got: "+str(svm(self.badVarTypeData[0]))) #Ver 0.3
        self.assertEqual(svm(self.badVarTypeData[1]),'NEG',"This example could still be predicted. Got: "+str(svm(self.badVarTypeData[1])))
        self.assertEqual(svm(self.badVarNameData[0]),None,"This example should NOT be predicted. Got: "+str(svm(self.badVarNameData[0])))
        self.assertEqual(svm(self.badVarCountData[0]),None,"This example should NOT be predicted. Got: "+str(svm(self.badVarCountData[0])))

    def test_ImputeTrain(self):
        """
        Assure that imputation works for the svm models. Test on data with missing values
        This test just assures the the model is trained. The correct imputation test is made on testImpute
        """

        svm = AZorngCvSVM.CvSVMLearner(self.missingTrain)

        Acc = evalUtilities.getClassificationAccuracy(self.missingTest, svm)
        self.assertEqual(round(0.59999999999999998,5),round(Acc,5))# Ver 0.3


    def test_Impute(self):
        """Test missing values imputation
        Assure that imputation works for the svm models. Test on data with missing values
        """
        ex1=self.contTest[1]
        ex2=self.contTest[6]
        self.assert_(ex1["DiscAttr2"]!="?","The var DiscAttr2 shouldn't be missing!")
        self.assert_(ex2["Level"]!="?","The var Level shouldn't be missing!")

        imputer = orange.ImputerConstructor_average(self.contTrain)
        svmL = AZorngCvSVM.CvSVMLearner(p=0.2)  #Here if p=2 and scaleClass is False, is ok, but with p=2 and also scale calss, the model will have no support vectors. So, with p=0.2 and also scale class, it goes right.
        svmL.svm_type = 103
        svm = svmL(self.contTrain)

        # Prediction for data as it is
        P1=svm(ex1)
        P2=svm(ex2)
       
        # Predictions changing one continuous and one discrete variable to 0
        ex1["DiscAttr2"]=0
        ex2["Level"]=0
        P1_0=svm(ex1)
        P2_0=svm(ex2)

        # Predictions changing the same continuous and discrete variable to it's correspondent imputation value
        ex1["DiscAttr2"]=imputer.defaults["DiscAttr2"]
        ex2["Level"]=imputer.defaults["Level"]
        P1_imp=svm(ex1)
        P2_imp=svm(ex2)
 
        # Predictions changing the same continuous and discrete variable to '?' wich means that the same imputation
        # as in the last case will have to be made inside the classifier. So, the predicted value must be the same
        ex1["DiscAttr2"]="?"
        ex2["Level"]="?"
        self.assert_(ex1["DiscAttr2"]=="?","The var DiscAttr2 should be missing now!")
        self.assert_(ex2["Level"]=="?","The var Level should be missing now!")
    
        P1Miss=svm(ex1)
        P2Miss=svm(ex2)

        # Test if the prediction made for the example with mising value is the same as the one 
        # for the example which missing values were substituted using the same method as the classifier does.
        self.assert_(round(P1_imp,4)==round(P1Miss,4),"Imputation was not made correctly inside the classifier")
        self.assert_(round(P2_imp,4)==round(P2Miss,4),"Imputation was not made correctly inside the classifier")

        # Assure that if other substitutions on those variables were made, the predicted value would be different, 
        # and so, this is a valid method for testing the imputation
        self.assert_(round(P1.value,4)!=round(P2.value,4))      # Just to assure that we are not comaring equal examples
        self.assert_(round(P1.value,4)!=round(P1_imp.value,4))
        self.assert_(round(P1_0.value,4)!=round(P1_imp.value,4))
        self.assert_(round(P2.value,4)!=round(P2_imp.value,4))
        self.assert_(round(P2_0.value,4)!=round(P2_imp.value,4))


        #Test the imputer for saved models
        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirSVMtest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"CvSVMModel")
        svm.write(modelPath)

        # Read in the model
        svmM = AZorngCvSVM.CvSVMread(modelPath)
        # Predict the ex1 and ex2 which are still the examples with missing values '?'
        self.assert_( ex1["DiscAttr2"]=="?","Value of Var DiscAttr2 should be missing!")
        self.assert_( ex2["Level"]=="?","Value of Var Level should be missing!")
        self.assert_(round(svmM(ex1),4)==round(P1Miss,4),"Imputation on loaded model is not correct")
        self.assert_(round(svmM(ex2),4)==round(P2Miss,4),"Imputation on loaded model is not correct")
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


    def test_MetaDataHandle(self):
        """Test the handling of Data with Meta Atributes
        """
        # Create an svm model

        svm = AZorngCvSVM.CvSVMLearner(self.NoMetaTrain)

        # Calculate classification accuracy (NoMetaTest and WMeta are the same appart from the meta atribute) 
        AccNoMeta = evalUtilities.getClassificationAccuracy(self.NoMetaTest, svm)
        AccWMeta = evalUtilities.getClassificationAccuracy(self.WMetaTest, svm)

        self.assertEqual(AccNoMeta,AccWMeta,"Predictions with and without meta data were different!")
        self.assertEqual(round(AccNoMeta,9), round(0.7,9),"Accuracy was not the expected value! Got: ") #Ver 0.3
        

    def test_MetaDataHandleForSavingModel(self):
        """Test the handling of SaveModel for Data with Meta Atributes
        """

        #Test the save of a model created from a train data with meta attributes
        self.assert_(len(self.WMetaTest.domain.getmetas())>=1,"The dataset WMetaTest should have Meta Attributes")
        svmM = AZorngCvSVM.CvSVMLearner(self.WMetaTest)
        AccNoMetaBefore = evalUtilities.getClassificationAccuracy(self.NoMetaTrain,svmM) 
        AccWMetaBefore = evalUtilities.getClassificationAccuracy(self.WMetaTest,svmM)

        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirSVMtest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"CvSVMModel")
        svmM.write(modelPath)

        # Read in the model
        svmR = AZorngCvSVM.CvSVMread(modelPath)
        self.assert_(len(svmR.imputer.defaults.domain.getmetas())==0,"There shouldn't be any Meta data now!")

        # Calculate classification accuracy 
        AccNoMetaAfter = evalUtilities.getClassificationAccuracy(self.NoMetaTrain, svmR)
        AccWMetaAfter = evalUtilities.getClassificationAccuracy(self.WMetaTest, svmR)

        # Test that the accuracy of the model before and after saved
        self.assertEqual(AccNoMetaBefore, AccNoMetaAfter,"NoMeta: Predictions after loading saved model were different")
        self.assertEqual(AccWMetaBefore, AccWMetaAfter, "WMeta: Predictions after loading saved model were different")
        self.assertEqual(round(AccWMetaAfter,9), round(0.7,9),"Accuracy was not the expected value!")
        self.assertEqual(round(AccNoMetaAfter,9), round(0.6,9),"Accuracy was not the expected value!")
 
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)




if __name__ == "__main__":
    scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirSVMtest"+str(time.time()))
    os.mkdir(scratchdir)

    suite = unittest.TestLoader().loadTestsFromTestCase(SVMClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    if os.path.isdir(scratchdir):
        os.system("/bin/rm -rf "+scratchdir)



