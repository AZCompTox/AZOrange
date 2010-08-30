from AZutilities import dataUtilities
import unittest
import os
import time

import orange
from trainingMethods import AZorngPLS
from AZutilities import evalUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngImpute


class PLSClassifierTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        """Creates the training and testing data set attributes. """
        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris.tab")
        contDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")
        contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")
        dataNoMetaTrainPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyTrain.tab")
        missingTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy_missing.tab")

        #These 2 datasets are equal apart from the meta atribute
        dataNoMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyTest.tab")
        dataWMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummywMeta.tab")

        # Read in the data
        missingInData = dataUtilities.DataTable(missingTestDataPath)
        contTrainData = dataUtilities.DataTable(contTrainDataPath)
        contData = dataUtilities.DataTable(contDataPath)
        inData = dataUtilities.DataTable(testDataPath)
        self.NoMetaTrain = dataUtilities.DataTable(dataNoMetaTrainPath)
        self.NoMetaTest = dataUtilities.DataTable(dataNoMetaTestPath)
        self.WMetaTest = dataUtilities.DataTable(dataWMetaTestPath)

        # Random sampling
        self.missingTrain = missingInData
        self.missingTest = missingInData
        self.train_data, self.test_data = self.randSamp(inData, 0.7)
        self.contTrain = contTrainData
        self.contTest = contData
        
        #Data for domain fix handling
        
        badVarTypePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyBadVarType.tab")
        badVarNamePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyBadVarName.tab")
        badVarOrderPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyBadVarOrder.tab")
        badVarCountPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyBadVarCount.tab")
        # Read in the data
        self.noBadDataTrain = self.NoMetaTrain
        self.noBadDataTest = self.NoMetaTest
        self.badVarTypeData = dataUtilities.DataTable(badVarTypePath)
        self.badVarNameData = dataUtilities.DataTable(badVarNamePath)
        self.badVarOrderData = dataUtilities.DataTable(badVarOrderPath)
        self.badVarCountData = dataUtilities.DataTable(badVarCountPath)   #One less example 


    def testTwoWays(self):
        """Test two ways pls creation
        Test that an pls created in one or two steps give the same results
        """

        # One step pls creation
        pls = AZorngPLS.PLSLearner(self.train_data)

        # Calculate classification accuracy for the classifier trained in one step
        oneStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, pls)

        # Two step pls creation
        learner = AZorngPLS.PLSLearner()
        pls = learner(self.train_data)
        
        # Calculate classification accuracy for the classifier trained in two steps
        twoStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, pls) 

        # Test that the accuracy of the classifiers created in different ways is the exact same
        self.assertEqual(oneStepAcc, twoStepAcc)


    def testPersistentClassAcc(self):
        """Test the persistence of the learner as Classifier
        Assure that the accuracy is perserved for models trained in the same way. 
        """
        PLSlearner = AZorngPLS.PLSLearner(k="5", method="kernel", precision="1e-6")
        PLSClassifier = PLSlearner(self.NoMetaTrain)

        # Calculate classification accuracy 
        ClassifierAcc = evalUtilities.getClassificationAccuracy(self.NoMetaTest, PLSClassifier)

        # Check that the accuracy is what it used to be
        self.assertEqual(round(0.851851851852,9),round(ClassifierAcc,9))

    def testPersistentRegAcc(self):
        """Test the persistence of the learner as Regressor
        Assure that the RMSE is perserved for models trained in the same way. 
        """
        PLSlearner = AZorngPLS.PLSLearner(k="3", method="simpls", precision="1e-6")
        PLSRegressor = PLSlearner(self.contTrain)

        # Calculate RMSE 
        RegressorRMSE = evalUtilities.getRMSE(self.contTest, PLSRegressor)

        # Check that RMSE is what it used to be
        self.assertEqual(round(0.67276457999999995,8),round(RegressorRMSE,8))


    def testSavedModel(self):
        """Test PLS model saving
        Test to assure that a saved pls model gives the same predictions as before saving."""

        # Create a pls model
        pls = AZorngPLS.PLSLearner(self.train_data)

        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.test_data, pls)

        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdir"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"PLSModel")
        pls.write(modelPath)
        
        # Read in the model
        plsM = AZorngPLS.PLSread(modelPath)

        # Calculate classification accuracy 
        savedAcc = evalUtilities.getClassificationAccuracy(self.test_data, plsM)

        # Test that the accuracy of the two classifiers is the exact same
        self.assertEqual(Acc, savedAcc)

        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


#Datasets for testing Bad Data format
# self.noBadDataTrain
# self.noBadDataTest
# self.badVarTypeData
# self.badVarNameData 
# self.badVarOrderData
# self.badVarCountData

    def testPredictionWithDiffVarType(self):
        """Test prediction with diff. VarType
        Test the prediction of examples with different varType
        """
        expectedAcc = 0.875
        # Create a pls model
        pls = AZorngPLS.PLSLearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc2 = evalUtilities.getClassificationAccuracy(self.noBadDataTest[3:],pls)
        Acc1 = evalUtilities.getClassificationAccuracy(self.badVarTypeData[3:],pls)
        self.assertEqual(Acc1,expectedAcc,"The Accuracy is not the expected")
        self.assertEqual(Acc2,expectedAcc,"The Accuracy is not the expected")    
        self.assert_(('Fixed Types of variables' in pls.examplesFixedLog) and (pls.examplesFixedLog['Fixed Types of variables']==24), "No report of fixing in classifier class")
        self.assert_(('Vars needing type fix' in pls.examplesFixedLog) and (pls.examplesFixedLog['Vars needing type fix']['SELMA_Max_pos_chrg_GH']=="EnumVariable to FloatVariable"), "No report of fixing in classifier class")


    def testPredictionWithDiffVarOrder(self):
        """Test Prediction with diff. VarOrder
        Test the prediction  examples with different varOrder
        """
        expectedAcc = 0.851851851852
        # Create a pls model
        pls = AZorngPLS.PLSLearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,pls)
        Acc2 = evalUtilities.getClassificationAccuracy(self.badVarOrderData,pls)
        self.assertEqual(round(Acc1,9),round(expectedAcc,9),"The Accuracy is not the expected")
        self.assertEqual(round(Acc2,9),round(expectedAcc,9),"The Accuracy is not the expected")
        #we do not report order fix anymore!         
        #self.assertEqual(str(pls.examplesFixedLog),"{'Fixed Order of variables': 27}", "No report of fixing in classifier class")

    def testPredictionWithIncompatibleDomain(self):
        """Test prediction with uncompatible domain
        Test the non-prediction of examples with an incompatible domain  
        """
        expectedAcc1 = 0.851851851852
        expectedAcc2 = 0.875
        # Create a pls model
        pls = AZorngPLS.PLSLearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,pls)
        self.assertEqual(round(Acc1,9),round(expectedAcc1,9),"The Accuracy is not the expected")
        self.assertEqual(pls(self.badVarTypeData[0]),"POS","This example could still be predicted")
        self.assertEqual(pls(self.badVarTypeData[1]),"NEG","This example could still be predicted")
        print pls(self.badVarNameData[0])
        self.assertEqual(pls(self.badVarNameData[0]),None,"This example should NOT  be predicted")
        self.assertEqual(pls(self.badVarCountData[0]),None,"This example should NOT be predicted")

    def testImputeTrain(self):
        """
        Assure that imputation works for the pls models. Test on data with missing values
        This test just assures the the model is trained. The correct imputation test is made on testImpute
        """

        pls = AZorngPLS.PLSLearner(self.missingTrain)

        Acc = evalUtilities.getClassificationAccuracy(self.missingTest, pls)

        self.assertEqual(round(0.84614999999999996,3),round(Acc,3))


    def testImpute(self):
        """Test missing values imputation
        Assure that imputation works for the pls models. Test on data with missing values
        """
        ex1=self.contTest[0]
        ex2=self.contTest[3]
        self.assert_(ex1["ClogP"]!="?","The var ClogP shouldn't be missing!")
        self.assert_(ex2["Lipinski"]!="?","The var Lipinski shouldn't be missing!")

        imputer = orange.ImputerConstructor_average(self.contTrain)

        pls = AZorngPLS.PLSLearner(self.contTrain)

        # Prediction for data as it is
        P1=pls(ex1)
        P2=pls(ex2)
       
        # Predictions changing one continuous and one discrete variable to 0
        ex1["ClogP"]=0
        ex2["Lipinski"]=0
        P1_0=pls(ex1)
        P2_0=pls(ex2)

        # Predictions changing the same continuous and discrete variable to it's correspondent imputation value
        ex1["ClogP"]=imputer.defaults["ClogP"]
        ex2["Lipinski"]=imputer.defaults["Lipinski"]
        P1_imp=pls(ex1)
        P2_imp=pls(ex2)
 
        # Predictions changing the same continuous and discrete variable to '?' wich means that the same imputation
        # as in the last case will have to be made inside the classifier. So, the predicted value must be the same
        ex1["ClogP"]="?"
        ex2["Lipinski"]="?"
        self.assert_(ex1["ClogP"]=="?","The var ClogP should be missing now!")
        self.assert_(ex2["Lipinski"]=="?","The var Lipinski should be missing now!")
    
        P1Miss=pls(ex1)
        P2Miss=pls(ex2)


        # Test if the prediction made for the example with mising value is the same as the one 
        # for the example which missing values were substituted using the same method as the classifier does.
        self.assert_(P1_imp==P1Miss,"Imputation was not made correctly inside the classifier")
        self.assert_(P2_imp==P2Miss,"Imputation was not made correctly inside the classifier")

        # Assure that if other substitutions on those variables were made, the predicted value would be different, 
        # and so, this is a valid method for testing the imputation
        self.assert_(P1.value!=P2.value)      # Just to assure that we are not comaring equal examples
        self.assert_(P1.value!=P1_imp.value)
        self.assert_(P1_0.value!=P1_imp.value)
        self.assert_(P2.value!=P2_imp.value)
        self.assert_(P2_0.value!=P2_imp.value)


        #Test the imputer for saved models
        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdir"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"PLSModel")
        pls.write(modelPath)

        # Read in the model
        plsM = AZorngPLS.PLSread(modelPath)
        # Predict the ex1 and ex2 which are still the examples with missing values '?'
        self.assert_( ex1["ClogP"]=="?","Value of Var ClogP should be missing!")
        self.assert_( ex2["Lipinski"]=="?","Value of Var Lipinski should be missing!")
        self.assert_(plsM(ex1)==P1Miss,"Imputation on loaded model is not correct")
        self.assert_(plsM(ex2)==P2Miss,"Imputation on loaded model is not correct")
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


    def testMetaDataHandle(self):
        """Test the handling of Data with Meta Atributes
        """
        # Create an pls model

        pls = AZorngPLS.PLSLearner(self.NoMetaTrain)

        # Calculate classification accuracy (NoMetaTest and WMeta are the same appart from the meta atribute) 
        AccNoMeta = evalUtilities.getClassificationAccuracy(self.NoMetaTest, pls)
        AccWMeta = evalUtilities.getClassificationAccuracy(self.WMetaTest, pls)

        self.assertEqual(AccNoMeta,AccWMeta,"Predictions with and without meta data were different!")
        self.assertEqual(round(AccNoMeta,9), round(0.851851851852,9),"Accuracy was not the expected value!")
        

    def testMetaDataHandleForSavingModel(self):
        """Test the handling of SaveModel for Data with Meta Atributes
        """

        #Test the save of a model created from a train data with meta attributes
        self.assert_(len(self.WMetaTest.domain.getmetas())>=1,"The dataset WMetaTest should have Meta Attributes")
        plsM = AZorngPLS.PLSLearner(self.WMetaTest)
        AccNoMetaBefore = evalUtilities.getClassificationAccuracy(self.NoMetaTrain,plsM) 
        AccWMetaBefore = evalUtilities.getClassificationAccuracy(self.WMetaTest,plsM)


        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdir"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"PLSModel")
        plsM.write(modelPath)

        # Read in the model
        plsR = AZorngPLS.PLSread(modelPath)
        self.assert_(len(plsR.imputer.defaults.domain.getmetas())==0,"There shouldn't be any Meta data now!")

        # Calculate classification accuracy 
        AccNoMetaAfter = evalUtilities.getClassificationAccuracy(self.NoMetaTrain, plsR)
        AccWMetaAfter = evalUtilities.getClassificationAccuracy(self.WMetaTest, plsR)

        # Test that the accuracy of the model before and after saved
        self.assertEqual(AccNoMetaBefore, AccNoMetaAfter,"NoMeta: Predictions after loading saved model were different")
        self.assertEqual(AccWMetaBefore, AccWMetaAfter, "WMeta: Predictions after loading saved model were different")
        self.assertEqual(round(AccWMetaAfter,9), round(0.888888888889,9),"Accuracy was not the expected value!")
        self.assertEqual(round(AccNoMetaAfter,9), round(0.605769230769,9),"Accuracy was not the expected value!")
 
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)




if __name__ == "__main__":
        suite = unittest.TestLoader().loadTestsFromTestCase(PLSClassifierTest)
        unittest.TextTestRunner(verbosity=2).run(suite)


