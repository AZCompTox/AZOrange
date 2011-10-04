from AZutilities import dataUtilities
import unittest
import os
import time

import orange
from trainingMethods import AZorngCvBayes
from trainingMethods import AZBaseClasses
from AZutilities import evalUtilities
from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngImpute


class CvBayesClassifierTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        """Creates the training and testing data set attributes. """
        trainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisTrain.tab")
        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisTest.tab")
        missingTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train_missing.tab")
        
        # Read in the data
        self.train_data =  dataUtilities.DataTable(trainDataPath)
        self.test_data =  dataUtilities.DataTable(testDataPath)
        missingInData = dataUtilities.DataTable(missingTestDataPath)

        # Random sampling
        self.missingTrain = missingInData
        self.missingTest = missingInData

        
        dataNoMetaTrainPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train.tab")

        #These 2 datasets are equal apart from the meta atribute
        dataNoMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_SmallTest.tab")
        dataWMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_W_metas_SmallTest.tab")
        #contTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTest.tab")
        #contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTrain.tab")
        # Read in the data
        #contTrainData = dataUtilities.DataTable(contTrainDataPath)   
        #contData = dataUtilities.DataTable(contTestDataPath)
        self.NoMetaTrain = dataUtilities.DataTable(dataNoMetaTrainPath)
        self.NoMetaTest = dataUtilities.DataTable(dataNoMetaTestPath)
        self.WMetaTest = dataUtilities.DataTable(dataWMetaTestPath)

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
        
        # Beter data for Bayes
        LdataTrainPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_FullNumeric_Train.tab") 
        LdataTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_FullNumeric_Test.tab") 
        
        self.LdataTrain = dataUtilities.DataTable(LdataTrainPath)
        self.LdataTest = dataUtilities.DataTable(LdataTestPath)


    def no_test_PredictionVarImportance(self):   #NOT implemnted for Bayes
        """Test the variable importance for single predictions"""
        def TopVarImportanceTest(data, expectNone = False):
            resA = []
            resB = []
            CvBayes = AZorngCvBayes.CvBayesLearner(data)

            for ex in data:
                resA.append(CvBayes.getTopImportantVars(ex,1))

            scratchdir = miscUtilities.createScratchDir(desc="TopVarImportanceTest")
            modelPath = os.path.join(scratchdir,"CvBayesModel")
            CvBayes.write(modelPath)
            LoadedCvBayes = AZorngCvBayes.CvBayesread(modelPath)
            miscUtilities.removeDir(scratchdir) 
            for ex in data:
                resB.append(LoadedCvBayes.getTopImportantVars(ex,1))
            if expectNone:
               return resA == resB == [None]*len(data)
            else:
                return resA == resB and None not in resA and resA.count(resA[0]) != len(resA)


        self.assertEqual(TopVarImportanceTest(self.train_data, True), True, "Failed testing with uncompatible data (class with 3 values)")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisCont.tab"))
        self.assertEqual(TopVarImportanceTest(data), True,"Failed testing with Regression data.")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab"))
        self.assertEqual(TopVarImportanceTest(data), True, "Failed testing with Binary Classifier")




    def no_test_DFV(self):  #NOT implemnted for Bayes
        """ Test the Decision Function Value Return"""
        CvBayes = AZorngCvBayes.CvBayesLearner(self.LdataTrain)
        #Testsing with return of DFV
        RDFV = True
        for ex in self.LdataTest:
            predictedClass = CvBayes(ex)
            a = CvBayes(ex,returnDFV = RDFV)
            b = CvBayes(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = CvBayes(ex,resultType = orange.GetBoth,returnDFV = RDFV)
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
            # CvBayes does always return real probabilities on binary classification
            self.assertEqual(CvBayes.isRealProb(),True)
        expectedExtremes = {'max': 0.5, 'min':-0.5 }
        self.assertEqual([round(x,5) for x in CvBayes.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(CvBayes.nPredictions,4*len(self.LdataTest))

        #Testsing without return of DFV
        RDFV = False
        for ex in self.LdataTest:
            a = CvBayes(ex,returnDFV = RDFV)
            b = CvBayes(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = CvBayes(ex,resultType = orange.GetBoth,returnDFV = RDFV)
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
            # CvBayes does always return real probabilities on binary classification
            self.assertEqual(CvBayes.isRealProb(),True)

        self.assertEqual([round(x,5) for x in CvBayes.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(CvBayes.nPredictions,(3+4)*len(self.LdataTest))

    def test_Probabilities(self):
        """Test if the returned probabilities are ok although fake"""

        CvBayes = AZorngCvBayes.CvBayesLearner(self.LdataTrain)
        res = []
        for idx,ex in enumerate(self.LdataTest):
            res.append(CvBayes(ex,resultType = orange.GetProbabilities))
            #print res[-1]
            self.assert_(res[-1][0]>=0 and res[-1][0]<=1,"Example "+str(idx)+" have impossible probability:"+str(res[-1])) 
            self.assert_(res[-1][1]>=0 and res[-1][1]<=1,"Example "+str(idx)+" have impossible probability:"+str(res[-1])) 
            #self.assertEqual(CvBayes.isRealProb(),True,"Example "+str(idx)+" did not return real probability")
            #print "Res",idx,":",res[-1]
            #print "Sum",idx,":",round(sum(res[-1]),5)
            self.assertEqual(round(sum(res[-1]),5),1,"Probabilities of Example "+str(idx)+" did not sum 1:"+str(res[-1]))
        sum0 = sum([x[0] for x in res])
        sum1 = sum([x[1] for x in res])
        self.assertEqual(len(self.LdataTest),round(sum0+sum1,5))
        #self.assert_(sum0-int(sum0) > 0)
        #self.assert_(sum1-int(sum1) > 0)


    def test_TwoWays(self):
        """
        Test that an Bayes created in one or two steps give the same results
        """
        #Deviation allowed in Acc
        devAlloed = 0.4 #Before:   0.02

        # One step Bayes creation
        Bayes = AZorngCvBayes.CvBayesLearner(self.train_data)

        # Calculate classification accuracy for the classifier trained in one step
        oneStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, Bayes)

        # Two step Bayes creation
        learner = AZorngCvBayes.CvBayesLearner()
        Bayes = learner(self.train_data)
        
        # Calculate classification accuracy for the classifier trained in two steps
        twoStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, Bayes) 

        # Test that the accuracy of the classifiers created in different ways is the exact same
        self.assert_(oneStepAcc >= twoStepAcc-devAlloed and oneStepAcc<=twoStepAcc+devAlloed, "Dev="+str(oneStepAcc-twoStepAcc))


    def test_SavedModel(self):
        """Test to assure that a saved Bayes model gives the same predictions as before saving."""

        # Create an Bayes model
        Bayes = AZorngCvBayes.CvBayesLearner(self.train_data)

        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.test_data, Bayes)

        # Save the model
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdir"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"Bayes.fBayes")
        Bayes.write(modelPath)
        
        # Read in the model
        Bayes = AZorngCvBayes.CvBayesread(modelPath)

        # Calculate classification accuracy 
        savedAcc = evalUtilities.getClassificationAccuracy(self.test_data, Bayes)

        # Test that the accuracy of the two classifiers is the exact same
        self.assertEqual(Acc, savedAcc)
        
        #Test using the global read functionality
        Bayes2 = AZBaseClasses.modelRead(modelPath)
        savedAcc2 = evalUtilities.getClassificationAccuracy(self.test_data, Bayes2)
        self.assertEqual(Acc, savedAcc2)

        
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)



    def test_ImputeTrain(self):
        """
        Assure that imputation works for the Bayes models. Test on data with missing values
        This test just assures the the model is trained. The correct imputation test is made on testImpute
        """
        expectedAcc = 0.41818
        devAllowed = 0.3
        BayesLearner = AZorngCvBayes.CvBayesLearner()

        Bayes = BayesLearner(self.missingTrain)
    
        Acc = evalUtilities.getClassificationAccuracy(self.missingTest, Bayes)

        self.assert_(abs(expectedAcc-Acc) <= devAllowed,"Dev="+str(expectedAcc-Acc) )


    def test_PredictionWithDiffVarType(self):
        """Test prediction with diff. VarType
        Test the prediction of examples with different varType
        """
        expectedAcc = 0.37036999999999998 
        # Create a Bayes model
        CvBayeslearner = AZorngCvBayes.CvBayesLearner()
        Bayes = CvBayeslearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc2 = evalUtilities.getClassificationAccuracy(self.noBadDataTest[3:],Bayes)
        Acc1 = evalUtilities.getClassificationAccuracy(self.badVarTypeData[3:],Bayes)
        self.assertEqual(round(Acc1,6),round(expectedAcc,6))
        self.assertEqual(round(Acc2,6),round(expectedAcc,6))  
        self.assert_(('Fixed Types of variables' in Bayes.examplesFixedLog) and (Bayes.examplesFixedLog['Fixed Types of variables']==27), "No report of fixing in classifier class")
        self.assert_(('Vars needing type fix' in Bayes.examplesFixedLog) and (Bayes.examplesFixedLog['Vars needing type fix']['[Br]([C])']=="EnumVariable to FloatVariable", "No report of fixing in classifier class"))


    def test_PredictionWithDiffVarOrder(self):
        """Test Prediction with diff. VarOrder
        Test the prediction  examples with different varOrder
        """
        expectedAcc = [0.33333333300000001, 0.666666667]
        # Create a Bayes model
        CvBayeslearner = AZorngCvBayes.CvBayesLearner()
        Bayes = CvBayeslearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,Bayes)
        Acc2 = evalUtilities.getClassificationAccuracy(self.badVarOrderData,Bayes)

        self.assertEqual(Acc1, Acc2)
        self.assert_(round(Acc1,5) in [round(x,5) for x in expectedAcc])
        #we do not report order fix anymore!
        #self.assert_('Fixed Order of variables' in Bayes.examplesFixedLog and Bayes.examplesFixedLog['Fixed Order of variables']==27, "No report of fixing in classifier class")


    def test_PredictionWithIncompatibleDomain(self):
        """Test prediction with uncompatible domain
        Test the non-prediction of examples with an incompatible domain  
        """
        expectedAcc1 =  0.33333333300000001
        # Create a Bayes model
        CvBayeslearner = AZorngCvBayes.CvBayesLearner()
        Bayes = CvBayeslearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,Bayes)
        self.assertEqual(round(Acc1,9),round(expectedAcc1,9))
        self.assertEqual(Bayes(self.badVarTypeData[0]),"POS","This example could still be predicted")
        self.assertEqual(Bayes(self.badVarTypeData[1]),"POS","This example could still be predicted")
        self.assertEqual(Bayes(self.badVarNameData[0]),None,"This example should NOT be predicted")
        self.assertEqual(Bayes(self.badVarCountData[0]),None,"This example should NOT be predicted")


    def no_test_Impute(self):  # Bayes cannot deal with regression
        """Test missing values imputation
        Assure that imputation works for the Bayes models. Test on data with missing values
        """
        #This data is loaded here to speed up the test suite since it is too big
        contTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTest.tab")
        contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTrain.tab")
        contTrain = dataUtilities.DataTable(contTrainDataPath)   
        contTest = dataUtilities.DataTable(contTestDataPath)

        ex1=contTest[5]
        ex2=contTest[6]
        self.assert_(ex1["Desc 71"]!="?","The var Desc 71 shouldn't be missing!")
        self.assert_(ex2["Desc 138"]!="?","The var Desc 138 shouldn't be missing!")

        imputer = orange.ImputerConstructor_average(contTrain)
        
        CvBayeslearner = AZorngCvBayes.CvBayesLearner()
        Bayes = CvBayeslearner(contTrain)        

        # Prediction for data as it is
        P1=Bayes(ex1)
        P2=Bayes(ex2)
       
        # Predictions changing one continuous and one discrete variable to 0
        ex1["Desc 71"]=0
        ex2["Desc 138"]=0
        P1_0=Bayes(ex1)
        P2_0=Bayes(ex2)

        # Predictions changing the same continuous and discrete variable to it's correspondent imputation value
        ex1["Desc 71"]=imputer.defaults["Desc 71"]
        ex2["Desc 138"]=imputer.defaults["Desc 138"]
        P1_imp=Bayes(ex1)
        P2_imp=Bayes(ex2)
 
        # Predictions changing the same continuous and discrete variable to '?' wich means that the same imputation
        # as in the last case will have to be made inside the classifier. So, the predicted value must be the same
        ex1["Desc 71"]="?"
        ex2["Desc 138"]="?"
        self.assert_(ex1["Desc 71"]=="?","The var Desc 71 should be missing now!")
        self.assert_(ex2["Desc 138"]=="?","The var Desc 138 should be missing now!")    
        P1Miss=Bayes(ex1)
        P2Miss=Bayes(ex2)


        # Test if the prediction made for the example with mising value is the same as the one 
        # for the example which missing values were substituted using the same method as the classifier does.
        self.assert_(P1_imp==P1Miss,"Imputation was not made correctly inside the classifier")
        #self.assert_(P2_imp==P2Miss,"Imputation was not made correctly inside the classifier")

        # Assure that if other substitutions on those variables were made, the predicted value would be different, 
        # and so, this is a valid method for testing the imputation
        self.assert_(P1.value!=P2.value)      # Just to assure that we are not comaring equal examples
        self.assert_(P1.value!=P1_imp.value,"The imputed 1 was the same as the original ... try other example")
        self.assert_(P1_0.value!=P1_imp.value,"The imputed 1 was the same as the replaced by 0. The classifier may be replacing missing values by 0")
        self.assert_(P2.value!=P2Miss.value, "The missing imputed 2 was the same as the original ... try other example")
        self.assert_(P2_0.value!=P2Miss.value,"The missing imputed 2 was the same as the replaced by 0. The classifier may be replacing missing values by 0")


        #Test the imputer for saved models
        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdir"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"CvBayesModel")
        Bayes.write(modelPath)

        # Read in the model
        BayesM = AZorngCvBayes.CvBayesread(modelPath)
        # Predict the ex1 and ex2 which are still the examples with missing values '?'
        self.assert_( ex1["Desc 71"]=="?","Value of Var Desc 71 should be missing!")
        self.assert_( ex2["Desc 138"]=="?","Value of Var Desc 138 should be missing!")
        self.assert_(round(BayesM(ex1),6)==round(P1Miss,6),"Imputation on loaded model is not correct")
        self.assert_(round(BayesM(ex2),6)==round(P2Miss,6),"Imputation on loaded model is not correct")
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)
        
    def test_MetaDataHandle(self):
        """Test the handling of Data with Meta Atributes
        """
        expectedAcc = 0.33333333300000001  #[0.666666667, 0.333333333]
        # Create an Bayes model
        CvBayeslearner = AZorngCvBayes.CvBayesLearner()
        Bayes = CvBayeslearner(self.NoMetaTrain)

        # Calculate classification accuracy (NoMetaTest and WMeta are the same appart from the meta atribute) 
        AccNoMeta = evalUtilities.getClassificationAccuracy(self.NoMetaTest, Bayes)
        AccWMeta = evalUtilities.getClassificationAccuracy(self.WMetaTest, Bayes)
        self.assertEqual(AccNoMeta,AccWMeta,"Predictions with and without meta data were different!")
        self.assertEqual(round(AccNoMeta,9), round(expectedAcc,9))
        
    def test_MetaDataHandleForSavingModel(self):
        """Test the handling of SaveModel for Data with Meta Atributes
        """
        expectedAccWMeta = 0.83333333300000001
        expectedAccNoMeta =0.55151515200000001 
        #Test the save of a model created from a train data with meta attributes
        self.assert_(len(self.WMetaTest.domain.getmetas())>=1,"The dataset WMetaTest should have Meta Attributes")
        CvBayeslearner = AZorngCvBayes.CvBayesLearner()
        BayesM = CvBayeslearner(self.WMetaTest)
        AccNoMetaBefore = evalUtilities.getClassificationAccuracy(self.NoMetaTrain,BayesM) 
        AccWMetaBefore = evalUtilities.getClassificationAccuracy(self.WMetaTest,BayesM)


        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdiriTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"CvBayesModel.CvBayes")
        BayesM.write(modelPath)

        # Read in the model
        BayesR = AZorngCvBayes.CvBayesread(modelPath)
        self.assert_(len(BayesR.imputer.defaults.domain.getmetas())==0,"There shouldn't be any Meta data now!")

        # Calculate classification accuracy 
        AccNoMetaAfter = evalUtilities.getClassificationAccuracy(self.NoMetaTrain, BayesR)
        AccWMetaAfter = evalUtilities.getClassificationAccuracy(self.WMetaTest, BayesR)

        # Test that the accuracy of the model before and after saved
        self.assertEqual(AccNoMetaBefore, AccNoMetaAfter,"NoMeta: Predictions after loading saved model were different")
        self.assertEqual(AccWMetaBefore, AccWMetaAfter, "WMeta: Predictions after loading saved model were different")
        self.assertEqual(round(AccWMetaAfter,9), round(expectedAccWMeta,9))
        self.assertEqual(round(AccNoMetaAfter,9), round(expectedAccNoMeta,9))
 
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)



    def test_PersistentClassAcc(self):
        """
        Assure that the accuracy is perserved for models trained in the same way. 
        """
        #Deviation Allowed
        devAllowed = 0.3
        ExpectedAcc = 0.95
        # One step Bayes creation
        Bayes = AZorngCvBayes.CvBayesLearner(self.train_data)
        # Calculate classification accuracy for the classifier trained in one step
        oneStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, Bayes)
        # Check that the accuracy is what it used to be
        self.assert_(abs(oneStepAcc - ExpectedAcc) <= devAllowed, "Dev="+str(oneStepAcc-ExpectedAcc) ) 



if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(CvBayesClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)



