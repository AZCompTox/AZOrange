from AZutilities import dataUtilities
import unittest
import os
import time

import orange
from trainingMethods import AZorngCvANN
from AZutilities import evalUtilities
from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngImpute


class CvANNClassifierTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        """Creates the training and testing data set attributes. """
        missingTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train_missing.tab")
        
        # Read in the data
        #IrisData

        irisDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris.tab")
        self.irisData = dataUtilities.DataTable(irisDataPath)

        trainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisTrain.tab")
        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisTest.tab")
        self.train_data =  dataUtilities.DataTable(trainDataPath)
        self.test_data =  dataUtilities.DataTable(testDataPath)

        missingInData = dataUtilities.DataTable(missingTestDataPath)
        self.missingTrain = missingInData
        self.missingTest = missingInData

        ##scPA
        
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
        
        # Beter data for ANN
        LdataTrainPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_FullNumeric_Train.tab") 
        LdataTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_FullNumeric_Test.tab") 
        
        self.LdataTrain = dataUtilities.DataTable(LdataTrainPath)
        self.LdataTest = dataUtilities.DataTable(LdataTestPath)

        ##ecPA

    def test_Priors(self):
        """Test to assure that priors are set correcly."""

        # Create a CvANN model
        CvANNlearner = AZorngCvANN.CvANNLearner(stopUPs = 0, priors = {"Iris-versicolor":0.35, "Iris-virginica":0.13, "Iris-setosa":0.52})
        CvANNmodel = CvANNlearner(self.irisData)
        #Model with No Priors
        CvANNlearnerNoP = AZorngCvANN.CvANNLearner(stopUPs=0)
        CvANNmodelNoP = CvANNlearnerNoP(self.irisData)


        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.irisData, CvANNmodel)

        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"modelPriors.CvANN")
        CvANNmodel.write(modelPath)

        # Read in the model
        newCvANNmodel = AZorngCvANN.CvANNread(modelPath)

        # Calculate classification accuracy 
        savedAcc = evalUtilities.getClassificationAccuracy(self.irisData, CvANNmodel)
        NoPAcc = evalUtilities.getClassificationAccuracy(self.irisData, CvANNmodelNoP)

        # Test that the accuracy of the two classifiers is the exact same
        self.assertEqual(Acc, savedAcc)
        self.assert_(Acc != NoPAcc)


        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


    def test_PredictionVarImportance(self):
        """Test the variable importance for single predictions"""
        def TopVarImportanceTest(data, expectNone = False):
            resA = []
            resB = []
            CvANN = AZorngCvANN.CvANNLearner(data, stopUPs=0)

            for ex in data:
                resA.append(CvANN.getTopImportantVars(ex,1))

            scratchdir = miscUtilities.createScratchDir(desc="TopVarImportanceTest")
            modelPath = os.path.join(scratchdir,"CvANNModel")
            CvANN.write(modelPath)
            LoadedCvANN = AZorngCvANN.CvANNread(modelPath)
            miscUtilities.removeDir(scratchdir) 
            for ex in data:
                resB.append(LoadedCvANN.getTopImportantVars(ex,1))
            if expectNone:
               return resA == resB == [None]*len(data)
            else:
                return resA == resB and None not in resA and resA.count(resA[0]) != len(resA)


        self.assertEqual(TopVarImportanceTest(self.train_data, True), True, "Failed testing with uncompatible data (class with 3 values)")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisCont.tab"))
        self.assertEqual(TopVarImportanceTest(data), True,"Failed testing with Regression data.")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab"))
        self.assertEqual(TopVarImportanceTest(data), True, "Failed testing with Binary Classifier")




    def test_DFV(self):
        """ Test the Decision Function Value Return"""
        CvANN = AZorngCvANN.CvANNLearner(self.LdataTrain,stopUPs=0)
        #Testsing with return of DFV
        RDFV = True
        for ex in self.LdataTest:
            predictedClass = CvANN(ex)
            a = CvANN(ex,returnDFV = RDFV)
            b = CvANN(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = CvANN(ex,resultType = orange.GetBoth,returnDFV = RDFV)
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
            # CvANN does always return real probabilities on binary classification
            self.assertEqual(CvANN.isRealProb(),True)
        expectedExtremes = {'max': 0.5, 'min':-0.5 }
        self.assertEqual([round(x,5) for x in CvANN.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(CvANN.nPredictions,4*len(self.LdataTest))

        #Testsing without return of DFV
        RDFV = False
        for ex in self.LdataTest:
            a = CvANN(ex,returnDFV = RDFV)
            b = CvANN(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = CvANN(ex,resultType = orange.GetBoth,returnDFV = RDFV)
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
            # CvANN does always return real probabilities on binary classification
            self.assertEqual(CvANN.isRealProb(),True)

        self.assertEqual([round(x,5) for x in CvANN.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(CvANN.nPredictions,(3+4)*len(self.LdataTest))

    def test_Probabilities(self):
        """Test if the returned probabilities are not fake"""

        CvANN = AZorngCvANN.CvANNLearner(self.LdataTrain,stopUPs=0)
        res = []
        for idx,ex in enumerate(self.LdataTest):
            res.append(CvANN(ex,resultType = orange.GetProbabilities))
            #print res[-1]
            self.assert_(res[-1][0]>=0 and res[-1][0]<=1,"Example "+str(idx)+" have impossible probability:"+str(res[-1])) 
            self.assert_(res[-1][1]>=0 and res[-1][1]<=1,"Example "+str(idx)+" have impossible probability:"+str(res[-1])) 
            self.assertEqual(CvANN.isRealProb(),True,"Example "+str(idx)+" did not return real probability")
            #print "Res",idx,":",res[-1]
            #print "Sum",idx,":",round(sum(res[-1]),5)
            self.assertEqual(round(sum(res[-1]),5),1,"Probabilities of Example "+str(idx)+" did not sum 1:"+str(res[-1]))
        sum0 = sum([x[0] for x in res])
        sum1 = sum([x[1] for x in res])
        self.assertEqual(len(self.LdataTest),round(sum0+sum1,5))
        self.assert_(sum0-int(sum0) > 0)
        self.assert_(sum1-int(sum1) > 0)


    def test_TwoWays(self):
        """
        Test that an ann created in one or two steps give the same results
        """

        # One step ann creation
        ann = AZorngCvANN.CvANNLearner(self.train_data,stopUPs=0)

        # Calculate classification accuracy for the classifier trained in one step
        oneStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, ann)

        # Two step ann creation
        learner = AZorngCvANN.CvANNLearner(randomWeights = False ,stopUPs=0)
        ann = learner(self.train_data)
        
        # Calculate classification accuracy for the classifier trained in two steps
        twoStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, ann) 

        # Test that the accuracy of the classifiers created in different ways is the exact same
        self.assertEqual(oneStepAcc, twoStepAcc)


    def test_SavedModel(self):
        """Test to assure that a saved ann model gives the same predictions as before saving."""

        # Create an ann model
        ann = AZorngCvANN.CvANNLearner(self.train_data,stopUPs=0)

        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.test_data, ann)

        # Save the model
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdir"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"ann.fann")
        ann.write(modelPath)
        
        # Read in the model
        ann = AZorngCvANN.CvANNread(modelPath)

        # Calculate classification accuracy 
        savedAcc = evalUtilities.getClassificationAccuracy(self.test_data, ann)

        # Test that the accuracy of the two classifiers is the exact same
        self.assertEqual(Acc, savedAcc)

        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)



    def test_ImputeTrain(self):
        """
        Assure that imputation works for the ann models. Test on data with missing values
        This test just assures the the model is trained. The correct imputation test is made on testImpute
        """
        annLearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)

        ann = annLearner(self.missingTrain)
    
        Acc = evalUtilities.getClassificationAccuracy(self.missingTest, ann)

        self.assertEqual(round(0.75758000000000003,5),round(Acc,5)) #opencv1.1: 0.95191999999999999


    def test_PredictionWithDiffVarType(self):
        """Test prediction with diff. VarType
        Test the prediction of examples with different varType
        """
        expectedAcc = 0.66666700000000001 # Ver 0.3 
        # Create a ann model
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)
        ann = CvANNlearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc2 = evalUtilities.getClassificationAccuracy(self.noBadDataTest[3:],ann)
        Acc1 = evalUtilities.getClassificationAccuracy(self.badVarTypeData[3:],ann)
        self.assertEqual(round(Acc1,6),round(expectedAcc,6))
        self.assertEqual(round(Acc2,6),round(expectedAcc,6))  
        self.assert_(('Fixed Types of variables' in ann.examplesFixedLog) and (ann.examplesFixedLog['Fixed Types of variables']==27), "No report of fixing in classifier class")
        self.assert_(('Vars needing type fix' in ann.examplesFixedLog) and (ann.examplesFixedLog['Vars needing type fix']['[Br]([C])']=="EnumVariable to FloatVariable", "No report of fixing in classifier class"))


    def test_PredictionWithDiffVarOrder(self):
        """Test Prediction with diff. VarOrder
        Test the prediction  examples with different varOrder
        """
        expectedAcc = 0.69999999999999996 # Ver 0.3
        # Create a ann model
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)
        ann = CvANNlearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,ann)
        Acc2 = evalUtilities.getClassificationAccuracy(self.badVarOrderData,ann)

        self.assertEqual(round(Acc1,9),round(expectedAcc,9),)
        self.assertEqual(round(Acc2,9),round(expectedAcc,9),)
        #we do not report order fix anymore!
        #self.assert_('Fixed Order of variables' in ann.examplesFixedLog and ann.examplesFixedLog['Fixed Order of variables']==27, "No report of fixing in classifier class")


    def test_PredictionWithIncompatibleDomain(self):
        """Test prediction with uncompatible domain
        Test the non-prediction of examples with an incompatible domain  
        """
        expectedAcc1 =  0.69999999999999996 # Ver 0.3 
        # Create a ann model
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)
        ann = CvANNlearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,ann)
        self.assertEqual(round(Acc1,9),round(expectedAcc1,9))
        self.assertEqual(ann(self.badVarTypeData[0]),"NEG","This example could still be predicted")
        self.assertEqual(ann(self.badVarTypeData[1]),"NEG","This example could still be predicted")
        self.assertEqual(ann(self.badVarNameData[0]),None,"This example should NOT be predicted")
        self.assertEqual(ann(self.badVarCountData[0]),None,"This example should NOT be predicted")


    def test_Impute(self):
        """Test missing values imputation
        Assure that imputation works for the ann models. Test on data with missing values
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
        
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)
        ann = CvANNlearner(contTrain)        

        # Prediction for data as it is
        P1=ann(ex1)
        P2=ann(ex2)
       
        # Predictions changing one continuous and one discrete variable to 0
        ex1["Desc 71"]=0
        ex2["Desc 138"]=0
        P1_0=ann(ex1)
        P2_0=ann(ex2)

        # Predictions changing the same continuous and discrete variable to it's correspondent imputation value
        ex1["Desc 71"]=imputer.defaults["Desc 71"]
        ex2["Desc 138"]=imputer.defaults["Desc 138"]
        P1_imp=ann(ex1)
        P2_imp=ann(ex2)
 
        # Predictions changing the same continuous and discrete variable to '?' wich means that the same imputation
        # as in the last case will have to be made inside the classifier. So, the predicted value must be the same
        ex1["Desc 71"]="?"
        ex2["Desc 138"]="?"
        self.assert_(ex1["Desc 71"]=="?","The var Desc 71 should be missing now!")
        self.assert_(ex2["Desc 138"]=="?","The var Desc 138 should be missing now!")    
        P1Miss=ann(ex1)
        P2Miss=ann(ex2)


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
        modelPath = os.path.join(scratchdir,"CvANNModel")
        ann.write(modelPath)

        # Read in the model
        annM = AZorngCvANN.CvANNread(modelPath)
        # Predict the ex1 and ex2 which are still the examples with missing values '?'
        self.assert_( ex1["Desc 71"]=="?","Value of Var Desc 71 should be missing!")
        self.assert_( ex2["Desc 138"]=="?","Value of Var Desc 138 should be missing!")
        self.assert_(round(annM(ex1),6)==round(P1Miss,6),"Imputation on loaded model is not correct")
        self.assert_(round(annM(ex2),6)==round(P2Miss,6),"Imputation on loaded model is not correct")
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)
        
    def test_MetaDataHandle(self):
        """Test the handling of Data with Meta Atributes
        """
        expectedAcc = 0.69999999999999996 # Ver 0.3
        # Create an ann model
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)
        ann = CvANNlearner(self.NoMetaTrain)

        # Calculate classification accuracy (NoMetaTest and WMeta are the same appart from the meta atribute) 
        AccNoMeta = evalUtilities.getClassificationAccuracy(self.NoMetaTest, ann)
        AccWMeta = evalUtilities.getClassificationAccuracy(self.WMetaTest, ann)
        self.assertEqual(AccNoMeta,AccWMeta,"Predictions with and without meta data were different!")
        self.assertEqual(round(AccNoMeta,9), round(expectedAcc,9))
        
    def test_MetaDataHandleForSavingModel(self):
        """Test the handling of SaveModel for Data with Meta Atributes
        """
        expectedAccWMeta = 1.0 # Ver 0.3 
        expectedAccNoMeta = 0.63333333300000005 # Ver 0.3
        #Test the save of a model created from a train data with meta attributes
        self.assert_(len(self.WMetaTest.domain.getmetas())>=1,"The dataset WMetaTest should have Meta Attributes")
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)
        annM = CvANNlearner(self.WMetaTest)
        AccNoMetaBefore = evalUtilities.getClassificationAccuracy(self.NoMetaTrain,annM) 
        AccWMetaBefore = evalUtilities.getClassificationAccuracy(self.WMetaTest,annM)


        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdiriTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"CvANNModel.CvANN")
        annM.write(modelPath)

        # Read in the model
        annR = AZorngCvANN.CvANNread(modelPath)
        self.assert_(len(annR.imputer.defaults.domain.getmetas())==0,"There shouldn't be any Meta data now!")

        # Calculate classification accuracy 
        AccNoMetaAfter = evalUtilities.getClassificationAccuracy(self.NoMetaTrain, annR)
        AccWMetaAfter = evalUtilities.getClassificationAccuracy(self.WMetaTest, annR)

        # Test that the accuracy of the model before and after saved
        self.assertEqual(AccNoMetaBefore, AccNoMetaAfter,"NoMeta: Predictions after loading saved model were different")
        self.assertEqual(AccWMetaBefore, AccWMetaAfter, "WMeta: Predictions after loading saved model were different")
        self.assertEqual(round(AccWMetaAfter,9), round(expectedAccWMeta,9))
        self.assertEqual(round(AccNoMetaAfter,9), round(expectedAccNoMeta,9))
 
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)

##ecPA

    def test_Scale(self):
        """
        Assure that scaling works. Also assures that the scale did not have the same resuts as in the test above (without scaling)
        """
       
        expectedAccWMeta = 1.0 # Ver 0.3
        expectedAccNoMeta = 0.63333333300000005 # Ver 0.3
        #Test the save of a model created from a train data with meta attributes
        self.assert_(len(self.WMetaTest.domain.getmetas())>=1,"The dataset WMetaTest should have Meta Attributes")
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100, scale = True,stopUPs=0)
        annM = CvANNlearner(self.WMetaTest)
        AccNoMetaBefore = evalUtilities.getClassificationAccuracy(self.NoMetaTrain,annM)
        AccWMetaBefore = evalUtilities.getClassificationAccuracy(self.WMetaTest,annM)


        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdiriTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"CvANNModel.CvANN")
        annM.write(modelPath)

        # Read in the model
        annR = AZorngCvANN.CvANNread(modelPath)
        self.assert_(len(annR.imputer.defaults.domain.getmetas())==0,"There shouldn't be any Meta data now!")

        # Calculate classification accuracy 
        AccNoMetaAfter = evalUtilities.getClassificationAccuracy(self.NoMetaTrain, annR)
        AccWMetaAfter = evalUtilities.getClassificationAccuracy(self.WMetaTest, annR)

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
        # One step ann creation
        ann = AZorngCvANN.CvANNLearner(self.train_data,nHidden = [3],stopUPs=0)
        # Calculate classification accuracy for the classifier trained in one step
        oneStepAcc = evalUtilities.getClassificationAccuracy(self.test_data, ann)
        # Check that the accuracy is what it used to be
        self.assertEqual(round(0.92381000000000002,5),round(oneStepAcc,5)) #orange1.0  0.95555999999999996,5


    def test_PersistentRegAcc(self): 
        """
        Assure that the accuracy is perserved for models trained in the same way. 
        """
        #This data is loaded here to speed up the test suite since it is too big
        contTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTest.tab")
        contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTrain.tab")
        contTrain = dataUtilities.DataTable(contTrainDataPath)
        contTest = dataUtilities.DataTable(contTestDataPath)

        # Create a CvANN model
        CvANNlearner = AZorngCvANN.CvANNLearner(randomWeights = False, nHidden = [3], nEpochs = 100,stopUPs=0)
        CvANNmodel = CvANNlearner(contTrain)
        # Calculate classification accuracy 
        Acc = evalUtilities.getRMSE(contTest, CvANNmodel)

        # Check that the accuracy is what it used to be
        self.assertEqual(round(0.109667,6),round(Acc,6))  #opencv1.1: 0.168131




if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(CvANNClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)



