from AZutilities import dataUtilities
import unittest
import os
import time

import orange
from AZutilities import evalUtilities
from AZutilities import miscUtilities
from trainingMethods import AZorngRF
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngImpute

class RFClassifierTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        """Creates the training and testing data set attributes. """
        trainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train.tab")
        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Test.tab")
        trainDataRegPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")
        testDataRegPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Test.tab")
        missingTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train_missing.tab")
        irisPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris.tab")
        # Read in the data
        missingInData = dataUtilities.DataTable(missingTestDataPath)
        self.trainData = dataUtilities.DataTable(trainDataPath)
        self.testData = dataUtilities.DataTable(testDataPath)
        self.trainDataReg = dataUtilities.DataTable(trainDataRegPath)
        self.testDataReg = dataUtilities.DataTable(testDataRegPath)
        self.irisData = dataUtilities.DataTable(irisPath)

        ##scPA
        dataNoMetaTrainPath = trainDataPath

        #These 2 datasets are equal apart from the meta atribute
        dataNoMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_SmallTest.tab")
        dataWMetaTestPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_W_metas_SmallTest.tab")

        # Read in the data   
        #contData = self.testDataReg 
        self.missingTrain = missingInData
        self.missingTest = missingInData
        self.NoMetaTrain = dataUtilities.DataTable(dataNoMetaTrainPath)
        self.NoMetaTest = dataUtilities.DataTable(dataNoMetaTestPath)
        self.WMetaTest = dataUtilities.DataTable(dataWMetaTestPath)

        #Data for domain fix handling        
        badVarTypePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarType.tab")
        badVarNamePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarName.tab")
        badVarOrderPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarOrder.tab")
        badVarCountPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarCount.tab")
        RegDAttrPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Imp_Train.tab")
        # Read in the data
        self.noBadDataTrain = self.NoMetaTrain
        self.noBadDataTest = self.NoMetaTest
        self.badVarTypeData = dataUtilities.DataTable(badVarTypePath)
        self.badVarNameData = dataUtilities.DataTable(badVarNamePath)
        self.badVarOrderData = dataUtilities.DataTable(badVarOrderPath)
        self.badVarCountData = dataUtilities.DataTable(badVarCountPath)   #One less example 
        self.RegDAttr = dataUtilities.DataTable(RegDAttrPath)
        # Random sampling
        #self.contTrain = self.trainDataReg
        #self.contTest = self.testDataReg
        ##ecPA
##scPA

    def test_PredictionVarImportance(self):
        """Test the variable importance for single predictions"""
        def TopVarImportanceTest(data, expectNone = False):
            resA = []
            resB = []
            RF = AZorngRF.RFLearner(data)

            for ex in data:
                resA.append(RF.getTopImportantVars(ex,1))

            scratchdir = miscUtilities.createScratchDir(desc="TopVarImportanceTest")
            modelPath = os.path.join(scratchdir,"CvRFModel")
            RF.write(modelPath)
            LoadedRF = AZorngRF.RFread(modelPath)
            miscUtilities.removeDir(scratchdir)

            for ex in data:
                resB.append(LoadedRF.getTopImportantVars(ex,1))
            if expectNone:
               return resA == resB == [None]*len(data)
            else:
                return resA == resB and None not in resA and resA.count(resA[0]) != len(resA)


        self.assertEqual(TopVarImportanceTest(self.irisData, True), True, "Failed testing with uncompatible data (class with 3 values)")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisCont.tab"))
        self.assertEqual(TopVarImportanceTest(data), True,"Failed testing with Regression data.")
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab"))
        self.assertEqual(TopVarImportanceTest(data), True, "Failed testing with Binary Classifier")


    def test_GlobalVarImportance(self):
        """Test the Global Variable Importance results"""        
        RF = AZorngRF.RFLearner(self.NoMetaTrain, getVarVariance = True)
        sum = 0
        for var in RF.domain.attributes:
            sum += RF.varImportance[var.name]
        self.assertEqual(round(sum,6),1.0)
                

    def test_DFV(self):
        """ Test the Decision Function Value Return"""
        RF = AZorngRF.RFLearner(self.NoMetaTrain)
        #Testsing with return of DFV
        RDFV = True
        for ex in self.NoMetaTest:
            predictedClass = RF(ex)
            a = RF(ex,returnDFV = RDFV)
            b = RF(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = RF(ex,resultType = orange.GetBoth,returnDFV = RDFV)
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
            # RF does always return real probabilities on binary classification
            self.assertEqual(RF.isRealProb(),True)
        expectedExtremes = {'max': 0.46000000000000002, 'min':-0.40999999999999998 } # Ver 0.3
        self.assertEqual([round(x,5) for x in RF.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(RF.nPredictions,4*len(self.NoMetaTest))

        #Testsing without return of DFV
        RDFV = False
        for ex in self.NoMetaTest:
            a = RF(ex,returnDFV = RDFV)
            b = RF(ex,resultType = orange.GetProbabilities,returnDFV = RDFV)
            c = RF(ex,resultType = orange.GetBoth,returnDFV = RDFV)
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
            # RF does always return real probabilities on binary classification
            self.assertEqual(RF.isRealProb(),True)

        self.assertEqual([round(x,5) for x in RF.getDFVExtremes().values()],[round(x,5) for x in expectedExtremes.values()])
        self.assertEqual(RF.nPredictions,(3+4)*len(self.NoMetaTest))

    def test_Probabilities(self):
        """Test if the returned probabilities are not fake"""

        RF = AZorngRF.RFLearner(self.trainData, nTrees = 200, nActVars = 155, maxDepth = 100)
        res = []
        for idx,ex in enumerate(self.testData):
            res.append(RF(ex,resultType = orange.GetProbabilities))
            self.assertEqual(RF.isRealProb(),True,"Example "+str(idx)+" did not return real probability")
            self.assert_(res[-1][0]>=0 and res[-1][0]<=1,"Example "+str(idx)+" have impossible probability:"+str(res[-1]))
            self.assert_(res[-1][1]>=0 and res[-1][1]<=1,"Example "+str(idx)+" have impossible probability:"+str(res[-1]))
            #print "Res",idx,":",res[-1]
            #print "Sum",idx,":",round(sum(res[-1]),5)
            self.assertEqual(round(sum(res[-1]),5),1,"Probabilities of Example "+str(idx)+" did not sum 1")
        sum0 = sum([x[0] for x in res])
        sum1 = sum([x[1] for x in res])
        self.assertEqual(len(self.testData),round(sum0+sum1,5))
        self.assert_(sum0-int(sum0) > 0)
        self.assert_(sum1-int(sum1) > 0)

    def test_save_load_Regression_D_Attr(self):
        """ Test Save/Load Regression model with Discrete Attribute"""

        #Create a selector to select just the correct attributes
        selector = range(len(self.RegDAttr.domain))
        #Remove the second attribute (idx=1)
        selector.pop(1)
        #Apply the selector to the self.RegDAttr
        data = self.RegDAttr.select(selector)

        RFsign = AZorngRF.RFLearner(data, nTrees = 200, nActVars = 155, maxDepth = 100)

        res1 = []
        for ex in self.RegDAttr:
            res1.append(str(RFsign(ex)))

        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"RFModel")
        RFsign.write(modelPath)

        loadedRFmodel = AZorngRF.RFread(modelPath)

        res2 = []
        for ex in self.RegDAttr:
            res2.append(str(loadedRFmodel(ex)))

        self.assertEqual(res1,res2)
        self.assertEqual(res1,['5.404782', '2.568249', '2.979486', '4.287185', '5.335753', '4.439877', '3.682451', '8.054751', '6.511803', '5.760388', '7.771009', '2.328262', '6.062288', '5.577081', '3.639579', '6.862591', '3.793468', '2.865258', '3.531777', '6.833398', '6.376686', '3.338588', '7.002612', '7.137580', '7.258987', '6.899173', '7.547265', '8.708020', '6.262212', '7.563741', '8.166364', '6.614120', '7.865033', '9.060866', '8.057292', '4.877943', '7.993115', '9.198319', '9.428467', '8.537990', '9.130789', '6.328936', '8.247712', '7.605743', '8.755456', '6.983065', '7.712387', '9.972745', '9.763152', '7.934700', '8.447981', '7.272462', '8.824869', '7.654151', '7.795481', '7.229007', '8.680950', '9.439033', '9.130064', '8.505672', '8.082146', '6.086042', '7.493593', '8.981513', '8.880632', '6.548739'])
   
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)



    def test_PredictionWithDiffVarType(self):
        """Test prediction with diff. VarType
        Test the prediction of examples with different varType
        """
        expectedAcc =0.96296296296296291 # Ver 0.3 
        # Create a rf model
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        rf = RFlearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc2 = evalUtilities.getClassificationAccuracy(self.noBadDataTest[3:],rf)
        Acc1 = evalUtilities.getClassificationAccuracy(self.badVarTypeData[3:],rf)
        self.assertEqual(Acc1,expectedAcc)
        self.assertEqual(Acc2,expectedAcc)   
        self.assert_(('Fixed Types of variables' in rf.examplesFixedLog) and (rf.examplesFixedLog['Fixed Types of variables']==27), "No report of fixing in classifier class")
        self.assert_(('Vars needing type fix' in rf.examplesFixedLog) and (rf.examplesFixedLog['Vars needing type fix']['[Br]([C])']=="EnumVariable to FloatVariable", "No report of fixing in classifier class"))


    def test_PredictionWithDiffVarOrder(self):
        """Test Prediction with diff. VarOrder
        Test the prediction  examples with different varOrder
        """
        expectedAcc = 0.96666666700000003 # ver 0.3
        # Create a rf model
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        rf = RFlearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,rf)
        Acc2 = evalUtilities.getClassificationAccuracy(self.badVarOrderData,rf)
        self.assertEqual(round(Acc1,9),round(expectedAcc,9))
        self.assertEqual(round(Acc2,9),round(expectedAcc,9))
        #we do not report order fix anymore!        
        #self.assert_('Fixed Order of variables' in rf.examplesFixedLog and rf.examplesFixedLog['Fixed Order of variables']==27, "No report of fixing in classifier class")

    def test_PredictionWithIncompatibleDomain(self):
        """Test prediction with uncompatible domain
        Test the non-prediction of examples with an incompatible domain  
        """
        expectedAcc1 = 0.96666666700000003 # Ver 0.3
        # Create a rf model
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        rf = RFlearner(self.noBadDataTrain)
        #using from index 3 o the end of data, because we know that from 0 to 2 the examples are not compatible
        Acc1 = evalUtilities.getClassificationAccuracy(self.noBadDataTest,rf)
        self.assertEqual(round(Acc1,9),round(expectedAcc1,9))
        self.assertEqual(rf(self.badVarTypeData[0]),"NEG","This example could still be predicted")
        self.assertEqual(rf(self.badVarTypeData[1]),"NEG","This example could still be predicted")
        self.assertEqual(rf(self.badVarNameData[0]),None,"This example should NOT be predicted")
        self.assertEqual(rf(self.badVarCountData[0]),None,"This example should NOT be predicted")
       
    def test_ImputeTrain(self):
        """
        Assure that imputation works for the rf models. Test on data with missing values
        This test just assures the the model is trained. The correct imputation test is made on testImpute
        """
        expected_Acc = [0.95757999999999999, 0.95455000000000001] #Ver 0.3 - Artifact: The second value can be expected on other Systems 
        rfLearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.001", termCrit = "0")
    
        rf = rfLearner(self.missingTrain)
    
        Acc = evalUtilities.getClassificationAccuracy(self.missingTest, rf)
        
        self.assert_(round(Acc,5) in [round(x,5) for x in expected_Acc])   #Ver 0.3  


    def test_Impute(self):
        """Test AZ missing values imputation
        Assure that imputation works for the rf models. Test on data with missing values
        """
        #This data is loaded here to speed up the test suite since it is too big
        contTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTest.tab")
        contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTrain.tab")
        contTrain = dataUtilities.DataTable(contTrainDataPath)
        contTest = dataUtilities.DataTable(contTestDataPath)

        ex1=contTest[5]
        ex2=contTest[6]
        self.assert_(ex1["Desc 71"]!="?","The var Desc 6 shouldn't be missing!")
        self.assert_(ex2["Desc 138"]!="?","The var Desc 71 shouldn't be missing!")

        imputer = orange.ImputerConstructor_average(contTrain)
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.001", termCrit = "0")
        rf = RFlearner(contTrain)

        # Prediction for data as it is
        P1=rf(ex1)
        P2=rf(ex2)
       
        # Predictions changing one continuous and one discrete variable to 0
        ex1["Desc 71"]=0
        ex2["Desc 138"]=0
        P1_0=rf(ex1)
        P2_0=rf(ex2)

        # Predictions changing the same continuous and discrete variable to it's correspondent imputation value
        ex1["Desc 71"]=imputer.defaults["Desc 71"]
        ex2["Desc 138"]=imputer.defaults["Desc 138"]
        P1_imp=rf(ex1)
        P2_imp=rf(ex2)
 
        # Predictions changing the same continuous and discrete variable to '?' wich means that the same imputation
        # as in the last case will have to be made inside the classifier. So, the predicted value must be the same
        ex1["Desc 71"]="?"
        ex2["Desc 138"]="?"
        self.assert_(ex1["Desc 71"]=="?","The var Desc 6 should be missing now!")
        self.assert_(ex2["Desc 138"]=="?","The var Desc 71 should be missing now!")    
        P1Miss=rf(ex1)
        P2Miss=rf(ex2)


        # Test if the prediction made for the example with mising value is the same as the one 
        # for the example which missing values were substituted using the same method as the classifier does.
        self.assert_(P1_imp==P1Miss,"Imputation was not made correctly inside the classifier")
        self.assert_(P2_imp==P2Miss,"Imputation was not made correctly inside the classifier")

        # Assure that if other substitutions on those variables were made, the predicted value would be different, 
        # and so, this is a valid method for testing the imputation
        self.assert_(P1.value!=P2.value)      # Just to assure that we are not comaring equal examples
        self.assert_(P1.value!=P1_imp.value,"The imputed 1 was the same as the original ... try other example")
        #opencv1.1:
        #self.assert_(P1_0.value!=P1_imp.value,"The imputed 1 was the same as the replaced by 0. The classifier may be replacing missing values by 0")
        #self.assert_(P2.value!=P2Miss.value, "The missing imputed 2 was the same as the original ... try other example")
        #self.assert_(P2_0.value!=P2Miss.value,"The missing imputed 2 was the same as the replaced by 0. The classifier may be replacing missing values by 0")


        #Test the imputer for saved models
        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"RFModel")
        rf.write(modelPath)

        # Read in the model
        rfM = AZorngRF.RFread(modelPath)
        # Predict the ex1 and ex2 which are still the examples with missing values '?'
        self.assert_( ex1["Desc 71"]=="?","Value of Var Desc 6 should be missing!")
        self.assert_( ex2["Desc 138"]=="?","Value of Var Desc 71 should be missing!")
        self.assert_(rfM(ex1)==P1Miss,"Imputation on loaded model is not correct")
        self.assert_(rfM(ex2)==P2Miss,"Imputation on loaded model is not correct")
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


    def test_BuiltIn_Impute(self):
        """Test RF BuiltIn missing values imputation
        Assure that imputation works for the rf models. Test on data with missing values
        """
        #This data is loaded here to speed up the test suite since it is too big
        contTestDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTest.tab")
        contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/linearTrain.tab")
        contTrain = dataUtilities.DataTable(contTrainDataPath)
        contTest = dataUtilities.DataTable(contTestDataPath)

        ex1=contTest[5]
        ex2=contTest[2]
        AttrEx1 = "Desc 71"
        AttrEx2 = "Desc 72"
        self.assert_(ex1[AttrEx1]!="?","The var Desc 671 shouldn't be missing!")
        self.assert_(ex2[AttrEx2]!="?","The var Desc 138 shouldn't be missing!")

        imputer = orange.ImputerConstructor_average(contTrain)
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.001", termCrit = "0",useBuiltInMissValHandling = True )
        rf = RFlearner(contTrain)

        # Prediction for data as it is
        P1=rf(ex1)
        P2=rf(ex2)
       
        # Predictions changing one continuous and one discrete variable to 0
        ex1[AttrEx1]=0
        ex2[AttrEx2]=0
        P1_0=rf(ex1)
        P2_0=rf(ex2)

        # Predictions changing the same continuous and discrete variable to it's correspondent imputation value
        #ex1["Desc 71"]=imputer.defaults["Desc 71"]
        #ex2["Desc 138"]=imputer.defaults["Desc 138"]
        #P1_imp=rf(ex1)
        #P2_imp=rf(ex2)
 
        # Predictions changing the same continuous and discrete variable to '?' wich means that the same imputation
        # as in the last case will have to be made inside the classifier. So, the predicted value must be the same
        ex1[AttrEx1]="?"
        ex2[AttrEx2]="?"
        self.assert_(ex1[AttrEx1]=="?","The var Desc 71 should be missing now!")
        self.assert_(ex2[AttrEx2]=="?","The var Desc 138 should be missing now!")    
        P1Miss=rf(ex1)
        P2Miss=rf(ex2)


        # Test if the prediction made for the example with mising value is the same as the one 
        # for the example which missing values were substituted using the same method as the classifier does.
        #self.assert_(P1_imp==P1Miss,"Imputation was not made correctly inside the classifier")
        #self.assert_(P2_imp==P2Miss,"Imputation was not made correctly inside the classifier")

        # Assure that if other substitutions on those variables were made, the predicted value would be different, 
        # and so, this is a valid method for testing the imputation


        self.assert_(P1.value!=P2.value)      # Just to assure that we are not comaring equal examples
        self.assert_(P1.value!=P1Miss.value,"The imputed 1 was the same as the original ... try other example")
        self.assert_(P1_0.value!=P1Miss.value,"The imputed 1 was the same as the replaced by 0. The classifier may be replacing missing values by 0")
        self.assert_(P2.value!=P2Miss.value, "The missing imputed 2 was the same as the original ... try other example")
        #self.assert_(P2_0.value!=P2Miss.value,"The missing imputed 2 was the same as the replaced by 0. The classifier may be replacing missing values by 0")

        self.assert_(rf.useBuiltInMissValHandling == True)
        #Test the imputer for saved models
        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"RFModel")
        rf.write(modelPath)

        # Read in the model
        rfM = AZorngRF.RFread(modelPath)
        self.assert_(rfM.useBuiltInMissValHandling == True)
        # Predict the ex1 and ex2 which are still the examples with missing values '?'
        self.assert_( ex1[AttrEx1]=="?","Value of Var Desc 6 should be missing!")
        self.assert_( ex2[AttrEx2]=="?","Value of Var Desc 71 should be missing!")
        self.assert_(rfM(ex1)==P1Miss,"Imputation on loaded model is not correct")
        self.assert_(rfM(ex2)==P2Miss,"Imputation on loaded model is not correct")
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


    def test_MetaDataHandle(self):
        """Test the handling of Data with Meta Atributes
        """
        expectedAcc = 0.96666666700000003 # Ver 0.3
        # Create an rf model
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        rf = RFlearner(self.NoMetaTrain)

        # Calculate classification accuracy (NoMetaTest and WMeta are the same appart from the meta atribute) 
        AccNoMeta = evalUtilities.getClassificationAccuracy(self.NoMetaTest, rf)
        AccWMeta = evalUtilities.getClassificationAccuracy(self.WMetaTest, rf)
        self.assertEqual(AccNoMeta,AccWMeta,"Predictions with and without meta data were different!")
        self.assertEqual(round(AccNoMeta,9), round(expectedAcc,9))
        
    def test_MetaDataHandleForSavingModel(self):
        """Test the handling of SaveModel for Data with Meta Atributes
        """
        expectedAccWMeta = 1.0 # VEr 0.3 
        expectedAccNoMeta = 0.56666666700000001 # Ver 0.3 
        #Test the save of a model created from a train data with meta attributes
        self.assert_(len(self.WMetaTest.domain.getmetas())>=1,"The dataset WMetaTest should have Meta Attributes")
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        rfM = RFlearner(self.WMetaTest)
        AccNoMetaBefore = evalUtilities.getClassificationAccuracy(self.NoMetaTrain,rfM) 
        AccWMetaBefore = evalUtilities.getClassificationAccuracy(self.WMetaTest,rfM)


        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"RFModel.RF")
        rfM.write(modelPath)

        # Read in the model
        rfR = AZorngRF.RFread(modelPath)
        self.assert_(len(rfR.domain.getmetas())==0,"There shouldn't be any Meta data now!")

        # Calculate classification accuracy 
        AccNoMetaAfter = evalUtilities.getClassificationAccuracy(self.NoMetaTrain, rfR)
        AccWMetaAfter = evalUtilities.getClassificationAccuracy(self.WMetaTest, rfR)

        # Test that the accuracy of the model before and after saved
        self.assertEqual(AccNoMetaBefore, AccNoMetaAfter,"NoMeta: Predictions after loading saved model were different")
        self.assertEqual(AccWMetaBefore, AccWMetaAfter, "WMeta: Predictions after loading saved model were different")
        self.assertEqual(round(AccWMetaAfter,9), round(expectedAccWMeta,9))
        self.assertEqual(round(AccNoMetaAfter,9), round(expectedAccNoMeta,9))
 
        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)
    
##ecPA
    def test_SavedModel(self):
        """Test to assure that a saved RF model gives the same predictions as before saving."""

        # Create a RF model
        RFlearner = AZorngRF.RFLearner(maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        RFmodel = RFlearner(self.trainData)

        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.testData, RFmodel)

        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"model.RF")
        RFmodel.write(modelPath)
        
        # Read in the model
        newRFmodel = AZorngRF.RFread(modelPath)

        # Calculate classification accuracy 
        savedAcc = evalUtilities.getClassificationAccuracy(self.testData, newRFmodel)

        # Test that the accuracy of the two classifiers is the exact same
        self.assertEqual(Acc, savedAcc)

        #Check the priors saved in the model
        file = open(os.path.join(modelPath,"model.RF"),"r")
        lines = file.readlines()
        file.close()
        priors = [round(x,2) for x in eval((lines[22].strip()).replace("data:",""))]
        self.assertEqual(len(priors),2)
        self.assertEqual(priors[self.testData.domain.classVar.values.index("POS")],0.50)
        self.assertEqual(priors[self.testData.domain.classVar.values.index("NEG")],0.50)

        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)


    def test_Priors(self):
        """Test to assure that priors are set correcly."""

        # Create a RF model
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0", priors = {"Iris-versicolor":0.35, "Iris-virginica":0.13, "Iris-setosa":0.52})
        RFmodel = RFlearner(self.irisData)

        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.irisData, RFmodel)

        # Save the model 
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirTest"+str(time.time()))
        os.mkdir(scratchdir)
        modelPath = os.path.join(scratchdir,"modelPriors.RF")
        RFmodel.write(modelPath)

        # Read in the model
        newRFmodel = AZorngRF.RFread(modelPath)

        # Calculate classification accuracy 
        savedAcc = evalUtilities.getClassificationAccuracy(self.irisData, newRFmodel)

        # Test that the accuracy of the two classifiers is the exact same
        self.assertEqual(Acc, savedAcc)

        #Check the priors saved in the model
        file = open(os.path.join(modelPath,"modelPriors.RF"),"r")
        lines = file.readlines()
        file.close()
        priors = [round(x,2) for x in eval((lines[22].strip()+lines[23].strip()).replace("data:",""))]
        self.assertEqual(len(priors),3)
        self.assertEqual(priors[self.irisData.domain.classVar.values.index("Iris-setosa")],0.52)
        self.assertEqual(priors[self.irisData.domain.classVar.values.index("Iris-versicolor")],0.35)
        self.assertEqual(priors[self.irisData.domain.classVar.values.index("Iris-virginica")],0.13)

        # Remove the scratch directory
        os.system("/bin/rm -rf "+scratchdir)



    def test_PersistentClassAcc(self): 
        """
        Assure that the accuracy is perserved for models trained in the same way. 
        """
        # Create a RF model
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        RFmodel = RFlearner(self.trainData)
        # Calculate classification accuracy 
        Acc = evalUtilities.getClassificationAccuracy(self.testData, RFmodel)
        # Check that the accuracy is what it used to be
        self.assertEqual(round(0.93869999999999998,5),round(Acc,5)) #opencv1.1: 0.77778000000000003

    def test_PersistentRegAcc(self): 
        """
        Assure that the accuracy is perserved for models trained in the same way. 
        """
        # Create a RF model
        RFlearner = AZorngRF.RFLearner(NumThreads = 1, maxDepth = "20", minSample = "5", useSurrogates = "false", getVarVariance = "false", \
                                        nActVars = "0", nTrees = "100", forestAcc = "0.1", termCrit = "0")
        RFmodel = RFlearner(self.trainDataReg)

        # Calculate classification accuracy 
        Acc = evalUtilities.getRMSE(self.testDataReg, RFmodel)

        # Check that the accuracy is what it used to be
        self.assertEqual(round(2.0158,5),round(Acc,5)) #opencv1.1:  0.32984999999999998,5


if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(RFClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

