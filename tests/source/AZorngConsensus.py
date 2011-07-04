from AZutilities import dataUtilities
import unittest
import os
import time
from trainingMethods import AZorngConsensus
from trainingMethods import AZorngCvSVM
from trainingMethods import AZorngCvANN
from trainingMethods import AZorngRF
from AZutilities import dataUtilities

import orange
from AZutilities import evalUtilities
from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngImpute


def float_compare(a, b):
    return abs(a-b)<0.0001

class ConsensusClassifierTest(AZorngTestUtil.AZorngTestUtil):
    def setUp(self):
        """Creates the training and testing data set attributes. """
        testDataRegPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")
        irisPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab")
        # Read in the data
        self.DataReg = dataUtilities.DataTable(testDataRegPath)
        self.irisData = dataUtilities.DataTable(irisPath)

    def getRegressionTrainingData(self):
        return self.DataReg

    def getClassificationTrainingData(self):
        return self.irisData

    def createTestLearners(self):
        return [AZorngRF.RFLearner(), AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner()]

    def test_CreateDefaultModel(self):
        """ Test the creation of Consensus Model using no learners """
        
        # Arrange

        # Act
        learner = AZorngConsensus.ConsensusLearner()

        # Assert
        self.assertEqual(learner.learners, None)
        self.assertEqual(learner.expression, None)
        self.assertEqual(learner.name, "Consensus learner")
        self.assertEqual(learner.verbose, 0)
        self.assertEqual(learner.imputeData, None)
        self.assertEqual(learner.NTrainEx, 0)
        self.assertEqual(learner.basicStat, None)
        self.assertEqual(learner.weights, None)


    def test_CreateModelWithLearnerList(self):
        """ Test the creation of Consensus Model using list of learners """
        
        # Arrange
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        
        # Act
        learner = AZorngConsensus.ConsensusLearner(learners = learners)

        # Assert
        for i, l in enumerate(learner.learners):
            self.assertEqual(learner.learners[i], learners[i])

        self.assertEqual(learner.expression, None)
        self.assertEqual(learner.name, "Consensus learner")
        self.assertEqual(learner.verbose, 0)
        self.assertEqual(learner.imputeData, None)
        self.assertEqual(learner.NTrainEx, 0)
        self.assertEqual(learner.basicStat, None)
        self.assertEqual(learner.weights, None)

    def test_CreateModelWithLearnerDictionary(self):
        """ Test the creation of Consensus Model using dictionary of learners """
        
        # Arrange
        learners = { 'a': AZorngCvSVM.CvSVMLearner(), 'b':AZorngCvANN.CvANNLearner(), 'c':AZorngRF.RFLearner() }
        expression = "a + b + c"
        
        # Act
        learner = AZorngConsensus.ConsensusLearner(learners = learners, expression = expression)

        # Assert
        for k,v in learner.learners.items():
            self.assertEqual(learner.learners[k], learners[k])
        self.assertEqual(learner.expression, expression)
        self.assertEqual(learner.name, "Consensus learner")
        self.assertEqual(learner.verbose, 0)
        self.assertEqual(learner.imputeData, None)
        self.assertEqual(learner.NTrainEx, 0)
        self.assertEqual(learner.basicStat, None)
        self.assertEqual(learner.weights, None)

    def test_CreateDefaultClassifierUsingTrainingData(self):
        """ Test the creation of default Classifier by calling learner with training data. """
        
        # Arrange
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        trainingData = self.getRegressionTrainingData()
        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        
        # Act
        classifier = learner(trainingData)

        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), len(learners))
        self.assertEqual(classifier.expression, None)
        self.assertEqual(classifier.name, "Consensus classifier")
        self.assertEqual(classifier.verbose, 0)
        self.assertNotEqual(classifier.imputeData, None)
        self.assertEqual(classifier.NTrainEx, len(trainingData))
        self.assertNotEqual(classifier.basicStat, None)
        self.assertEqual(classifier.weights, None)

    def test_CreateCustomRegressionClassifierUsingTrainingData(self):
        """ Test the creation of custom regression Classifier by calling learner with training data. """
        
        # Arrange
        learners = { 'a': AZorngCvSVM.CvSVMLearner(), 'b':AZorngCvANN.CvANNLearner(), 'c':AZorngRF.RFLearner() }
        expression = "a + b + c"
        trainingData = self.getRegressionTrainingData()
        learner = AZorngConsensus.ConsensusLearner(learners = learners, expression = expression)
        
        # Act
        classifier = learner(trainingData)

        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), len(learners))
        self.assertNotEqual(classifier.basicStat, None)
        self.assertNotEqual(classifier.classVar, None)
        self.assertNotEqual(classifier.domain, None)
        self.assertEqual(classifier.expression, expression)
        self.assertNotEqual(classifier.imputeData, None)
        self.assertEqual(classifier.NTrainEx, len(trainingData))
        self.assertEqual(classifier.name, "Consensus classifier")
        self.assertNotEqual(classifier.varNames, None)
        self.assertEqual(classifier.verbose, 0)
        self.assertEqual(classifier.weights, None)

    def test_CreateCustomClassificationClassifierUsingTrainingData(self):
        """ Test the creation of custom classification Classifier by calling learner with training data. """
        
        # Arrange
        learners = { 'a': AZorngCvSVM.CvSVMLearner(), 'b':AZorngCvANN.CvANNLearner(), 'c':AZorngRF.RFLearner() }
        expression = ["firstLearner == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        trainingData = self.getClassificationTrainingData()
        learner = AZorngConsensus.ConsensusLearner(learners = learners, expression = expression)
        
        # Act
        classifier = learner(trainingData)

        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), len(learners))
        self.assertNotEqual(classifier.basicStat, None)
        self.assertNotEqual(classifier.classVar, None)
        self.assertNotEqual(classifier.domain, None)
        self.assertEqual(classifier.expression, expression)
        self.assertNotEqual(classifier.imputeData, None)
        self.assertEqual(classifier.NTrainEx, len(trainingData))
        self.assertEqual(classifier.name, "Consensus classifier")
        self.assertNotEqual(classifier.varNames, None)
        self.assertEqual(classifier.verbose, 0)
        self.assertEqual(classifier.weights, None)


    def test_CreateDefaultClassifierUsingPreTrainedClassificationClassifiers(self):
        """ Test the creation of default Consensus Classifier using pre-trained classification classifiers. """

        # Arrange
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        classifiers = [l(self.getClassificationTrainingData()) for l in learners]

        # Act
        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers)

        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), len(learners))
        self.assertNotEqual(classifier.basicStat, None)
        self.assertNotEqual(classifier.classVar, None)
        self.assertNotEqual(classifier.domain, None)
        self.assertEqual(classifier.expression, None)
        self.assertNotEqual(classifier.imputeData, None)
        #self.assertEqual(classifier.NTrainEx, len(trainingData))
        self.assertEqual(classifier.name, "Consensus classifier")
        self.assertNotEqual(classifier.varNames, None)
        self.assertEqual(classifier.verbose, 0)
        self.assertEqual(classifier.weights, None)


    def test_CreateDefaultClassifierUsingPreTrainedRegressionClassifiers(self):
        """ Test the creation of default Consensus Classifier using pre-trained classification classifiers. """

        # Arrange
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        classifiers = [l(self.getRegressionTrainingData()) for l in learners]

        # Act
        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers)

        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), len(learners))
        self.assertNotEqual(classifier.basicStat, None)
        self.assertNotEqual(classifier.classVar, None)
        self.assertNotEqual(classifier.domain, None)
        self.assertEqual(classifier.expression, None)
        self.assertNotEqual(classifier.imputeData, None)
        #self.assertEqual(classifier.NTrainEx, len(trainingData))
        self.assertEqual(classifier.name, "Consensus classifier")
        self.assertNotEqual(classifier.varNames, None)
        self.assertEqual(classifier.verbose, 0)
        self.assertEqual(classifier.weights, None)


    def test_CreateCustomClassifierUsingPreTrainedClassificationClassifiers(self):
        """ Test the creation of custom Consensus Classifier using pre-trained classification classifiers. """

        # Arrange
        learners = { 'a':AZorngCvSVM.CvSVMLearner(), 'b':AZorngCvANN.CvANNLearner(), 'c':AZorngRF.RFLearner() }
        classifiers = {}
        for k,v in learners.items():
            classifiers[k] = v(self.getClassificationTrainingData())
        expression = ["firstLearner == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]

        # Act
        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers, expression = expression)

        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), len(classifiers))
        self.assertNotEqual(classifier.basicStat, None)
        self.assertNotEqual(classifier.classVar, None)
        self.assertNotEqual(classifier.domain, None)
        self.assertEqual(classifier.expression, expression)
        self.assertNotEqual(classifier.imputeData, None)
        #self.assertEqual(classifier.NTrainEx, len(trainingData))
        self.assertEqual(classifier.name, "Consensus classifier")
        self.assertNotEqual(classifier.varNames, None)
        self.assertEqual(classifier.verbose, 0)
        self.assertEqual(classifier.weights, None)


    def test_CreateDefaultClassifierUsingPreTrainedRegressionClassifiers(self):
        """ Test the creation of custom Consensus Classifier using pre-trained regression classifiers. """

        # Arrange
        learners = { 'a':AZorngCvSVM.CvSVMLearner(), 'b':AZorngCvANN.CvANNLearner(), 'c':AZorngRF.RFLearner() }
        classifiers = {}
        for k,v in learners.items():
            classifiers[k] = v(self.getRegressionTrainingData())
        expression = "a + b + c"

        # Act
        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers, expression = expression)

        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), len(learners))
        self.assertNotEqual(classifier.basicStat, None)
        self.assertNotEqual(classifier.classVar, None)
        self.assertNotEqual(classifier.domain, None)
        self.assertEqual(classifier.expression, expression)
        self.assertNotEqual(classifier.imputeData, None)
        #self.assertEqual(classifier.NTrainEx, len(trainingData))
        self.assertEqual(classifier.name, "Consensus classifier")
        self.assertNotEqual(classifier.varNames, None)
        self.assertEqual(classifier.verbose, 0)
        self.assertEqual(classifier.weights, None)


    def test_CanPersistClassificationModelMajority(self):
        """Test the save/load for a classification model - Using Majority"""

        """ Arrange """
        learners = self.createTestLearners()
        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.getClassificationTrainingData())

        """ Act """
        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        """ Assert """
        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain)) 
        self.assertEqual(Loaded.NTrainEx, 0.8*len(self.irisData))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        miscUtilities.removeDir(scratchdir) 


    def test_CanPersistClassificationModelProbabilities(self):
        """Test the save/load for a classification model - Using probabilities average"""

        # Arrange
        learners = [AZorngRF.RFLearner(), AZorngCvANN.CvANNLearner()]
        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.irisData)

        # Act
        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        # Assert
        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)
        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.irisData) - int(0.2 * len(self.irisData)))

        miscUtilities.removeDir(scratchdir)

       
    def test_CanPersistRegressionModelUsingClassifiers(self):
        """Test the save/load for a regression model - Using average of N classifiers"""

        # Arrange
        learners = [AZorngRF.RFLearner(), AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner()]
        learner = AZorngConsensus.ConsensusLearner(learners = learners)
        classifier = learner(self.DataReg)

        # Act
        predictions = []
        for ex in self.DataReg:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        # Assert
        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.DataReg:
            predictionsL.append(Loaded(ex))

        self.assertEqual([round(pred.value,4) for pred in predictions],
                         [round(pred.value,4) for pred in predictionsL],
                         "Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(self.DataReg.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.DataReg)*0.8)

        miscUtilities.removeDir(scratchdir)
 
        
    def test_CreateLearnerWithObjectMapping(self):
        """ Test the creation of learners with an object map """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}

        # Act
        learner = AZorngConsensus.ConsensusLearner(learners = learners)

        # Assert
        self.assertEqual(len(learner.learners), len(learners))

    def test_CreateLearnerWithObjectMappingWithoutExpression(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        learner = AZorngConsensus.ConsensusLearner(learners = learners)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertEqual(classifier, None)


    def test_CanCreateClassifierUsingObjMapping(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        discreteExpression = ""
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 2"
        
        learner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), 3)
        self.assertEqual(classifier.expression, regressionExpression)


    def test_AverageNRegressionExpressionUsingObjMap(self):
        """ Test regular expression using average N regression with object map """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        # Construct expression learner/classifier
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression)
        expressionClassifier = expressionLearner(self.DataReg)

        # Construct default learner/classifier
        defaultLearners = [AZorngRF.RFLearner(), AZorngCvANN.CvANNLearner(), AZorngCvSVM.CvSVMLearner()]
        defaultLearner = AZorngConsensus.ConsensusLearner(learners = defaultLearners)
        defaultClassifier = defaultLearner(self.DataReg)
        
        # Act
        expressionPredictions = []
        for ex in self.DataReg:
            expressionPredictions.append(expressionClassifier(ex))

        defaultPredictions = []
        for ex in self.DataReg:
            defaultPredictions.append(defaultClassifier(ex))

        # Assert
        for index in range(len(expressionPredictions)):
            self.assertEqual(True, float_compare(expressionPredictions[index], defaultPredictions[index]))

    def test_CreateLogicalExpressionConsensusLearner(self):
        """ Test creation of logical expression consensus learner """
        # Arrange
        
        # Construct expression learner/classifier
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        discreteExpression = ["firstLearner == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)

        verifiedLearner = AZorngCvSVM.CvSVMLearner()
        verifiedClassifier = verifiedLearner(self.irisData)
        
        # Act
        result = []
        verifiedResult = []
        for ex in self.irisData:
            result.append(discreteClassifier(ex))
            verifiedResult.append(verifiedClassifier(ex))

        # Assert
        for index, item in enumerate(result):
            if not result[index].value == verifiedResult[index].value:
                print "Not equal on index: ", index
            self.assertEqual(result[index].value, verifiedResult[index].value)

    def test_SaveLoadCustomLogicalExpression(self):
        """ Test save/load functionality with a custom logical expression """
        # Arrange
        
        # Construct expression learner/classifier
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        discreteExpression = ["firstLearner == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)

        result = []
        for ex in self.irisData:
            result.append(discreteClassifier(ex))
        
        # Act
        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        discreteClassifier.write(os.path.join(scratchdir,"./CM.model"))

        resultLoaded = []
        loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        self.assertNotEqual(loaded, None)
        for ex in self.irisData:
            resultLoaded.append(loaded(ex))

        # Assert
        for index, item in enumerate(result):
            if not result[index].value == resultLoaded[index].value:
                print "Not equal on index: ", index
            self.assertEqual(result[index].value, resultLoaded[index].value)

        self.assertEqual(len(loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(loaded.imputeData) , len(loaded.domain))
        self.assertEqual(len(loaded.basicStat), len(loaded.domain))
        self.assertEqual(loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)

    def test_SaveLoadCustomRegressionExpression(self):
        """ Test save/load custom expression using average N regression with object map """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        # Construct expression learner/classifier
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression)
        expressionClassifier = expressionLearner(self.DataReg)

        # Construct default learner/classifier
        result = []
        for ex in self.DataReg:
            result.append(expressionClassifier(ex))

        # Act
        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        expressionClassifier.write(os.path.join(scratchdir,"./CM.model"))

        resultLoaded = []
        loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        self.assertNotEqual(loaded, None)
        for ex in self.DataReg:
            resultLoaded.append(loaded(ex))

        # Assert
        for index, item in enumerate(result):
            if not float_compare(result[index].value, resultLoaded[index].value):
                print "Not equal on index: ", index
            self.assertEqual(float_compare(result[index].value, resultLoaded[index].value), True)

        self.assertEqual(len(loaded.domain),len(self.DataReg.domain))
        self.assertEqual(len(loaded.imputeData) , len(loaded.domain))
        self.assertEqual(len(loaded.basicStat), len(loaded.domain))
        self.assertEqual(loaded.NTrainEx, len(self.DataReg))

        miscUtilities.removeDir(scratchdir)


    def test_CustomLogicalExpressionUsingAndStatement(self):
        """ Test logical expression using AND statements """
        # Arrange

        # Construct verification learners
        a = AZorngCvSVM.CvSVMLearner()
        a = a(self.irisData)
        b = AZorngCvANN.CvANNLearner()
        b = b(self.irisData)
        c = AZorngRF.RFLearner()
        c = c(self.irisData)
        
        # Construct expression learner/classifier
        learners = {'a':AZorngCvSVM.CvSVMLearner(),
                    'b':AZorngCvANN.CvANNLearner(),
                    'c':AZorngRF.RFLearner()}
        discreteExpression = ["a == Iris-setosa and b == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)

        # Act
        result = []
        for ex in self.irisData:
            result.append(discreteClassifier(ex))

        verifiedResult = []
        for ex in self.irisData:
            if a(ex).value == "Iris-setosa" and b(ex).value == "Iris-setosa":
                verifiedResult.append("Iris-setosa")
            else:
                verifiedResult.append("Iris-virginica")


        # Assert
        for index, item in enumerate(result):
            if not result[index].value == verifiedResult[index]:
                print "Not equal on index: ", index, " Predicted: ", result[index].value, " Real: ", verifiedResult[index]
            self.assertEqual(result[index].value, verifiedResult[index])


    def test_CustomLogicalExpressionUsingOrStatement(self):
        """ Test logical expression using OR statements """
        # Arrange

        # Construct verification learners
        a = AZorngCvSVM.CvSVMLearner()
        a = a(self.irisData)
        b = AZorngCvANN.CvANNLearner()
        b = b(self.irisData)
        c = AZorngRF.RFLearner()
        c = c(self.irisData)
        
        # Construct expression learner/classifier
        learners = {'a':AZorngCvSVM.CvSVMLearner(),
                    'b':AZorngCvANN.CvANNLearner(),
                    'c':AZorngRF.RFLearner()}
        discreteExpression = ["a == Iris-setosa or b == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)

        # Act
        result = []
        for ex in self.irisData:
            result.append(discreteClassifier(ex))

        verifiedResult = []
        for ex in self.irisData:
            if a(ex).value == "Iris-setosa" or b(ex).value == "Iris-setosa":
                verifiedResult.append("Iris-setosa")
            else:
                verifiedResult.append("Iris-virginica")


        # Assert
        for index, item in enumerate(result):
            if not result[index].value == verifiedResult[index]:
                print "Not equal on index: ", index, " Predicted: ", result[index].value, " Real: ", verifiedResult[index]
            self.assertEqual(result[index].value, verifiedResult[index])

    def test_CustomLogicalExpressionUsingOrAndStatement(self):
        """ Test logical expression using OR/AND statements """
        # Arrange

        # Construct verification learners
        a = AZorngCvSVM.CvSVMLearner()
        a = a(self.irisData)
        b = AZorngCvANN.CvANNLearner()
        b = b(self.irisData)
        c = AZorngRF.RFLearner()
        c = c(self.irisData)
        
        # Construct expression learner/classifier
        learners = {'a':AZorngCvSVM.CvSVMLearner(),
                    'b':AZorngCvANN.CvANNLearner(),
                    'c':AZorngRF.RFLearner()}
        discreteExpression = ["a == Iris-setosa and c == Iris-virginica or b == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)

        # Act
        result = []
        for ex in self.irisData:
            result.append(discreteClassifier(ex))

        verifiedResult = []
        for ex in self.irisData:
            if a(ex).value == "Iris-setosa" and c(ex).value == "Iris-virginica" or b(ex).value == "Iris-setosa":
                verifiedResult.append("Iris-setosa")
            else:
                verifiedResult.append("Iris-virginica")


        # Assert
        for index, item in enumerate(result):
            if not result[index].value == verifiedResult[index]:
                print "Not equal on index: ", index, " Predicted: ", result[index].value, " Real: ", verifiedResult[index]
            self.assertEqual(result[index].value, verifiedResult[index])

    def test_CustomRegressionExpressionUsingWeights(self):
        """ Test regression expression using weights """
        # Arrange
        learners = {'a':AZorngCvSVM.CvSVMLearner(),
                    'b':AZorngCvANN.CvANNLearner(),
                    'c':AZorngRF.RFLearner()}

        weights = { 'a': lambda x: 1,
                    'b': lambda x: 2,
                    'c': lambda x: 3 }
        
        regressionExpression = "(a + b + c) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression, weights = weights)
        classifier = expressionLearner(self.DataReg)

        # Act
        result = []
        for ex in self.DataReg:
            result.append(classifier(ex))

        verifiedResult = []
        for ex in self.DataReg:
            a_value = classifier.classifiers['a'](ex)
            a_weight_value = weights['a'](a_value)
            b_value = classifier.classifiers['b'](ex)
            b_weight_value = weights['b'](b_value)
            c_value = classifier.classifiers['c'](ex)
            c_weight_value = weights['c'](c_value) 

            prediction = (a_value * a_weight_value +
                          b_value * b_weight_value +
                          c_value * c_weight_value)/3
            
            verifiedResult.append(prediction)
            
        # Assert
        for index, item in enumerate(result):
            if float_compare(result[index].value, verifiedResult[index]) == False:
                print "Not equal on index: ", index
                print "Result: ", result[index].value, " Verified: ", verifiedResult[index]
                print "Delta: ", abs(result[index].value - verifiedResult[index])
            self.assertEqual(float_compare(result[index].value, verifiedResult[index]), True)

    def test_InvalidCustomRegressionExpression(self):
        """ Test invalid custom expression """
        # Arrange
        learners = {'a':AZorngCvSVM.CvSVMLearner(),
                    'b':AZorngCvANN.CvANNLearner(),
                    'c':AZorngRF.RFLearner()}

        regressionExpression = "(a + b + 3cd45 + c) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = regressionExpression)

        # Act
        classifier = expressionLearner(self.DataReg)

        # Assert
        self.assertEqual(classifier(self.DataReg[0]), None)

    def test_InvalidCustomClassificationExpression(self):
        """ Test invalid custom expression """
        # Arrange
        learners = {'a':AZorngCvSVM.CvSVMLearner(),
                    'b':AZorngCvANN.CvANNLearner(),
                    'c':AZorngRF.RFLearner()}

        discreteExpression = ["a == Iris-setosa and or c == Iris-virginica or b == Iris-setosa -> Iris-setosa", "-> Iris-virginica"]
        expressionLearner = AZorngConsensus.ConsensusLearner(learners = learners, expression = discreteExpression)

        # Act
        classifier = expressionLearner(self.getClassificationTrainingData())

        # Assert
        self.assertEqual(classifier(self.getClassificationTrainingData()[0]), None)



if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(ConsensusClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

