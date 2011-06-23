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
    return abs(a-b)<0.00001

class ConsensusClassifierTest(AZorngTestUtil.AZorngTestUtil):
    def setUp(self):
        """Creates the training and testing data set attributes. """
        testDataRegPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")
        irisPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab")
        # Read in the data
        self.DataReg = dataUtilities.DataTable(testDataRegPath)
        self.irisData = dataUtilities.DataTable(irisPath)

    def test_saveloadClass(self):
        """Test the save/load for a classification model - Using Majority"""
        learnersNames = ["RF","CvANN","CvSVM"]        

        learner = AZorngConsensus.ConsensusLearner(learnersNames = learnersNames)
        classifier = learner(self.irisData)

        self.assertEqual(len(classifier.domain),len(self.irisData.domain))
        self.assertEqual(len(classifier.imputeData) , len(classifier.domain))
        self.assertEqual(len(classifier.basicStat), len(classifier.domain))
        self.assertEqual(classifier.NTrainEx, len(self.irisData))
        

        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))
        
        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain)) 
        self.assertEqual(Loaded.NTrainEx, len(self.irisData))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        miscUtilities.removeDir(scratchdir) 


    def test_saveloadClass2(self):
        """Test the save/load for a classification model - Using probabilities average"""
        learnersNames = ["CvANN","RF"]

        learner = AZorngConsensus.ConsensusLearner(learnersNames = learnersNames)
        classifier = learner(self.irisData)
        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)
        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.irisData)-int(0.2 * len(self.irisData)))

        miscUtilities.removeDir(scratchdir)

       
    def test_saveloadReg(self):
        """Test the save/load for a regression model - Using average of N classifiers"""
        learnersNames = ["RF","CvSVM","CvANN"]

        learner = AZorngConsensus.ConsensusLearner(learnersNames = learnersNames)
        classifier = learner(self.DataReg)
        predictions = []
        for ex in self.DataReg:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

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
        self.assertEqual(Loaded.NTrainEx, len(self.DataReg))

        miscUtilities.removeDir(scratchdir)
 

    def test_FeedLearnersReg(self):
        """Test the creation of Consensus feeding Learners for regression"""
        #The Learners can be individualy costumized before passing them to the Consensus
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]

        #Passing now the learnersObj instead
        learner = AZorngConsensus.ConsensusLearner(learnersObj = learners)
        classifier = learner(self.DataReg)
        predictions = []
        for ex in self.DataReg:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.DataReg:
            predictionsL.append(Loaded(ex))

        self.assertEqual([round(pred.value,4) for pred in predictions],[round(pred.value,4) for pred in predictionsL],"Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(self.DataReg.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.DataReg))

        miscUtilities.removeDir(scratchdir)


    def test_FeedClassifiersClass(self):
        """Test the creation of Consensus feeding Classifiers"""

        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        classifiers = [l(self.irisData) for l in learners]

        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers)

        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)


    def test_FeedClassifiersReg(self):
        """Test the feeding of regression classifiers """
        #DataSet = dataUtilities.DataTable("/home/palmeida/dev/OpenAZOTesteInstall/tests/source/data/linearTrain.tab")
        DataSet = self.DataReg
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]
        classifiers = [l(DataSet) for l in learners]

        classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers)
        predictions = []
        for ex in DataSet:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in DataSet:
            predictionsL.append(Loaded(ex))
        self.assertEqual([round(pred.value,4) for pred in predictions],
                         [round(pred.value,4) for pred in predictionsL],
                         "Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(DataSet.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(DataSet))

        miscUtilities.removeDir(scratchdir)


    def test_FeedLearnersClass(self):
        """Test the creation of Consensus feeding Learners for classification"""
        #The Learners can be individualy costumized before passing them to the Consensus
        learners = [AZorngCvSVM.CvSVMLearner(), AZorngCvANN.CvANNLearner(), AZorngRF.RFLearner()]

        #Passing now the learnersObj instead
        learner = AZorngConsensus.ConsensusLearner(learnersObj = learners)
        classifier = learner(self.irisData)
        predictions = []
        for ex in self.irisData:
            predictions.append(classifier(ex))

        scratchdir = miscUtilities.createScratchDir(desc="ConsensusSaveLoadTest")
        classifier.write(os.path.join(scratchdir,"./CM.model"))

        predictionsL = []
        Loaded = AZorngConsensus.Consensusread(os.path.join(scratchdir,"./CM.model"))
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        self.assertEqual(len(Loaded.domain),len(self.irisData.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)

        
    def test_CreateLearnerWithObjectMapping(self):
        """ Test the creation of learners with an object map """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}

        # Act
        learner = AZorngConsensus.ConsensusLearner(learnerObjMap = learners)

        # Assert
        self.assertEqual(len(learner.learnerObjMap), len(learners))

    def test_CreateLearnerWithoutLearnersDefined(self):
        """ Test ConsensusLearner without required constructor arguments """

        # Test with no learner names or objects defined
        # Arrange
        learners = {'firstLearner':'CvSVM',
                    'secondLearner':'CvANN',
                    'thirdLearner':'RF'}

        # Act
        learner = AZorngConsensus.ConsensusLearner()
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertEqual(classifier, None)
        

    def test_CreateRegressionLearnerWithoutVariableMapping(self):
        """ Test with regression expression defined, but not variable mapping """
        # Arrange
        learner = AZorngConsensus.ConsensusLearner(regressionExpression = "testing")

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertEqual(classifier, None)

        
    def test_CreateDiscreteLearnerWithoutVariableMapping(self):
        """ Test with discrete expression defined, but not variable mapping """
        # Arrange
        learner = AZorngConsensus.ConsensusLearner(discreteExpression = "testing")

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertEqual(classifier, None)


    def test_CreateLearnerWithNameMappingWithoutExpression(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        learners = {'firstLearner':'CvSVM',
                    'secondLearner':'CvANN',
                    'thirdLearner':'RF'}
        learner = AZorngConsensus.ConsensusLearner(learnerNameMap = learners)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertEqual(classifier, None)


    def test_CreateLearnerWithObjectMappingWithoutExpression(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        learners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        learner = AZorngConsensus.ConsensusLearner(learnerObjMap = learners)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertEqual(classifier, None)


    def test_CanCreateClassifierUsingNameMapping(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        objLearners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        nameLearners = {'firstLearner':'CvSVM',
                    'secondLearner':'CvANN',
                    'thirdLearner':'RF'}
        discreteExpression = ""
        regressionExpression = "(firstLearner + secondLearner) / 2"
        
        learner = AZorngConsensus.ConsensusLearner(learnerNameMap = nameLearners, regressionExpression = regressionExpression)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), 3)
        self.assertEqual(classifier.regressionExpression, regressionExpression)

    def test_CanCreateClassifierUsingObjMapping(self):
        """ Test with name variable mapping defined but not expression given """
        # Arrange
        objLearners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        
        nameLearners = {'firstLearner':'CvSVM',
                    'secondLearner':'CvANN',
                    'thirdLearner':'RF'}
        discreteExpression = ""
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 2"
        
        learner = AZorngConsensus.ConsensusLearner(learnerObjMap = objLearners, regressionExpression = regressionExpression)

        # Act
        classifier = learner(self.DataReg)
        
        # Assert
        self.assertNotEqual(classifier, None)
        self.assertEqual(len(classifier.classifiers), 3)
        self.assertEqual(classifier.regressionExpression, regressionExpression)

    def test_AverageNRegressionExpressionUsingObjMap(self):
        """ Test regular expression using average N regression with object map """
        # Arrange

        # Construct expression learner/classifier
        objLearners = {'firstLearner':AZorngCvSVM.CvSVMLearner(),
                    'secondLearner':AZorngCvANN.CvANNLearner(),
                    'thirdLearner':AZorngRF.RFLearner()}
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learnerObjMap = objLearners, regressionExpression = regressionExpression)
        expressionClassifier = expressionLearner(self.DataReg)

        # Construct default learner/classifier
        learnersNames = ["RF","CvANN","CvSVM"]
        defaultLearner = AZorngConsensus.ConsensusLearner(learnersNames = learnersNames)
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

    def test_AverageNRegressionExpressionUsingNameMap(self):
        """ Test regular expression using average N regression with name map """
        # Arrange

        # Construct expression learner/classifier
        nameLearners = {'firstLearner':'CvSVM',
                    'secondLearner':'CvANN',
                    'thirdLearner':'RF'}
        regressionExpression = "(firstLearner + secondLearner + thirdLearner) / 3"
        expressionLearner = AZorngConsensus.ConsensusLearner(learnerNameMap = nameLearners, regressionExpression = regressionExpression)
        expressionClassifier = expressionLearner(self.DataReg)

        # Construct default learner/classifier
        learnersNames = ["RF","CvANN","CvSVM"]
        defaultLearner = AZorngConsensus.ConsensusLearner(learnersNames = learnersNames)
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
        nameLearners = {'firstLearner':'CvSVM',
                        'secondLearner':'CvANN',
                        'thirdLearner':'RF'}
        discreteExpression = ["firstLearner == Iris-setosa -> Iris-setosa"]
        discreteLearner = AZorngConsensus.ConsensusLearner(learnerNameMap = nameLearners, discreteExpression = discreteExpression)
        discreteClassifier = discreteLearner(self.irisData)
        
        # Act
        result = discreteClassifier(self.irisData[0])

        # Assert
        self.assertEqual(result.value, 'Iris-setosa')

if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(ConsensusClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

