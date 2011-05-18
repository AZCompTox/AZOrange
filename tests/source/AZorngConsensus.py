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
        learnersNames = ["CvANN","CvSVM","RF"]        

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
        self.assertEqual(Loaded.NTrainEx, len(self.irisData)-30)
        for ex in self.irisData:
            predictionsL.append(Loaded(ex))

        self.assertEqual(predictions,predictionsL)

        miscUtilities.removeDir(scratchdir) 


    def test_saveloadClass2(self):
        """Test the save/load for a classification model - Using probabilities average"""
        learnersNames = ["RF","CvANN"]

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
        self.assertEqual(Loaded.NTrainEx, len(self.irisData))

        miscUtilities.removeDir(scratchdir)

       
    def test_saveloadReg(self):
        """Test the save/load for a regression model - Using average of N classifiers"""
        learnersNames = ["CvANN","CvSVM","RF"]

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

        self.assertEqual([round(pred.value,4) for pred in predictions],[round(pred.value,4) for pred in predictionsL],"Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(self.DataReg.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(self.DataReg)-66)

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
        self.assertEqual([round(pred.value,4) for pred in predictions],[round(pred.value,4) for pred in predictionsL],"Loaded model predictions differ: Pred. 1 (saved/loaded):"+str(predictions[0])+" / "+str(predictionsL[0]))

        self.assertEqual(len(Loaded.domain),len(DataSet.domain))
        self.assertEqual(len(Loaded.imputeData) , len(Loaded.domain))
        self.assertEqual(len(Loaded.basicStat), len(Loaded.domain))
        self.assertEqual(Loaded.NTrainEx, len(DataSet))

        miscUtilities.removeDir(scratchdir)


    def teste_FeelLearnersClass(self):
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




if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(ConsensusClassifierTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

