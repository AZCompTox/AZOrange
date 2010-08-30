import unittest
import os
import time

import orange
import orngTest
from AZutilities import dataUtilities
from AZutilities import evalUtilities
from AZutilities import similarityMetrics
from trainingMethods import AZorngRF
import AZOrangeConfig as AZOC
import AZorngTestUtil


class evalUtilitiesTest(unittest.TestCase):

    def setUp(self):

        trainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train.tab")
        self.trainData = dataUtilities.DataTable(trainDataPath,createNewOn=orange.Variable.MakeStatus.OK)

        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Test.tab")
        self.testData = dataUtilities.DataTable(testDataPath,createNewOn=orange.Variable.MakeStatus.OK)

        self.regDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")

    def testQuantileCalc(self):

        MD = similarityMetrics.calcMahalanobis(self.trainData, self.testData)
        quantiles = similarityMetrics.calcMahalanobisDistanceQuantiles(MD)

        self.assertEqual(round(quantiles[0],3), round(1.216,3))
        self.assertEqual(round(quantiles[1],3), round(1.717,3))
        self.assertEqual(round(quantiles[2],3), round(2.359,3))

    def testRMSEstdCalc(self):

        data = dataUtilities.DataTable(self.regDataPath)
        RFlearner = AZorngRF.RFLearner()
        learners = [RFlearner]
        nFolds = 5
        res = orngTest.crossValidation(learners, data, strat=orange.MakeRandomIndices.StratifiedIfPossible, folds = nFolds) 
        RMSEstd = evalUtilities.getRMSEstd(res, nFolds)[0]
        self.assertEqual(round(RMSEstd,3), round(0.101, 3))



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

