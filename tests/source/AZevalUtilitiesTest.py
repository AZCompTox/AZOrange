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

        trainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/trainDataMDcalc.tab")
        self.trainData = dataUtilities.DataTable(trainDataPath,createNewOn=orange.Variable.MakeStatus.OK)

        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/testDataMDcalc.tab")
        self.testData = dataUtilities.DataTable(testDataPath,createNewOn=orange.Variable.MakeStatus.OK)

        self.regDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")

    def testQuantileCalc(self):

        MD = similarityMetrics.calcMahalanobis(self.trainData, self.testData)
        quantiles = similarityMetrics.calcMahalanobisDistanceQuantiles(MD)

        self.assertEqual(round(quantiles[0],3), round(6.0507258966973936,3))
        self.assertEqual(round(quantiles[1],3), round(9.5167312813901184,3))
        self.assertEqual(round(quantiles[2],3), round(13.764618017361794,3))

    def testRMSEstdCalc(self):

        data = dataUtilities.DataTable(self.regDataPath)
        RFlearner = AZorngRF.RFLearner()
        learners = [RFlearner]
        nFolds = 5
        res = orngTest.crossValidation(learners, data, strat=orange.MakeRandomIndices.StratifiedIfPossible, folds = nFolds) 
        RMSEstd = evalUtilities.getRMSEstd(res, nFolds)[0]
        self.assertEqual(round(RMSEstd,3), round(0.0434402, 3))



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

