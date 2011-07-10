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

        trainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_FullNumeric_Train.tab")
        self.trainData = dataUtilities.DataTable(trainDataPath,createNewOn=orange.Variable.MakeStatus.OK)

        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_FullNumeric_Test.tab")
        self.testData = dataUtilities.DataTable(testDataPath,createNewOn=orange.Variable.MakeStatus.OK)

        self.regDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")


    def testMahalanobis(self):
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Mahalanobis/trainData.tab"))
        testData =dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Mahalanobis/testData.tab"))
        invCovMatFile =os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Mahalanobis/SqrtInvCovMatrix.npy")
        #invCovMatFile =os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Mahalanobis/invCovMatrix.npy")
        centerFile =os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Mahalanobis/center.npy")
        MahalanobisData =os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Mahalanobis/mahalanobisDataTable.npy")
        #MahalanobisData =os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Mahalanobis/dataTable.npy")
        domain = data.domain

        
        MD1 = similarityMetrics.calcMahalanobis(data, testData)
        MD2 = similarityMetrics.calcMahalanobis(data, testData, invCovMatFile, centerFile, MahalanobisData, domain)
        expected = [{'_train_id_near1': 'XX', '_train_dist_near1': 8.5529724254000528e-15, '_train_dist_near2': 3.2469494012727824, '_train_id_near2': 'XX', '_train_dist_near3': 3.2731657154209062, '_train_av3nearest': 2.1733717055645658, '_train_SMI_near3': 'XXX', '_train_id_near3': 'XX', '_train_SMI_near1': 'XXX', '_train_SMI_near2': 'XXX', '_MD': 17.204748057165233}, {'_train_id_near1': 'XX', '_train_dist_near1': 5.115384476164845, '_train_dist_near2': 5.5191228149142137, '_train_id_near2': 'XX', '_train_dist_near3': 5.5944864806135319, '_train_av3nearest': 5.4096645905641978, '_train_SMI_near3': 'XXX', '_train_id_near3': 'XX', '_train_SMI_near1': 'XXX', '_train_SMI_near2': 'XXX', '_MD': 32.561505764913186}, {'_train_id_near1': 'XX', '_train_dist_near1': 4.2608827352786338e-15, '_train_dist_near2': 2.7885310210397281, '_train_id_near2': 'XX', '_train_dist_near3': 2.7885310210397281, '_train_av3nearest': 1.8590206806931533, '_train_SMI_near3': 'XXX', '_train_id_near3': 'XX', '_train_SMI_near1': 'XXX', '_train_SMI_near2': 'XXX', '_MD': 14.565865561707858}, {'_train_id_near1': 'XX', '_train_dist_near1': 5.1696950502912085e-15, '_train_dist_near2': 3.2469494012428579, '_train_id_near2': 'XX', '_train_dist_near3': 3.273165715391364, '_train_av3nearest': 2.1733717055447421, '_train_SMI_near3': 'XXX', '_train_id_near3': 'XX', '_train_SMI_near1': 'XXX', '_train_SMI_near2': 'XXX', '_MD': 17.204748056995847}, {'_train_id_near1': 'XX', '_train_dist_near1': 2.3478624106042396, '_train_dist_near2': 3.9525803005509639, '_train_id_near2': 'XX', '_train_dist_near3': 3.9997921263766512, '_train_av3nearest': 3.4334116125106182, '_train_SMI_near3': 'XXX', '_train_id_near3': 'XX', '_train_SMI_near1': 'XXX', '_train_SMI_near2': 'XXX', '_MD': 22.070533877619617}, {'_train_id_near1': 'XX', '_train_dist_near1': 3.2549006987559131e-15, '_train_dist_near2': 3.2469494012420288, '_train_id_near2': 'XX', '_train_dist_near3': 3.2731657153902853, '_train_av3nearest': 2.1733717055441057, '_train_SMI_near3': 'XXX', '_train_id_near3': 'XX', '_train_SMI_near1': 'XXX', '_train_SMI_near2': 'XXX', '_MD': 17.204748057000234}, {'_train_id_near1': 'XX', '_train_dist_near1': 1.887904134152689, '_train_dist_near2': 3.0871126343782667, '_train_id_near2': 'XX', '_train_dist_near3': 3.0950749249075327, '_train_av3nearest': 2.6900305644794962, '_train_SMI_near3': 'XXX', '_train_id_near3': 'XX', '_train_SMI_near1': 'XXX', '_train_SMI_near2': 'XXX', '_MD': 16.634324652887784}]


        dists = ["_MD",'_train_dist_near1','_train_dist_near2','_train_dist_near3','_train_av3nearest']
        for idx,x in enumerate(expected):
            for d in dists:
                self.assert_(abs(MD2[idx][d]-x[d]) < 0.00001, "MD2: idx "+str(idx) +"    diff = "+str(MD2[idx][d]-x[d]))

        for idx,x in enumerate(expected):
            for d in dists:
                self.assert_(abs(MD1[idx][d]-x[d]) < 0.00001, "MD1: idx "+str(idx) +"    diff = "+str(MD1[idx][d]-x[d]))


    def testQuantileCalc(self):

        MD = similarityMetrics.calcMahalanobis(self.trainData, self.testData)
        quantiles = similarityMetrics.calcMahalanobisDistanceQuantiles(MD)
        self.assertEqual(round(quantiles[0],3), round(1.1802994935528328,3))
        self.assertEqual(round(quantiles[1],3), round(1.6345715553057074,3))
        self.assertEqual(round(quantiles[2],3), round(2.316536184143136,3))

    def testgetRsqrt(self):
        data = dataUtilities.DataTable(self.regDataPath)
        RFlearner = AZorngRF.RFLearner()
        trainData = data[0:int(len(data)/2)]
        testData = data[int(len(data)/2)+1:]
        classifier = RFlearner(data)
        Rsqrt = evalUtilities.getRsqrt(testData,classifier)
        self.assert_(Rsqrt-0.684011336894 < 0.05, "Got:"+str(Rsqrt))

    def testgetRMSE(self):
        data = dataUtilities.DataTable(self.regDataPath)
        RFlearner = AZorngRF.RFLearner()
        trainData = data[0:int(len(data)/2)]
        testData = data[int(len(data)/2)+1:]
        classifier = RFlearner(data)
        RMSE = evalUtilities.getRMSE(testData,classifier)
        self.assert_(RMSE-2.07396535555 < 0.05, "Got:"+str(RMSE))
        

    def testRMSEstdCalc(self):

        data = dataUtilities.DataTable(self.regDataPath)
        RFlearner = AZorngRF.RFLearner()
        learners = [RFlearner]
        nFolds = 5
        res = orngTest.crossValidation(learners, data, strat=orange.MakeRandomIndices.StratifiedIfPossible, folds = nFolds) 
        RMSEstd = evalUtilities.getRMSEstd(res, nFolds)[0]
        self.assertEqual(round(RMSEstd,3), round(0.141, 3))



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

