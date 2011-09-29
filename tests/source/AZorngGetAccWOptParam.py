import logging
import os
import unittest

from AZutilities import  getUnbiasedAccuracy
from trainingMethods import AZorngRF
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil

class GetAccWOptParam(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        logging.basicConfig(level=logging.ERROR)
        self.log = logging.getLogger("GetAccWOptParam")

        irisPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris.tab")
        iris2Path = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris2.tab")
        irisContPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/irisCont.tab")
        # Read in the data
        self.irisData = dataUtilities.DataTable(irisPath)
        self.iris2Data = dataUtilities.DataTable(iris2Path)
        self.irisContData = dataUtilities.DataTable(irisContPath)

    def test_Classification3CVal(self):
        """Testing Classification problem with 3 class values"""
        learner = AZorngRF.RFLearner()
        paramList = ["nActVars"]
        evaluator = getUnbiasedAccuracy.UnbiasedAccuracyGetter(data = self.irisData, learner = learner, paramList = paramList, nExtFolds = 3, nInnerFolds = 3)
        res = evaluator.getAcc()
        self.assert_(abs(res["CA"]-0.96666666666666667) < 0.01)
        expected =  [[50.0, 0.0, 0.0], [0.0, 48.0, 2.0], [0.0, 3.0, 47.0]]
        for i,c in enumerate(res["CM"]):
            for j,l in enumerate(c):
                self.assert_(abs(l-expected[i][j]) < 3)


    def test_Classification2CVal(self):
        """Testing Classification problem with 2 class values"""
        learner = AZorngRF.RFLearner()
        paramList = ["nActVars"]
        evaluator = getUnbiasedAccuracy.UnbiasedAccuracyGetter(data = self.iris2Data, learner = learner, paramList = paramList, nExtFolds = 3, nInnerFolds = 3)
        res = evaluator.getAcc()
        self.assertEqual(round(res["CA"],5),round(0.96666666666666667,5))
        self.assertEqual(res["CM"],  [[98.0, 2.0], [3.0, 47.0]])


    def test_RegressionRMSE(self):
        """Testing Regression RMSE problem"""
        learner = AZorngRF.RFLearner()
        paramList = ["nActVars"]
        evaluator = getUnbiasedAccuracy.UnbiasedAccuracyGetter(data = self.irisContData, learner = learner, paramList = paramList ,nExtFolds = 3, nInnerFolds = 3)
        res = evaluator.getAcc()
        expectedRes = [0.27741430697239661, 0.27945999999999999, 0.276116805384, 0.277488734272, 0.276164200118]  # [InHouse, Ubuntu, Ubuntu 64 bits]
        self.log.info("")
        self.log.info("res.keys()" + str(res.keys()))
        expected = [round(x,5) for x in expectedRes]
        actual = round(res["RMSE"],5)
        self.log.info("expected=" + str(expected))
        self.log.info("actual  =" + str(actual))
        self.assert_( actual in expected,"Got: "+str(res["RMSE"]))


    def test_RegressionR2(self):
        """Testing Regression R2 problem"""
        learner = AZorngRF.RFLearner()
        paramList = ["nActVars"]
        evaluator = getUnbiasedAccuracy.UnbiasedAccuracyGetter(data = self.irisContData, learner = learner, paramList = paramList ,nExtFolds = 3, nInnerFolds = 3)
        res = evaluator.getAcc()
        expectedRes = [
                       0.97464488216654444, 0.97420405774, 0.974887867109, 0.97510044677,        # [InHouse, Ubuntu, Ubuntu 64 bits]
                       0.97533758502003
                      ] 
        self.log.info("")
        self.log.info("res.keys()" + str(res.keys()))
        if "R2" in res.keys():
            acctual = res["R2"]
        elif "Q2" in res.keys():
            acctual = res["Q2"]
        else:
            acctual = None
        self.assertNotEqual(None, acctual, "Could not find either 'R2' nor 'Q2' in the result.")            
        self.assertRoundedToExpectedArray(acctual, expectedRes, 5)



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(GetAccWOptParam)
    unittest.TextTestRunner(verbosity=2).run(suite)


