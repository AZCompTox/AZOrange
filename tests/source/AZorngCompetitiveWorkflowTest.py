from AZutilities import dataUtilities
from AZutilities import competitiveWorkflow
import unittest
import os
import time
import AZOrangeConfig as AZOC
import AZorngTestUtil
import pprint

class competitiveWFTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        """Creates the training and testing data set attributes. """
        contTestDataPath  = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Test.tab")
        discTestDataPath  = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Test.tab")
        contTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_Train.tab")
        discTrainDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_Train.tab")

        # Read in the data
        self.Ctrain_data =  dataUtilities.DataTable(contTrainDataPath)
        self.Ctest_data =  dataUtilities.DataTable(contTestDataPath)

        self.Dtrain_data =  dataUtilities.DataTable(discTrainDataPath)
        self.Dtest_data =  dataUtilities.DataTable(discTestDataPath)

        
    def no_testClass_Serial(self):
        """Test classification in serial mode
        """
        res = competitiveWorkflow.competitiveWorkflow(self.Dtrain_data)
        print "Results :  ",res
        self.assert_("statistics" in res)
        self.assert_("model" in res)
        self.assert_("selectedML" in res["statistics"])
        self.assert_(res["model"][res["model"].keys()[0]] is not None)
        self.assertEqual(res["statistics"]["selectedML"]["responseType"], "Classification")
        self.assert_(res["statistics"]["selectedML"]["CA"] > 0 )


    def testReg_Serial(self):
        """Test regression in serial mode
        """ 
        res = competitiveWorkflow.competitiveWorkflow(self.Ctrain_data)
        print "Results :  ",res
        self.assert_("statistics" in res)
        self.assert_("model" in res)
        self.assert_("selectedML" in res["statistics"])
        self.assert_(res["model"][res["model"].keys()[0]] is not None)
        self.assertEqual(res["statistics"]["selectedML"]["responseType"], "Regression")
        self.assert_(res["statistics"]["selectedML"]["Q2"] > 0 )


    def no_testClass_SGE(self):
        """Test classification in serial mode
        """
        res = competitiveWorkflow.competitiveWorkflow(self.Dtrain_data, queueType = "batch.q")
        print "Results :  ",res
        self.assert_("statistics" in res)
        self.assert_("model" in res)
        self.assert_("selectedML" in res["statistics"])
        self.assert_(res["model"][res["model"].keys()[0]] is not None)
        self.assertEqual(res["statistics"]["selectedML"]["responseType"], "Classification")
        self.assert_(res["statistics"]["selectedML"]["CA"] > 0 )


    def no_testReg_SGE(self):
        """Test regression in serial mode
        """
        res = competitiveWorkflow.competitiveWorkflow(self.Ctrain_data, queueType = "batch.q")
        print "Results :  ",res
        self.assert_("statistics" in res)
        self.assert_("model" in res)
        self.assert_("selectedML" in res["statistics"])
        self.assert_(res["model"][res["model"].keys()[0]] is not None)
        self.assertEqual(res["statistics"]["selectedML"]["responseType"], "Regression")
        self.assert_(res["statistics"]["selectedML"]["Q2"] > 0 )



if __name__ == "__main__":
        suite = unittest.TestLoader().loadTestsFromTestCase(competitiveWFTest)
        unittest.TextTestRunner(verbosity=2).run(suite)


