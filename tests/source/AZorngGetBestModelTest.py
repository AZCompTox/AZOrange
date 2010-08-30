import unittest
import os
import time
import sys

import orange
from AZutilities import getBestModel
from AZutilities import miscUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil


class getBestModelTest(AZorngTestUtil.AZorngTestUtil):

    def setUp(self):
        
        #These 2 datasets are equal apart from the meta atribute
        self.trainPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyTrain.tab")
        self.testPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummyTest.tab")

        self.trainPath2 = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/dummy.tab")


    def testTopRankedNoGrid(self):
        """
        Test the TopRanked method in getBestModel with a test set and without grid computing. 
        """
        # Fix input arg
        resultsDir = miscUtilities.createScratchDir(desc="GetBestModelTest") 
        descList = [5, 10]
        grid = False
        batchQueue = False
        optParam = True

        # Method to test
        getBestModel.getBestModelTopRank(self.trainPath, self.testPath, resultsDir, descList, grid, optParam, batchQueue)

        # Assert the existens of a results file
        resultsFile = resultsDir+"/batchResults.tex"
        self.assert_(os.path.exists(resultsFile), "No results file created with getBestModelTopRank")
            
        resultsFile = resultsDir+"/batchResults.pdf"
        self.assert_(os.path.exists(resultsFile), "No pdf file created with getBestModelTopRank")
        miscUtilities.removeDir(resultsDir)

    def testDescSetNoGrid(self):
        """
        Test the descSet method in getBestModel without a test set and without grid computing. 
        """
        # Fix input arg
        resultsDir = miscUtilities.createScratchDir(desc="GetBestModelTest")  
        descList = ["AZ_descriptors"]#, "SELMA"]
        grid = False
        batchQueue = False
        optParam = False

        # Method to test
        getBestModel.getBestModelDescSet(self.trainPath2, "noTest", resultsDir, descList, grid, optParam, batchQueue)

        # Assert the existens of a results file
        resultsFile = resultsDir+"/batchResults.tex"
        self.assert_(os.path.exists(resultsFile), "No results file created with getBestModelDescSet")
            
        resultsFile = resultsDir+"/batchResults.pdf"
        self.assert_(os.path.exists(resultsFile), "No pdf file created with getBestModelDescSet")
        
        miscUtilities.removeDir(resultsDir)



if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(getBestModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)



