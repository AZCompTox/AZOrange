import unittest
import os

import orange
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC
from trainingMethods import AZorngRF
import Orange
import os
import sys



def printResults(res, estIdx, methodIdx, methodDict):

    learnerIdx = 0
    fileName = os.path.join(resultsDir, "results"+str(Orange.evaluation.reliability.METHOD_NAME[methodIdx])+".txt")
    fid = open(fileName, "w")
    fid.write("Reliability\tPrediction Error\tAbsolute Prediction Error\tfoldIdx\n")
    for ex in res.results:
        predError = ex.actualClass - ex.classes[learnerIdx]
        fid.write(str(ex.probabilities[learnerIdx].reliability_estimate[estIdx].estimate)+"\t"+str(predError)+"\t"+str(abs(predError))+"\t"+str(ex.iterationNumber)+"\n")
    fid.close()



def getReliability(data, learner):
    """
    Evaluate all reliability metrics with a model type and a data set.
    """
    BGV =  Orange.evaluation.reliability.BaggingVariance(m=50)
    CNK = Orange.evaluation.reliability.CNeighbours(k=5)
    # Define reliability methods
    estimators = [Orange.evaluation.reliability.Mahalanobis(k=3),                
                  Orange.evaluation.reliability.LocalCrossValidation(k = 10),   
                  #BGV,                                                          # Included in BVCK
                  #CNK,                    
                  Orange.evaluation.reliability.BaggingVarianceCNeighbours(BGV, CNK)    
                  ]
    
    # Create reliability objects which can be used as learners
    reliability = Orange.evaluation.reliability.Learner(learner, estimators = estimators)

    # CV assessment of the reliability objects
    res = Orange.evaluation.testing.cross_validation([reliability], data, folds = 2)

    results = "RMSE\t"+str(Orange.evaluation.scoring.RMSE(res))+"\n"
    results += "R2\t"+str(Orange.evaluation.scoring.R2(res))+"\n"

    # Print calculate the correlation between reliability and prediction error
    reliability_res = Orange.evaluation.reliability.get_pearson_r(res)

    results += "Estimate               r       p\n"
    estIdx = 0
    for estimate in reliability_res:
        #print estimate
        results += "%-20s %7.3f %7.3f \n" % (Orange.evaluation.reliability.METHOD_NAME[estimate[3]], \
                                     estimate[0], estimate[1])
        #printResults(res, estIdx, estimate[3], resultsDir)
        estIdx = estIdx + 1
    return results



class ReliabilityTest(unittest.TestCase):

    def setUp(self):
        ACEdataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/ACE_RDK.tab")
        ACEresPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/reliabilityACETestRes.txt")
        self.data = orange.ExampleTable(ACEdataPath)
        fh = open(ACEresPath,"r")
        self.resultLines = fh.read()
        fh.close()

    def test_reliability(self):
        """ Test of reliability methods """
        learner = AZorngRF.RFLearner()
        res = getReliability(self.data, learner)
        self.assertEqual(res, self.resultLines)
                


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(ReliabilityTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

