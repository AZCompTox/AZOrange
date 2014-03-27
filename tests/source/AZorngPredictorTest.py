import unittest
import os

import orange
from AZutilities import AZOrangePredictor 
from AZutilities import dataUtilities



class AZOrangePredictorTest(unittest.TestCase):


    def test_AZOrangePredictor(self):
        """ Test of AZOrangePredictor 
        """
        self.modelPath = "data/QTcB_SVM_Sign_Model"  # Signatures hight 3
        train = dataUtilities.DataTable("data/QTcB_sign.txt")
        self.ex = dataUtilities.DataTable([train[0]])  # Test with precalc signatures 

        predictor = AZOrangePredictor.AZOrangePredictor(self.modelPath)

        #Needed for Regression
        predictor.significanceThreshold = 0.0

        predictor.setDescriptors(self.ex)
        prediction = predictor.predict()
        significance = predictor.getSDs(None, prediction, exWithDesc = self.ex)

        print prediction
        print significance
        self.assertEqual(round(prediction, 5), 0.01787)
        self.assertEqual(significance, None)    # All gradients are 0


    def test_signHeight1(self):

        self.modelPath = "data/QTcB_SVM_Sign1_Model"  # Signatures hight 1
        train = dataUtilities.DataTable("data/QTcB_sign1.txt")
        self.smiles = "OC(=O)C1=CN(C2CC2)c3cc(N4CCNCC4)c(F)cc3C1=O"  # Test with smiles

        predictor = AZOrangePredictor.AZOrangePredictor(self.modelPath)

        #Needed for Regression
        predictor.significanceThreshold = 0.0

        predictor.getDescriptors(self.smiles)
        prediction = predictor.predict()
        significance = predictor.getSDs(self.smiles, prediction, topN = 10)
        
        print prediction
        print significance
        self.assert_(round(prediction,5) == 0.02973, "Got: "+str(prediction))   
        self.assert_(str(significance['signature']) =="['[Cac]', '[Car]([Car][Car][O2])', '[Cac]([Car][O.co2][O.co2])', '[C3]([C3][C3][Nar])', '[Car]([Car][Nar][O2])', '[O.co2]', '[C3]([C3][Nar])', '[C3]([C3][C3][N3])', '[Nar]([C3][Car][Car])', '[C3]([Car][N3])']"  , "Got: "+str(significance['signature']))   


		
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(AZOrangePredictorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
