import unittest
import os

import orange
from AZutilities import AZOrangePredictor 



class AZOrangePredictorTest(unittest.TestCase):

    def setUp(self):
        #self.modelPath = "data/DescModel.model"  # Just RDK descriptors and RDK Fingerprints
        self.modelPath = "data/BBRC.model"
        self.smi = "C123C5C(O)C=CC2C(N(C)CC1)Cc(ccc4O)c3c4O5"  # NEG
        #self.smi = "CCC"  # POS

    def test_AZOrangePredictor(self):
        """ Test of AZOrangePredictor 
        """
        predictor = AZOrangePredictor.AZOrangePredictor(self.modelPath)
        #Needed for classification
        predictor.predictionOutcomes = ["1", "2"]

        #Needed for Regression
        #predictor.significanceThreshold = 0.4

        predictor.getDescriptors(self.smi)
        prediction = predictor.predict()
        significance = predictor.getSDs(self.smi, prediction)
        self.assert_(prediction == "1", "Got:"+str(prediction))   

        # Expecting: {'color': (0, 0, 0), 'imgPath': '', 'signature': '', 'atoms': [], 'non-signature': 'rdk.FP_1510328189', 'molStr': ''}
        self.assert_(significance['non-signature'] ==  'Continuous: \n   [#7&A]-[#6&A]-[#6&A]-[#7&A]\n\n', "Got: "+str(significance))   



		
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(AZOrangePredictorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
