import unittest
import os
import time

from AZutilities import dataUtilities
from AZutilities import SimBoostedQSAR
import AZOrangeConfig as AZOC

class evalUtilitiesTest(unittest.TestCase):

    def setUp(self):
        
        smiDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mol.smi")
        activesDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/ace_fmc.smi")
        self.smiData = dataUtilities.loadSMI(smiDataPath)
        self.activesData = dataUtilities.loadSMI(activesDataPath)
        self.methods = SimBoostedQSAR.methods.values()[0:1]

    def test_SimDesc_standard_smiles(self):
        actives = []
        for ex in self.smiData:
            actives.append(str(ex["SMILES"].value))
        newData = SimBoostedQSAR.getSimDescriptors(actives, self.smiData, self.methods)
        self.assertEqual(len(newData),len(self.smiData))
        self.assertEqual(newData.has_missing_values(),0)

    def test_SimDesc_non_standard_smiles(self):
        actives = []
        for ex in self.activesData:
        #for ex in self.smiData:
            actives.append(str(ex["SMILES"].value))
        newData = SimBoostedQSAR.getSimDescriptors(actives, self.smiData, self.methods)
        self.assertEqual(len(newData),len(self.smiData))
        self.assertEqual(newData.has_missing_values(),0)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

