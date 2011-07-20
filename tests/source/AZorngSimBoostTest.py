import unittest
import os

import orange
from AZutilities import SimBoostedQSAR
import AZOrangeConfig as AZOC

class SimBoostedQsarTest(unittest.TestCase):

    def setUp(self):
        qsarDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/QSAR_10mols.tab")
        self.data = orange.ExampleTable(qsarDataPath)

    
    def test_simBoost(self):
        smilesActive = ['C1CCN(CC1)S(=O)(=O)C2=CC=C(C=C2)C3=NC(=NN3)SCC(=O)NC4=NC=CS4',
            'CC1=CC(=NC2=NC(=NN12)C(=O)N3CCOCC3)C']
        methods = ['rdk_topo_fps', 'rdk_MACCS_keys', 'rdk_morgan_fps', 'rdk_atompair_fps', 'rdk_morgan_features_fps']
        
        c = 0
        for me in methods:
            result = SimBoostedQSAR.getSimDescriptors(smilesActive, self.data, [me])
            expected_atts = 5  # 3 existing + 2*1similarities
            self.assertEqual(expected_atts,len(result.domain.attributes))
            
            expected_values = [0.286,0.530,0.270,0.244,0.369]
            
            for i in range(1):
                self.assertEqual(expected_values[c], round(float(result[i][4]), 3))
            print
            c += 1


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(SimBoostedQsarTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
