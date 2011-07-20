import unittest
import os

import orange
from AZutilities import getStructuralDesc
import AZOrangeConfig as AZOC

class StructuralDescTest(unittest.TestCase):

    def setUp(self):
        qsarDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/QSAR_10mols.tab")
        self.data = orange.ExampleTable(qsarDataPath)

    def test_bbrc_feature_gen(self):
        """ Test if bbrc freatures are generated correctly """
        import bbrc
        expected_no_features = 6
        features = getStructuralDesc.getBBRCsmartsList(self.data,2)
        self.assertEqual(expected_no_features,len(features))

    def test_bbrc(self):
        """ Test if bbrc is running properly """
        import bbrc
        # the minimum support for BBRC is absolute not relative like for FTM
        result = getStructuralDesc.getBBRCDescResult(self.data, 2)
        expected_atts = 9
        expected_data_length = 10
        self.assertEqual(expected_atts,len(result.domain.attributes))
        self.assertEqual(expected_data_length, len(result[0]))
                
    def test_mkTMPfile(self):
        """Test if temporary SDF files can be written """
        import tempfile
        from subprocess import Popen, PIPE
        from cinfony import rdk
				
        sdf_mols, temp_occ = getStructuralDesc.makeTempFilesFTM(self.data)
        cmd = 'cat ' + sdf_mols.name + ' | grep \'\$\$\$\$\' | wc'
        p = Popen(cmd, shell=True, close_fds=True, stdout=PIPE)
        stdout = p.communicate()
        self.assertEqual(stdout[0].strip(),"10      10      50")
	
		
    def test_ftm(self):
        """Test if ftm is running properly"""
        result = getStructuralDesc.getFTMDescResult(self.data, 0.9)
        expected_atts = 29
        expected_data_length = 30
        self.assertEqual(expected_atts,len(result.domain.attributes))
        self.assertEqual(expected_data_length, len(result[0]))


    def test_smartsRecalc(self):
        """Test structural feature recalculation"""
        result = getStructuralDesc.getFTMDescResult(self.data, 0.9)
        smarts = result.domain.attributes[len(self.data.domain.attributes):]
        
        result2 = getStructuralDesc.getSMARTSrecalcDesc(self.data,smarts)
                
        #expected_atts = 3
        self.assertEqual(len(result),len(result2))
	

		
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(StructuralDescTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
