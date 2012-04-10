import unittest
import os
import time

from AZutilities import dataUtilities
from AZutilities import getCinfonyDesc
import AZOrangeConfig as AZOC

class evalUtilitiesTest(unittest.TestCase):

    def setUp(self):
        
        smiDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mol.smi")
        self.smiData = dataUtilities.loadSMI(smiDataPath)

    def test_getCinfonyDesc(self):
        allDescs = getCinfonyDesc.getAvailableDescs()
        print "\nN Retrieved desc: ",len(allDescs)
        self.assert_(len(allDescs) > 200)  # Actualy 220
        descs = ['cdk.BCUT', 'cdk.weight','webel.AtomCountDescriptor', 'webel.AutocorrelationDescriptorCharge', 'webel.AutocorrelationDescriptorMass',"obabel.TPSA","obabel.HBA2" , "obabel.MW","rdk.TPSA","rdk.Chi0n","rdk.MolWt"]
        resD = getCinfonyDesc.getCinfonyDescResults(self.smiData,descs)
        self.assertEqual(len(resD),len(self.smiData))

        
    def test_webel(self):
        from cinfony import webel
        mol = webel.readstring("smi", "CCC")
        desc = mol.calcdesc(webel.getdescs()[5:7])
        for d in desc:
            if desc[d] != desc[d]:
                desc[d] = '?'
        expectedDesc = {'AtomCountDescriptor_nAtom': 11.0, 'AutocorrelationDescriptorCharge_ATSc5': 0.0, 'AutocorrelationDescriptorCharge_ATSc4': 0.0, 'AutocorrelationDescriptorCharge_ATSc3': 0.0, 'AutocorrelationDescriptorCharge_ATSc2': 0.0, 'AutocorrelationDescriptorCharge_ATSc1': 0.0}
        self.assertEqual(desc,expectedDesc)



    def test_cdk(self):
        from cinfony import cdk
        mol = cdk.readstring("smi", "CCC")
        desc = mol.calcdesc(cdk.descs[1:3])
        for d in desc:
            if desc[d] != desc[d]:
                desc[d] = '?'
        expectedDesc = {'ip': 0.0, 'weight': 44.062600320000001}
        self.assertEqual(desc,expectedDesc)



    def test_openbabel(self):
        from cinfony import obabel
        mol = obabel.readstring("smi", "CCC")
        desc = mol.calcdesc(obabel.descs[0:3])
        for d in desc:
            if desc[d] != desc[d]:
                desc[d] = '?'
        #expectedDesc = {'HBA1': 0.0, 'HBD': 0.0, 'HBA2': 0.0}
        expectedDesc = {'abonds': 0.0, 'atoms': 3.0, 'bonds': 2.0}
        self.assertEqual(desc,expectedDesc)

    def test_RDKit(self):
        from cinfony import rdk
        mol = rdk.readstring("smi", "CCC")
        desc = mol.calcdesc(rdk.descs[0:3])
        for d in desc:
            if desc[d] != desc[d]:
                desc[d] = '?'
        expectedDesc = {'fr_Ar_COO': 0, 'Chi4v': 0.0, 'fr_C_O_noCOO': 0}
        self.assertEqual(desc,expectedDesc)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

