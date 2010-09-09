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
        self.assert_(len(allDescs) > 230)
        descs = ['webel.AtomCountDescriptor', 'webel.AutocorrelationDescriptorCharge', 'webel.AutocorrelationDescriptorMass',"obabel.TPSA","obabel.HBA2" , "obabel.MW","rdk.TPSA","rdk.Chi0n","rdk.MolWt"]
        resD = getCinfonyDesc.getCinfonyDescResults(self.smiData,descs)
        self.assertEqual(len(resD),len(self.smiData))
        self.assertEqual(str(resD.domain),"[SMILES, TPSA, HBA2, MW, TPSA_1, Chi0n, MolWt, AtomCountDescriptor_nAtom, AutocorrelationDescriptorMass_ATSm5, AutocorrelationDescriptorMass_ATSm4, AutocorrelationDescriptorMass_ATSm3, AutocorrelationDescriptorMass_ATSm2, AutocorrelationDescriptorMass_ATSm1, AutocorrelationDescriptorCharge_ATSc5, AutocorrelationDescriptorCharge_ATSc3, AutocorrelationDescriptorCharge_ATSc1, AutocorrelationDescriptorCharge_ATSc4]")
        


    def test_openbabel(self):
        from cinfony import obabel
        mol = obabel.readstring("smi", "CCC")
        desc = mol.calcdesc()
        for d in desc:
            if desc[d] != desc[d]:
                desc[d] = '?'
        expectedDesc = {'TPSA': 0.0, 'smarts': '?', 'logP': 1.4163000000000001, 'title': '?', 'HBD': 0.0, 's': '?', 'nF': 0.0, 'InChI': '?', 'MW': 44.095619999999997, 'L5': '?', 'MR': 16.535, 'HBA1': 0.0, 'spinMult': 1.0, 'HBA2': 0.0, 'nHal': 0.0}
        self.assertEqual(desc,expectedDesc)

    def test_RDKit(self):
        from cinfony import rdk
        mol = rdk.readstring("smi", "CCC")
        desc = mol.calcdesc()

        expectedDesc = {'BertzCT': 0.0, 'fr_C_O_noCOO': 0, 'Chi4v': 0.0, 'fr_Ar_COO': 0, 'Chi4n': 0.0, 'SMR_VSA4': 0.0, 'fr_urea': 0, 'fr_para_hydroxylation': 0, 'fr_barbitur': 0, 'fr_halogen': 0, 'fr_dihydropyridine': 0, 'fr_priamide': 0, 'fr_Al_COO': 0, 'fr_guanido': 0, 'fr_furan': 0, 'fr_morpholine': 0, 'fr_term_acetylene': 0, 'fr_amidine': 0, 'fr_benzodiazepine': 0, 'SlogP_VSA1': 0.0, 'MolWt': 44.096999999999994, 'fr_hdrzine': 0, 'fr_Ar_NH': 0, 'fr_quatN': 0, 'fr_benzene': 0, 'fr_phos_acid': 0, 'fr_sulfone': 0, 'VSA_EState10': 0.0, 'fr_aniline': 0, 'fr_N_O': 0, 'fr_sulfonamd': 0, 'fr_thiazole': 0, 'TPSA': 0.0, 'fr_piperzine': 0, 'SMR_VSA10': 0.0, 'PEOE_VSA13': 0.0, 'PEOE_VSA12': 0.0, 'PEOE_VSA11': 0.0, 'PEOE_VSA10': 0.0, 'BalabanJ': 1.6329931618554523, 'fr_lactone': 0, 'Chi3v': 0.0, 'Chi2n': 0.70710678118654746, 'EState_VSA10': 0.0, 'EState_VSA11': 0.0, 'HeavyAtomMolWt': 36.033000000000001, 'Chi0': 2.7071067811865475, 'Chi1': 1.4142135623730951, 'MolLogP': 1.4163000000000001, 'fr_nitro': 0, 'fr_Al_OH': 0, 'fr_azo': 0, 'fr_C_O': 0, 'fr_ether': 0, 'fr_phenol_noOrthoHbond': 0, 'RingCount': 0, 'fr_alkyl_halide': 0, 'NumValenceElectrons': 20.0, 'fr_aryl_methyl': 0, 'HallKierAlpha': 0.0, 'fr_C_S': 0, 'fr_thiocyan': 0, 'fr_NH0': 0, 'VSA_EState4': 0.0, 'VSA_EState5': 0.0, 'VSA_EState6': 0.0, 'VSA_EState7': 0.0, 'NumHDonors': 0, 'VSA_EState2': 0.0, 'VSA_EState3': 0.0, 'fr_HOCCN': 0, 'fr_NH2': 0, 'VSA_EState8': 0.0, 'VSA_EState9': 5.5, 'SlogP_VSA10': 0.0, 'SlogP_VSA11': 0.0, 'fr_COO': 0, 'NHOHCount': 0, 'fr_unbrch_alkane': 0, 'fr_methoxy': 0, 'fr_amide': 0, 'SlogP_VSA8': 0.0, 'SlogP_VSA9': 0.0, 'SlogP_VSA4': 0.0, 'SlogP_VSA5': 20.268296022307258, 'SlogP_VSA6': 0.0, 'SlogP_VSA7': 0.0, 'fr_Imine': 0, 'SlogP_VSA2': 0.0, 'SlogP_VSA3': 0.0, 'fr_phos_ester': 0, 'SMR_VSA3': 0.0, 'NumHeteroatoms': 0, 'fr_NH1': 0, 'fr_ketone_Topliss': 0, 'fr_SH': 0, 'LabuteASA': 21.469135256622767, 'fr_thiophene': 0, 'Chi3n': 0.0, 'fr_imidazole': 0, 'fr_nitrile': 0, 'SMR_VSA2': 0.0, 'SMR_VSA1': 0.0, 'fr_nitro_arom': 0, 'SMR_VSA6': 0.0, 'EState_VSA8': 13.847474399381248, 'EState_VSA9': 0.0, 'EState_VSA6': 0.0, 'SMR_VSA7': 0.0, 'SMR_VSA9': 0.0, 'SMR_VSA8': 0.0, 'EState_VSA2': 0.0, 'fr_Ndealkylation2': 0, 'fr_Ndealkylation1': 0, 'EState_VSA1': 0.0, 'PEOE_VSA14': 0.0, 'Kappa3': 0, 'fr_nitroso': 0, 'fr_diazo': 0, 'Kappa2': 2.0, 'fr_Ar_N': 0, 'fr_Nhpyrrole': 0, 'EState_VSA7': 0.0, 'MolMR': 15.965, 'EState_VSA4': 0.0, 'fr_COO2': 0, 'fr_prisulfonamd': 0, 'fr_oxime': 0, 'EState_VSA5': 6.4208216229260087, 'fr_isocyan': 0, 'EState_VSA3': 0.0, 'Chi2v': 0.70710678118654746, 'HeavyAtomCount': 3, 'fr_aldehyde': 0, 'SMR_VSA5': 20.268296022307258, 'NumHAcceptors': 0, 'fr_lactam': 0, 'fr_allylic_oxid': 0, 'SlogP_VSA12': 0.0, 'fr_oxazole': 0, 'fr_piperdine': 0, 'fr_Ar_OH': 0, 'fr_sulfide': 0, 'fr_alkyl_carbamate': 0, 'NOCount': 0, 'Chi1n': 1.4142135623730951, 'PEOE_VSA8': 0.0, 'PEOE_VSA7': 0.0, 'PEOE_VSA6': 20.268296022307258, 'PEOE_VSA5': 0.0, 'PEOE_VSA4': 0.0, 'PEOE_VSA3': 0.0, 'PEOE_VSA2': 0.0, 'PEOE_VSA1': 0.0, 'fr_imide': 0, 'Chi1v': 1.4142135623730951, 'fr_Al_OH_noTert': 0, 'fr_epoxide': 0, 'fr_hdrzone': 0, 'fr_isothiocyan': 0, 'VSA_EState1': 0.0, 'fr_bicyclic': 0, 'Kappa1': 3.0, 'Ipc': 2.7548875021634687, 'Chi0n': 2.7071067811865475, 'fr_phenol': 0, 'fr_ester': 0, 'PEOE_VSA9': 0.0, 'fr_azide': 0, 'fr_pyridine': 0, 'fr_tetrazole': 0, 'fr_ketone': 0, 'fr_nitro_arom_nonortho': 0, 'Chi0v': 2.7071067811865475, 'fr_ArN': 0, 'NumRotatableBonds': 0}  
        self.assertEqual(desc,expectedDesc)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

