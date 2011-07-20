""" Module for similarity boosed QSAR project AZ + TUM
	author: Tobias Girschick; tobias.girschick@in.tum.de
			TUM - I12 (wwwkramer.in.tum.de/girschic)
	dependencies:
"""
import orange
from cinfony import rdk
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs
import AZOrangeConfig as AZOC

def getSimDescriptors(actives, data, methods):
	""" calculates similarity descriptors for a training set (orange object) using the 
		given similarity methods against the given actives
		Possible method strings in methods are the names of the sim_* methods below,
		e.g. rdk_topo_fps for sim_rdk_topo_fps
	"""
	# adjust the header
	atts = []
	for m in methods:
		count = 1
		for a in actives:
			attname = m + '(active_'+ str(count)+ ')'
			atts.append(orange.FloatVariable(attname))
			count += 1	
			
	newdomain = orange.Domain(data.domain.attributes + atts, data.domain.classVar)
	newdata = orange.ExampleTable(newdomain, data)
	
	att_idx = 0	
	# fill up the data	
	for m in methods:
		if m == 'rdk_topo_fps':
			count = 1
			for a in actives:
				attname = m + '(active_'+ str(count)+ ')'
				for j in range(len(newdata)):
					instance = newdata[j]
					tmp = orange.Value(atts[att_idx], orng_sim_rdk_topo_fps(a, instance))
					instance[atts[att_idx]] = tmp
				att_idx += 1	
				
		elif m == 'rdk_MACCS_keys':
			count = 1
			for a in actives:
				attname = m + '(active_'+ str(count)+ ')'
				for j in range(len(newdata)):
					instance = newdata[j]
					tmp = orange.Value(atts[att_idx], orng_sim_rdk_MACCS_keys(a, instance))
					instance[atts[att_idx]] = tmp
				att_idx += 1		

		elif m == 'rdk_morgan_fps':
			count = 1
			for a in actives:
				attname = m + '(active_'+ str(count)+ ')'
				for j in range(len(newdata)):
					instance = newdata[j]
					tmp = orange.Value(atts[att_idx], orng_sim_rdk_morgan_fps(a, instance))
					instance[atts[att_idx]] = tmp		
				att_idx += 1	
				
		elif m == 'rdk_morgan_features_fps':
			count = 1
			for a in actives:
				attname = m + '(active_'+ str(count)+ ')'
				for j in range(len(newdata)):
					instance = newdata[j]
					tmp = orange.Value(atts[att_idx], orng_sim_rdk_morgan_features_fps(a, instance))
					instance[atts[att_idx]] = tmp		
				att_idx += 1	
					
		elif m == 'rdk_atompair_fps':
			count = 1
			for a in actives:
				attname = m + '(active_'+ str(count)+ ')'
				for j in range(len(newdata)):
					instance = newdata[j]
					tmp = orange.Value(atts[att_idx], orng_sim_rdk_atompair_fps(a, instance))
					instance[atts[att_idx]] = tmp		
				att_idx += 1	
							
	return newdata
	
	
	

def orng_sim_rdk_topo_fps(smile_active, train_instance):
	""" calculate the fingerprint similarity using the RDK topological fingerprints
		(The fingerprinting algorithm used is similar to that used in the Daylight fingerprinter)
		input are a smiles string and a orange data instance
		returned is a similaritie value
	"""
	smilesName = getSMILESAttr(train_instance)
	if not smilesName: return None
	smile_train = str(train_instance[smilesName].value)
	
	fp_A = FingerprintMols.FingerprintMol(rdk.Chem.MolFromSmiles(smile_active))
	fp_T = FingerprintMols.FingerprintMol(rdk.Chem.MolFromSmiles(smile_train))
	
	sim = DataStructs.FingerprintSimilarity(fp_A,fp_T)

	return sim
	

def orng_sim_rdk_MACCS_keys(smile_active, train_instance):
	""" calculate the fingerprint similarity using the RDK MACCS keys
		(SMARTS-based implementation of the 166 public MACCS keys)
		input are a smiles string and a orange data instance
		returned is a similaritie value
	"""
	smilesName = getSMILESAttr(train_instance)
	if not smilesName: return None
	smile_train = str(train_instance[smilesName].value)
	
	fp_A = rdk.Chem.MACCSkeys.GenMACCSKeys(rdk.Chem.MolFromSmiles(smile_active))
	fp_T = rdk.Chem.MACCSkeys.GenMACCSKeys(rdk.Chem.MolFromSmiles(smile_train))
	sim = DataStructs.FingerprintSimilarity(fp_A,fp_T)
	
	return sim
			
def orng_sim_rdk_morgan_fps(smile_active, train_instance):
	""" calculate the fingerprint similarity using the RDK morgan fingerprints
		(circular fingerprints, ECFP, connectivity-based invariant)
		input are a smiles string and a orange data instance
		returned is a similaritie value
	"""
	smilesName = getSMILESAttr(train_instance)
	if not smilesName: return None
	smile_train = str(train_instance[smilesName].value)
	
	fp_A = rdk.AllChem.GetMorganFingerprint(rdk.Chem.MolFromSmiles(smile_active),2)
	fp_T = rdk.AllChem.GetMorganFingerprint(rdk.Chem.MolFromSmiles(smile_train),2)
	sim = DataStructs.DiceSimilarity(fp_A,fp_T)
	
	return sim
	
def orng_sim_rdk_morgan_features_fps(smile_active, train_instance):
	""" calculate the fingerprint similarity using the RDK morgan fingerprints
		(circular fingerprints, FCFP, feature-based invariant)
		input are a smiles string and a orange data instance
		returned is a similaritie value
	"""
	smilesName = getSMILESAttr(train_instance)
	if not smilesName: return None
	smile_train = str(train_instance[smilesName].value)
	
	fp_A = rdk.AllChem.GetMorganFingerprint(rdk.Chem.MolFromSmiles(smile_active),2,useFeatures=True)
	fp_T = rdk.AllChem.GetMorganFingerprint(rdk.Chem.MolFromSmiles(smile_train),2,useFeatures=True)
	sim = DataStructs.DiceSimilarity(fp_A,fp_T)
	
	return sim
	
	
def orng_sim_rdk_atompair_fps(smile_active, train_instance):
	""" calculate the fingerprint similarity using the RDK atom pair fingerprints
		input are a smiles string and a orange data instance
		returned is a similaritie value
	"""
	smilesName = getSMILESAttr(train_instance)
	if not smilesName: return None
	smile_train = str(train_instance[smilesName].value)
	
	fp_A = Pairs.GetAtomPairFingerprint(rdk.Chem.MolFromSmiles(smile_active))
	fp_T = Pairs.GetAtomPairFingerprint(rdk.Chem.MolFromSmiles(smile_train))
	sim = DataStructs.DiceSimilarity(fp_A,fp_T)
	
	return sim
	

def get_similarity_matrix(actives, trainset, methods):
	""" calculates similarity descriptors for a training set (list of smiles) using the 
		given similarity methods against the given actives
		Possible method strings in methods are the names of the sim_* methods below,
		e.g. rdk_topo_fps for sim_rdk_topo_fps
	"""
	sim_matrix = []
	for m in methods:
		if m == 'rdk_topo_fps':
			for a in actives:
				sim_matrix.append(sim_rdk_topo_fps(a, trainset))
		elif m == 'rdk_MACCS_keys':
			for a in actives:
				sim_matrix.append(sim_rdk_MACCS_keys(a, trainset))
		elif m == 'rdk_morgan_fps':
			for a in actives:
				sim_matrix.append(sim_rdk_morgan_fps(a, trainset)) 
		elif m == 'rdk_atompair_fps':
			for a in actives:
				sim_matrix.append(sim_rdk_atompair_fps(a, trainset)) 
				
	return sim_matrix



def sim_rdk_topo_fps(smiA, smisT):
	""" calculate the fingerprint similarity using the RDK atompair fingerprints
		input are a smiles string and a list of smiles strings
		returned is a list of similarities
	"""
	fp_A = Pairs.GetAtomPairFingerprint(rdk.Chem.MolFromSmiles(smiA))
	fps_T = [Pairs.GetAtomPairFingerprint(rdk.Chem.MolFromSmiles(y)) for y in smisT]
	
	sim_vector = []
	for t in fps_T:
		sim_vector.append(DataStructs.DiceSimilarity(fp_A,t))

	return sim_vector	
	


def sim_rdk_topo_fps(smiA, smisT):
	""" calculate the fingerprint similarity using the RDK topological fingerprints
		(The fingerprinting algorithm used is similar to that used in the Daylight fingerprinter)
		input are a smiles string and a list of smiles strings
		returned is a list of similarities
	"""
	fp_A = FingerprintMols.FingerprintMol(rdk.Chem.MolFromSmiles(smiA))
	fps_T = [FingerprintMols.FingerprintMol(rdk.Chem.MolFromSmiles(y)) for y in smisT]
	
	sim_vector = []
	for t in fps_T:
		sim_vector.append(DataStructs.FingerprintSimilarity(fp_A,t))

	return sim_vector		
		
		
		
def sim_rdk_MACCS_keys(smiA, smisT):
	""" calculate the fingerprint similarity using the RDK MACCS keys
		(SMARTS-based implementation of the 166 public MACCS keys)
		input are a smiles string and a list of smiles strings
		returned is a list of similarities
	"""
	fp_A = rdk.Chem.MACCSkeys.GenMACCSKeys(rdk.Chem.MolFromSmiles(smiA))
	fps_T = [rdk.Chem.MACCSkeys.GenMACCSKeys(rdk.Chem.MolFromSmiles(y)) for y in smisT]
	
	sim_vector = []
	for t in fps_T:
		sim_vector.append(DataStructs.FingerprintSimilarity(fp_A,t))
	
	return sim_vector	
			
def sim_rdk_morgan_fps(smiA, smisT):
	""" calculate the fingerprint similarity using the RDK morgan fingerprints
		(circular fingerprints)
		input are a smiles string and a list of smiles strings
		returned is a list of similarities
	"""
	fp_A = rdk.AllChem.GetMorganFingerprint(rdk.Chem.MolFromSmiles(smiA),2)
	fps_T = [rdk.AllChem.GetMorganFingerprint(rdk.Chem.MolFromSmiles(y),2) for y in smisT]
	
	sim_vector = []
	for t in fps_T:
		sim_vector.append(DataStructs.DiceSimilarity(fp_A,t))
	
	return sim_vector


def getSMILESAttr(data):
    # Check that the data contains a SMILES attribute
    smilesName = None
    for attr in [a.name for a in  data.domain] + [a.name for a in data.domain.getmetas().values()]:
        if attr in AZOC.SMILESNAMES:
            smilesName = attr
    if not smilesName:
        print "Warning: The data set does not contain any known smiles attribute!"
        print "No Cinfony descriptors added!"
        return None
    else:       
        return smilesName
