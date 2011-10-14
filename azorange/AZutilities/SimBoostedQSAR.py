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

methods = { "Rdk Topo fps"               :'rdk_topo_fps',
            "Rdk MACCS keys"             :'rdk_MACCS_keys',
            "Rdk Morgam fps"             :'rdk_morgan_fps',
            "Rdk Morgan Features fps"    :'rdk_morgan_features_fps',
            "Rdk Atom Pair fps"          :'rdk_atompair_fps'
            #"AZO-pharmacophore fps"    :'azo_pharmacophore_fps'
} 

def getSimDescriptors(actives, data, methods, active_ids = None, pharmacophore_file = None, callBack = None):
        """ calculates similarity descriptors for a training set (orange object) using the 
                given similarity methods against the given actives
                Possible method strings in methods are the names of the sim_* methods below,
                e.g. rdk_topo_fps for sim_rdk_topo_fps
            callBack function, if defined, will be called on each step sending the pergentage done (0-100): 
                   e.g. callBack(25)
                the callBack function shall return True of False which will indicate to this method if the process it to be continued or Not.
                   e.g. if callBack(25) == False it indicates the caller want's to stop the process of calculating descriptors                 
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
        # if callBack is defined, it will be called with the percentage done, i.e. 0-100
        if active_ids:
            nTotalSteps = len(newdata) * ( (len(methods)-1) * len(actives) + len(active_ids) )
        else:
            nTotalSteps = len(methods) * len(actives) * len(newdata)
        stepsDone   = 0  
        
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
                                        if callBack: 
                                            stepsDone += 1
                                            if not callBack((100*stepsDone)/nTotalSteps): return None
                                att_idx += 1        
                                
                elif m == 'rdk_MACCS_keys':
                        count = 1
                        for a in actives:
                                attname = m + '(active_'+ str(count)+ ')'
                                for j in range(len(newdata)):
                                        instance = newdata[j]
                                        tmp = orange.Value(atts[att_idx], orng_sim_rdk_MACCS_keys(a, instance))
                                        instance[atts[att_idx]] = tmp
                                        if callBack: 
                                            stepsDone += 1
                                            if not callBack((100*stepsDone)/nTotalSteps): return None

                                att_idx += 1                

                elif m == 'rdk_morgan_fps':
                        count = 1
                        for a in actives:
                                attname = m + '(active_'+ str(count)+ ')'
                                for j in range(len(newdata)):
                                        instance = newdata[j]
                                        tmp = orange.Value(atts[att_idx], orng_sim_rdk_morgan_fps(a, instance))
                                        instance[atts[att_idx]] = tmp
                                        if callBack: 
                                            stepsDone += 1
                                            if not callBack((100*stepsDone)/nTotalSteps): return None
        
                                att_idx += 1        
                                
                elif m == 'rdk_morgan_features_fps':
                        count = 1
                        for a in actives:
                                attname = m + '(active_'+ str(count)+ ')'
                                for j in range(len(newdata)):
                                        instance = newdata[j]
                                        tmp = orange.Value(atts[att_idx], orng_sim_rdk_morgan_features_fps(a, instance))
                                        instance[atts[att_idx]] = tmp
                                        if callBack: 
                                            stepsDone += 1
                                            if not callBack((100*stepsDone)/nTotalSteps): return None

                                att_idx += 1        
                                        
                elif m == 'rdk_atompair_fps':
                        count = 1
                        for a in actives:
                                attname = m + '(active_'+ str(count)+ ')'
                                for j in range(len(newdata)):
                                        instance = newdata[j]
                                        tmp = orange.Value(atts[att_idx], orng_sim_rdk_atompair_fps(a, instance))
                                        instance[atts[att_idx]] = tmp                
                                        if callBack: 
                                            stepsDone += 1
                                            if not callBack((100*stepsDone)/nTotalSteps): return None

                                att_idx += 1        
        
                elif m == 'azo_pharmacophore_fps':
                        count = 1
                        for a in active_ids:
                                attname = m + '(active_'+ str(count)+ ')'
                                for j in range(len(newdata)):
                                            instance = newdata[j]
                                            tmp = orange.Value(atts[att_idx], azo_pharmacophore_az_inhouse(a, instance, pharmacophore_file))
                                            instance[atts[att_idx]] = tmp
                                            if callBack: 
                                                stepsDone += 1
                                                if not callBack((100*stepsDone)/nTotalSteps): return None

                                att_idx += 1
                                                                
        return newdata
        

def azo_pharmacophore_az_inhouse(active_id, train_instance, pharmacophore_file):
        """ calculate the pharmacophore fingerprint similarity using the AZ inhouse calculated pharmacophore fp
            (the fps are read from a text file for first implementation convenience)
            input are the smiles string and a orange data instance
            returned is a similarty value
        """
        cidName = getCIDAttr(train_instance)
        if not cidName: return None
        train_id = str(int(train_instance[cidName].value))
        
        #print "act"
        fp_A = getPharmacophoreFP(active_id, pharmacophore_file)
        #print "train " + str(train_id)
        fp_T = getPharmacophoreFP(train_id, pharmacophore_file)
        if (fp_A == None or fp_T == None):
                print "Couldn't calc both FPs"
        else:
                sim = getContinuousTanimoto(fp_A,fp_T)

        return sim



def getPharmacophoreFP(mol_id, pharmacophore_file):
        """ extracts the pharmacophore fingerprint in orange fingerprint format from 
            the AZ in-house pharmacophore text file via a mol id match
        """
        pf = open(pharmacophore_file, 'r')
        fp_vals = {}
        for line in pf:
                splitlist = str(line.strip()).split(" ")
                # remove CID, smiles string and bit count
                cid = splitlist.pop(0)
                splitlist.pop(0)
                if (cid.strip() == mol_id.strip()):
                        #print "mol found"
                        splitlist.pop(0)
                        for bit in range(len(splitlist)):
                                if bit % 2 == 0.0:
                                        fp_vals[splitlist[bit]] = splitlist[bit+1]
                        break
        #print fp_vals                        
            
        pf.close()    
        return fp_vals 



def getContinuousTanimoto(fp_A, fp_B):
        """ calculate the Tanimoto coefficient for countinuous valued fingerprints
            according to:
            sim(a,b) = sum_i^N x_a*x_b   /  sum x_a^2 + sum x_b^2 - sum x_a*x_b
            fp_A and fp_B are dictionaries with key = bit number and value = bit value
            if the bit is set to 0 no key value pair is assumed to be set        
        """
        sum_b = 0.0
        sum_a = 0.0
        sum_c = 0.0

        for bit,value in fp_A.iteritems():
                sum_a = sum_a + int(value)**2
        
        for bit_b,value_b in fp_B.iteritems():
                sum_b = sum_b + int(value_b)**2
                if (bit_b in fp_A):
                        sum_c = sum_c + (int(value_b) * int(fp_A[bit_b]))

        sim = sum_c / (sum_a + sum_b - sum_c)
        #print "A: " + str(sum_a) + " B: " + str(sum_b) + " C: " + str(sum_c) + " SIM: " + str(sim)

        return sim


def orng_sim_rdk_topo_fps(smile_active, train_instance):
        """ calculate the fingerprint similarity using the RDK topological fingerprints
                (The fingerprinting algorithm used is similar to that used in the Daylight fingerprinter)
                input are a smiles string and a orange data instance
                returned is a similarity value
        """
        smilesName = getSMILESAttr(train_instance)
        if not smilesName: return None
        smile_train = str(train_instance[smilesName].value)
        
        molAct = rdk.Chem.MolFromSmiles(smile_active)
        molTrain = rdk.Chem.MolFromSmiles(smile_train)
        
        if not molAct: return None
        if not molTrain: return None
        
        fp_A = FingerprintMols.FingerprintMol(molAct)
        fp_T = FingerprintMols.FingerprintMol(molTrain)
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
        
        molAct = rdk.Chem.MolFromSmiles(smile_active)
        molTrain = rdk.Chem.MolFromSmiles(smile_train)
    
        if not molAct: return None
        if not molTrain: return None
    
        fp_A = rdk.Chem.MACCSkeys.GenMACCSKeys(molAct)
        fp_T = rdk.Chem.MACCSkeys.GenMACCSKeys(molTrain)
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
        
        molAct = rdk.Chem.MolFromSmiles(smile_active)
        molTrain = rdk.Chem.MolFromSmiles(smile_train)
    
        if not molAct: return None
        if not molTrain: return None
        
        fp_A = rdk.AllChem.GetMorganFingerprint(molAct,2)
        fp_T = rdk.AllChem.GetMorganFingerprint(molTrain,2)
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
        
        molAct = rdk.Chem.MolFromSmiles(smile_active)
        molTrain = rdk.Chem.MolFromSmiles(smile_train)
    
        if not molAct: return None
        if not molTrain: return None
    
        fp_A = rdk.AllChem.GetMorganFingerprint(molAct,2,useFeatures=True)
        fp_T = rdk.AllChem.GetMorganFingerprint(molTrain,2,useFeatures=True)
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
        
        molAct = rdk.Chem.MolFromSmiles(smile_active)
        molTrain = rdk.Chem.MolFromSmiles(smile_train)
    
        if not molAct: return None
        if not molTrain: return None
    
        fp_A = Pairs.GetAtomPairFingerprint(molAct)
        fp_T = Pairs.GetAtomPairFingerprint(molTrain)
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


def getCIDAttr(data):
        cidName = None
        # "PUBCHEM_CID"
        for attr in [a.name for a in data.domain] + [a.name for a in data.domain.getmetas().values()]:
                if attr in ['"PUBCHEM_CID"',"PUBCHEM_CID", "CID", '"CID"']:
                        cidName = attr
        if not cidName:
                print "Warning: The data set does not contain any known compound identifier"
                print "No pharmacophoric descriptors added!"
                return None
                
        else:
                return cidName
                
                               
def getSMILESAttr(data):
    # Check that the data contains a SMILES attribute
    smilesName = None
    for attr in [a.name for a in  data.domain] + [a.name for a in data.domain.getmetas().values()]:
        if attr in AZOC.SMILESNAMES:
            smilesName = attr
    if not smilesName:
        print "Warning: The data set does not contain any known smiles attribute!"
        print "No similarity descriptors added!"
        return None
    else:       
        return smilesName
