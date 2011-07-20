""" Module for structural descriptor calculation
	author: Tobias Girschick; tobias.girschick@in.tum.de
			TUM - I12 (wwwkramer.in.tum.de/girschic)
	dependencies: FreeTreeMiner http://wwwkramer.in.tum.de/research/data_mining/pattern_mining/graph_mining
	to come: gSpan', BBRC and LAST-PM 
"""

import os
import sys
import time
import random
import string

import subprocess
from subprocess import Popen, PIPE
import tempfile
import orange
import orngTest, orngStat
from cinfony import rdk
import AZOrangeConfig as AZOC
from AZutilities import dataUtilities

#FTM = os.path.join(os.environ["HOME"], "OpenTox", "serverFiles", "FTM", "ftm")
#FTM = os.path.join("/home/kgvf414/projects/SimBoostedQSAR/FTM", "ftm")
FTM = os.path.join("", "ftm")


def getStructuralDescResult(data,algo,minsup):
    """ delegate to different algorithm methods 
    """
    if (algo == "FTM"):
        return getFTMDescResult(data,minsup)
    elif (algo == "BBRC"):
        return getBBRCDescResult(data,minsup)



def getBBRCDescResult(data,minSup):
	""" Calculates the structural BBRC descriptors for the data using Fminer with the BBRC plugin (python bindings)
		It expects a relative minimum frequency parameter, an optional chi-squared significance parameter
		and a data attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
		It returns a dataset with the same smiles input variable, and as many variables as the descriptors returned by the toolkit
	"""
	smarts = getBBRCsmartsList(data,minSup)

	newdomain = orange.Domain(data.domain.attributes + smarts, data.domain.classVar)
        newdata = orange.ExampleTable(newdomain, data)

	smilesName = getSMILESAttr(data)
        if not smilesName: return None
			
   	count = 0
	for a in newdata:
        	smile = str(a[smilesName].value)
		m = rdk.Chem.MolFromSmiles(smile)
	        if m is None:
        		count += 1
		        continue

	        for b in range(len(smarts)):
        		patt = rdk.Chem.MolFromSmarts(smarts[b].name)
		        if m.HasSubstructMatch(patt):
                		tmp = orange.Value(smarts[b],1.0)
		        else:
                		tmp = orange.Value(smarts[b],0.0)
		        a[smarts[b]] = tmp

	return newdata 


def getBBRCsmartsList(data,minSup):
	""" Calculates the BBRC class-correlated structural features using Fminer with libbbrc python bindings
	    returned is a list of SMARTS string that describe the features
	    Helper function for getBBRCDescResult()
	"""
	import bbrc

        smilesName = getSMILESAttr(data)
        if not smilesName: return None

        # Constructor for standard settings: 95% significance Bbrclevel, minimum frequency 2, type trees, dynamic upper bound, BBRC.
        Fminer = bbrc.Bbrc()
        #if (chisqSig):
         #   Fminer.SetChisqSig(chisqSig)

        # Pass 'true' here to disable usage of result vector and directly print each fragment to the console (saves memory).     
        Fminer.SetConsoleOut(0)
        # Pass 'true' here to enable aromatic rings and use Kekule notation. IMPORTANT! SET THIS BEFORE CALLING AddCompound()! Same as '-a'. 
        Fminer.SetAromatic(1)
        Fminer.SetMinfreq(minSup)

        # add compounds     IMPORTANT! Do not change settings after adding compounds!
        count = 1
        for a in data:
                smile = str(a[smilesName].value)
                activity = float(a.getclass())
                Fminer.AddCompound(smile, count)
		Fminer.AddActivity(activity, count)
                count += 1

        features = []
        # gather results for every root node in vector instead of immediate output
        for j in range(0, Fminer.GetNoRootNodes()-1):
                result = Fminer.MineRoot(j);
                for i in range(0, result.size()-1):
                        #print result[i]
                        # YAML
                        # - [ smarts,    p_chisq,    occ_list_class1,    occ_list_class2,    ... ]                       
                        start_idx = result[i].find('"') + 1
                        smarts = result[i][start_idx:result[i].rfind('"')]
                        # add new attributes to list
                        features.append(orange.FloatVariable(smarts, numberOfDecimals=1))
                        #m = rdk.Chem.MolFromSmarts(smarts)
	
	return features



def getFTMDescResult(data,minSup):
	""" Calculates the structural FTM descriptors for the data using FTM
		It expects a relative minimum frequency parameter and a data attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
		It returns a dataset with the same smiles input variable, and as many variables as the descriptors returned by the toolkit
	"""
	sdf_mols, temp_occ = makeTempFilesFTM(data)
		
	# start FTM subprocess
	ftm(temp_occ.name, minSup, sdf_mols.name)
	#ftm(temp_occ.name, minSup, 'testftm.sdf')
	
	# parse result file and create new data
	occ_flag = 0
	atts = []
	occ_list = []
	for line in temp_occ:

		if occ_flag == 1:
			line_split = line.split()
			#
			# line_split[0] is the "name" of the feature
			# line_split[1] is the occurrence string for all instances separated by comma
			#

			# add new attributes to list
			atts.append(orange.FloatVariable(line_split[0], numberOfDecimals=1))
			#atts.append(orange.EnumVariable(line_split[0], values=['0','1']))
									
			occ = line_split[1].split(",")
			occ_list.append(occ)
					
		if line.startswith('Occurrences:'):
			occ_flag = 1
	
	newdomain = orange.Domain(data.domain.attributes + atts, data.domain.classVar)
	newdata = orange.ExampleTable(newdomain, data)

	count1 = 1
	for j in range(len(newdata)):
		instance = newdata[j]
		count = 1	
				
		for k in range(len(occ_list)):
			if str(occ_list[k][j]) == '1':
				tmp = orange.Value(atts[k],1.0)
				instance[atts[k]] = tmp
			elif str(occ_list[k][j]) == '0':
				tmp = orange.Value(atts[k],0.0)
				instance[atts[k]] = tmp
							
			count += 1
		
	temp_occ.close()
	sdf_mols.close()
	return newdata


def makeTempFilesFTM(data):
	"""	create temporary files for FTM I/O
		returns to file objects that still have to be closed!
		temp_occ contains the occurrence list of substructures (result)
		sdf_mols is for the FTM input (SDF structure file)
	"""
	
	smilesName = getSMILESAttr(data)
	if not smilesName: return None
	
	temp_occ = tempfile.NamedTemporaryFile()
	sdf_mols = tempfile.NamedTemporaryFile(suffix='.sdf') 
    
    # write structures in sdf format to this file
	count = 0
	w = rdk.Chem.SDWriter(sdf_mols.name)
	#w = rdk.Chem.SDWriter('testftm.sdf')
	count_mols = 0
	for a in data:
		smile = str(a[smilesName].value)
		m = rdk.Chem.MolFromSmiles(smile)
		if m is None: 
			count += 1
			#print "ALERT"
			continue
		# set smiles as molname (just to have any name)	
		m.SetProp("_Name", "-->"+smile)
		w.write(m)

		count_mols+=1
	#print "Number of molecules that could not be read: ", count
	#print "Number of molecules read: ", count_mols

	return sdf_mols, temp_occ
	

def getSMARTSrecalcDesc(data,smarts):
    """ Calculates structural descriptors for test and training data.
		In other words, checks for the substructure occurrence (0/1) in the 
		test or prediction molecules. Uses RDK.
		Expects the test/prediction data and a list of SMARTS strings.
		Returns the data including the new features. 
    """
    smilesName = getSMILESAttr(data)
    if not smilesName: return None
    		
    atts = []
    for attr in smarts:
        atts.append(orange.FloatVariable(attr.name, numberOfDecimals=1))

    newdomain = orange.Domain(data.domain.attributes + atts, data.domain.classVar)
    newdata = orange.ExampleTable(newdomain, data)
    	
    count = 0
    for a in newdata:
        smile = str(a[smilesName].value)
        m = rdk.Chem.MolFromSmiles(smile)
        if m is None: 
            count += 1
            continue
		
        for b in range(len(smarts)):
            patt = rdk.Chem.MolFromSmarts(smarts[b].name)
            if m.HasSubstructMatch(patt):
                tmp = orange.Value(atts[b],1.0)
            else:
                tmp = orange.Value(atts[b],0.0)
            a[atts[b]] = tmp
				
    return newdata
	
	
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
	
	

def cross_validation_plusFTM(data, learners, k, f, att_list):
    """
    Perform k-fold cross validation and add FTM features (minsup = f) in each fold 
    The FTM features for each training fold are recalculated for the test fold (NO FTM run!)
    att_list - is the
    list of attributes that will be removed before learning 
    For reference see also:
    http://orange.biolab.si/doc/ofb/accuracy5.py
    http://orange.biolab.si/doc/ofb/c_performance.htm
    """
    acc = [0.0]*len(learners)
    roc = [0.0]*len(learners)
    selection = orange.MakeRandomIndicesCV(data, folds=k)
    for test_fold in range(k):
        train_data = data.select(selection, test_fold, negate=1)
#        print "len->train: ",
#        print len(train_data)
        # add ftm features 
        train_data_ftm = getFTMDescResult(train_data, f)
	minsupStr = str(f).replace(".","")
        filename = data.name + "_ftm_" + minsupStr + "_" + str(test_fold) + ".tab"
        train_data.save(filename)
        train_scaled = dataUtilities.attributeDeselectionData(train_data_ftm, att_list)
        
        # recalc and add ftm features to test fold
        test_data = data.select(selection, test_fold)
        smarts = train_data_ftm.domain.attributes[len(train_data.domain.attributes):]
        print "# FTM features: ",
        print len(smarts)
        test_data_ftm = getSMARTSrecalcDesc(test_data,smarts)
        test_scaled = dataUtilities.attributeDeselectionData(test_data_ftm, att_list)
                
        classifiers = []
        for l in learners:
            classifiers.append(l(train_scaled))
        acc1 = accuracy(test_scaled, classifiers)
        auroc1 = aroc(test_scaled, classifiers)
        print "%d: %s" % (test_fold+1, acc1)
        print "%d: %s" % (test_fold+1, auroc1)
        for j in range(len(learners)):
            acc[j] += acc1[j]
            roc[j] += auroc1[j]
    for j in range(len(learners)):
        acc[j] = acc[j]/k
        roc[j] = roc[j]/k
    return acc, roc
    

def accuracy(test_data, classifiers):
    """
    Taken from: 
    http://orange.biolab.si/doc/ofb/accuracy5.py
    TBD---other measures---reusable stuff??
    """
    correct = [0.0]*len(classifiers)
    for ex in test_data:
        for i in range(len(classifiers)):
            if classifiers[i](ex) == ex.getclass():
                correct[i] += 1
    for i in range(len(correct)):
        correct[i] = correct[i] / len(test_data)
    return correct	
    
    
def aroc(data, classifiers):
    """
    Taken from: 
    http://orange.biolab.si/doc/ofb/roc.py
    """
    ar = []
    for c in classifiers:
        p = []
        for d in data:
            p.append(c(d, orange.GetProbabilities)[0])
        correct = 0.0; valid = 0.0
        for i in range(len(data)-1):
            for j in range(i+1,len(data)):
                if data[i].getclass() <> data[j].getclass():
                    valid += 1
                    if p[i] == p[j]:
                        correct += 0.5
                    elif data[i].getclass() == 0:
                        if p[i] > p[j]:
                            correct += 1.0
                    else:
                        if p[j] > p[i]:
                            correct += 1.0
        ar.append(correct / valid)
    return ar
	
def ftm(temp_occ, freq, mols):
    """
    Actual method that calls FTM on the command line
    """
    ftm_opt = ' -f ' + str(freq) + ' -o ' + temp_occ	+ ' ' + mols
    cmd = FTM + ftm_opt
    #print cmd
    p = Popen(cmd, shell=True, close_fds=True, stdout=PIPE)
    stdout = p.communicate()
    #	p = subprocess.call(ftm_call, shell=True)
