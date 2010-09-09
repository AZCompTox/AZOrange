import os
import sys
import time
import random
import string

import orange
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC

from cinfony import rdk      # crashing with *** glibc detected *** python: double free or corruption (fasttop): 0x0a0f36a8 ***
from cinfony import obabel   # crashing with Floating exception starting azorange
from cinfony import webel

obabelTag = "obabel."
rdkTag = "rdk."
webelTag = "webel."


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

def getObabelDescResult(data,descList):
    """ Calculates the descriptors for the descList using obabel
        It expects an attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
        It returns a dataset with the same smiles input variable, and as many variables as the descriptors 
       returned by the toolkit
    """
    smilesName = getSMILESAttr(data)
    if not smilesName: return None

    myDescList = [desc.replace(obabelTag,"") for desc in descList if obabelTag in desc]
    if not myDescList: return None
       
    resData = orange.ExampleTable(orange.Domain([data.domain[smilesName]] + [orange.FloatVariable(name) for name in myDescList],0))
    for ex in data:
        newEx = orange.Example(resData.domain)
        newEx[smilesName] = ex[smilesName]
        mol = obabel.readstring("smi", str(newEx[smilesName].value)) 
        moldesc = mol.calcdesc(myDescList)
        for desc in myDescList:
            newEx[desc] = moldesc[desc]
        resData.append(newEx)
    return resData
   
def getWebelDescResult(data,descList):
    """ Calculates the descriptors for the descList using obabel
        It expects an attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
        It returns a dataset with the same smiles input variable, and as many variables as the descriptors 
       returned by the toolkit
    """
    smilesName = getSMILESAttr(data)
    if not smilesName: return None

    myDescList = [desc.replace(webelTag,"") for desc in descList if webelTag in desc]
    if not myDescList: return None

    #Compute the results
    results = {}
    for ex in data:
        smile = str(ex[smilesName].value)
        mol = webel.readstring("smi", smile)
        results[smile] = mol.calcdesc(myDescList)
    # Get all the different descriptor names returned by webel
    varNames = []
    for res in results:
        for desc in results[res]:
            if desc not in varNames:
                varNames.append(desc)
    # Generate the dataset assuring the same order of examples
    resData = orange.ExampleTable(orange.Domain([data.domain[smilesName]] + [orange.FloatVariable(name) for name in varNames],0))
    for ex in data:
        newEx = orange.Example(resData.domain)
        smile = str(ex[smilesName].value)
        newEx[smilesName] =smile
        for desc in results[smile]:
            newEx[desc] = results[smile][desc]
        resData.append(newEx)

    return resData
   
def getRdkDescResult(data,descList):
    """ Calculates the descriptors for the descList using obabel
        It expects an attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
        It returns a dataset with the same smiles input variable, and as many variables as the descriptors 
       returned by the toolkit
    """
    smilesName = getSMILESAttr(data) 
    if not smilesName: return None
    
    myDescList = [desc.replace(rdkTag,"") for desc in descList if rdkTag in desc]
    if not myDescList: return None

    resData = orange.ExampleTable(orange.Domain([data.domain[smilesName]] + [orange.FloatVariable(name) for name in myDescList],0))     
    for ex in data:
        newEx = orange.Example(resData.domain)
        newEx[smilesName] = ex[smilesName]
        mol = rdk.readstring("smi", str(newEx[smilesName].value))
        moldesc = mol.calcdesc(myDescList)
        for desc in myDescList:
            newEx[desc] = moldesc[desc]
        resData.append(newEx)
    return resData
 
def getCinfonyDescResults(data,descList):
    if not data or not descList: return None
    smilesName = getSMILESAttr(data)
    if not smilesName: return None
        
    results = []

    res = getObabelDescResult(data,descList)
    if res: results.append(res)

    res = getRdkDescResult(data,descList)
    if res: results.append(res)

    res = getWebelDescResult(data,descList)
    if res: results.append(res)

    # Convert any nan to a '?'
    if len(results):
        for res in results:
            for ex in res:
                for attr in ex.domain:
                    if ex[attr] != ex[attr]:   # Will fail if it is 'nan'
                        ex[attr] = '?'


    if len(results) == 1:
        return results[0]
        
    resData = results[0]
    for res in results[1:]:
        resData = dataUtilities.horizontalMerge(resData, res, smilesName, smilesName)
        
    return resData
      

def getAvailableDescs():
    #Get descs from obabel
    obabelDescs = [obabelTag+desc for desc in obabel.descs]
    #Get descs from RDKit
    rdkDescs = [rdkTag+desc for desc in rdk.descs]
    #Get descs from webel     
    try:  
        webelDescs = [webelTag+desc for desc in webel.getdescs()]
    except:
        print "WARNING: Could not get the descriptors from webel!"
        webelDescs = []

    cinfonyDesc = obabelDescs + rdkDescs + webelDescs

    return cinfonyDesc



