import os
import sys
import time
import random
import string

import orange
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC


#toolkitsEnabled = ["cdk","rdk","obabel","webel"]        # NOT so stable!!
toolkitsEnabled = ["rdk","obabel","webel"]        # Testing Stability!!
#toolkitsEnabled = ["rdk","webel"]                      # Stable

if "cdk" in toolkitsEnabled:
    try: 
        from cinfony import cdk  
    except:
        print "WARNING: cdk is not available for cinfony."
        toolkitsEnabled.remove("cdk")   # craching with :  A fatal error has been detected by the Java Runtime Environment:
                                        #                  SIGFPE (0x8) at pc=0x89de6617, pid=12239, tid=3078584000
                                        #  Seems that this is coused by openbabel
if "rdk" in toolkitsEnabled:
    try:
        from cinfony import rdk      # crashing with *** glibc detected *** python: double free or corruption (fasttop): 0x0a0f36a8 ***
    except:
        print "WARNING: rdk is not available for cinfony."
        toolkitsEnabled.remove("rdk")

if "obabel" in toolkitsEnabled:
    try:
        from cinfony import obabel   # crashing with Floating exception starting azorange and ** Open Babel Error  in SMARTSError
    except:
        print "WARNING: obabel is not available for cinfony."
        toolkitsEnabled.remove("obabel")

if "webel" in toolkitsEnabled:
    try:
        from cinfony import webel
    except:
        print "WARNING: webel is not available for cinfony."
        toolkitsEnabled.remove("webel")

obabelTag = "obabel."
rdkTag = "rdk."
cdkTag = "cdk."
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
    if "obabel" not in toolkitsEnabled:
        return None
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
    """ Calculates the descriptors for the descList using Webel
        It expects an attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
        It returns a dataset with the same smiles input variable, and as many variables as the descriptors 
       returned by the toolkit
    """
    if "webel" not in toolkitsEnabled:
        return None
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


def getCdkDescResult(data,descList):
    """ Calculates the descriptors for the descList using cdk
        It expects an attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
        It returns a dataset with the same smiles input variable, and as many variables as the descriptors 
       returned by the toolkit
    """
    if "cdk" not in toolkitsEnabled:
        return None
    smilesName = getSMILESAttr(data)
    if not smilesName: return None

    myDescList = [desc.replace(cdkTag,"") for desc in descList if cdkTag in desc]
    if not myDescList: return None
    #Compute the results
    results = {}
    for ex in data:
        smile = str(ex[smilesName].value)
        mol = cdk.readstring("smi", smile)
        results[smile] = mol.calcdesc(myDescList)
    # Get all the different descriptor names returned by cdk
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
    """ Calculates the descriptors for the descList using RDK
        It expects an attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
        It returns a dataset with the same smiles input variable, and as many variables as the descriptors 
       returned by the toolkit
    """
    if "rdk" not in toolkitsEnabled:
        return None
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

    res = getCdkDescResult(data,descList)
    if res: results.append(res)

    # Convert any nan to a '?'
    if len(results):
        for res in results:
            for ex in res:
                for attr in ex.domain:
                    if ex[attr] != ex[attr]:   # Will fail if it is 'nan'
                        ex[attr] = '?'

    resData = results[0]
    if len(results) > 1:
        for res in results[1:]:
            resData = dataUtilities.horizontalMerge(resData, res, smilesName, smilesName)
        
    return dataUtilities.horizontalMerge(data, resData, smilesName, smilesName)
      

def getAvailableDescs():
    #Get descs from obabel
    if "obabel" in toolkitsEnabled:
        obabelDescs = [obabelTag+desc for desc in obabel.descs]
    else:
        obabelDescs = []
    #Get descs from RDKit
    if "rdk" in toolkitsEnabled:
        rdkDescs = [rdkTag+desc for desc in rdk.descs]
    else:
        rdkDescs = [] 
    #Get cdk from CDK
    if "cdk" in toolkitsEnabled:
        cdkDescs = [cdkTag+desc for desc in cdk.descs]
    else:
        cdkDescs = []

    #Get descs from webel     
    try:
        if "webel" in toolkitsEnabled: 
            webelDescs = [webelTag+desc for desc in webel.getdescs()]
        else:
            webelDescs = []
    except:
        print "WARNING: Could not get the descriptors from webel!"
        if "webel" in toolkitsEnabled: toolkitsEnabled.remove("webel")
        webelDescs = []

    cinfonyDesc = obabelDescs + rdkDescs + cdkDescs + webelDescs

    return cinfonyDesc



