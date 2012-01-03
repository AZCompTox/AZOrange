import os
import sys
import time
import random
import string

import orange
from AZutilities import dataUtilities

#import userDefined Utilites if it exists
if os.path.isfile(os.path.join( os.environ["AZORANGEHOME"], "azorange","AZutilities","extraUtilities.py")):
    from AZutilities import extraUtilities
 
import AZOrangeConfig as AZOC


toolkitsEnabled = AZOC.cinfonyToolkits   

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

toolkitsDef = {
        "rdk":          {"tag" : "rdk."},
        "obabel":       {"tag" : "obabel."},
        "webel":        {"tag" : "webel."},
        "cdk":          {"tag" : "cdk."},
}


def getSMILESAttr(data):
    # Check that the data contains a SMILES attribute
    smilesName = dataUtilities.getSMILESAttr(data)
    if not smilesName:
        print "Warning: The data set does not contain any known smiles attribute!"
        print "         Expected SMILES attribute names: "+str(AZOC.SMILESNAMES)
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

    myDescList = [desc.replace(toolkitsDef["obabel"]["tag"],"") for desc in descList if toolkitsDef["obabel"]["tag"] in desc]
    if not myDescList: return None
       
    resData = orange.ExampleTable(orange.Domain([data.domain[smilesName]] + [orange.FloatVariable(toolkitsDef["obabel"]["tag"]+name) for name in myDescList],0))
    for ex in data:
        newEx = orange.Example(resData.domain)
        newEx[smilesName] = ex[smilesName]
        mol = obabel.readstring("smi", str(newEx[smilesName].value)) 
        moldesc = mol.calcdesc(myDescList)
        for desc in myDescList:
            newEx[toolkitsDef["obabel"]["tag"]+desc] = moldesc[desc]
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

    myDescList = [desc.replace(toolkitsDef["webel"]["tag"],"") for desc in descList if toolkitsDef["webel"]["tag"] in desc]
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
    resData = orange.ExampleTable(orange.Domain([data.domain[smilesName]] + [orange.FloatVariable(toolkitsDef["webel"]["tag"]+name) for name in varNames],0))
    for ex in data:
        newEx = orange.Example(resData.domain)
        smile = str(ex[smilesName].value)
        newEx[smilesName] =smile
        for desc in results[smile]:
            newEx[toolkitsDef["webel"]["tag"]+desc] = results[smile][desc]
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

    myDescList = [desc.replace(toolkitsDef["cdk"]["tag"],"") for desc in descList if toolkitsDef["cdk"]["tag"] in desc]
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
    resData = orange.ExampleTable(orange.Domain([data.domain[smilesName]] + [orange.FloatVariable(toolkitsDef["cdk"]["tag"]+name) for name in varNames],0))
    for ex in data:
        newEx = orange.Example(resData.domain)
        smile = str(ex[smilesName].value)
        newEx[smilesName] =smile
        for desc in results[smile]:
            newEx[toolkitsDef["cdk"]["tag"]+desc] = results[smile][desc]
        resData.append(newEx)

    return resData
  
 
def getRdkDescResult(data,descList, radius = 1):
    """ Calculates the descriptors for the descList using RDK
        It expects an attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
        It returns a dataset with the same smiles input variable, and as many variables as the descriptors 
       returned by the toolkit
    """
    if "rdk" not in toolkitsEnabled:
        return None
    FingerPrints = False
    smilesName = getSMILESAttr(data) 
    if not smilesName: return None
    
    FP_desc = []
    myDescList = [desc.replace(toolkitsDef["rdk"]["tag"],"") for desc in descList if toolkitsDef["rdk"]["tag"] in desc]
    if not myDescList: return None

    if "FingerPrints" in myDescList:
        FingerPrints = True
        myDescList.remove("FingerPrints")
    if sum(["FP_" in fp for fp in myDescList]):
        tmpDescList = []
        FingerPrints = True
        for attr in myDescList:
            if "FP_" not in attr:
                tmpDescList.append(attr)
            else:
                FP_desc.append(attr)
        myDescList = tmpDescList

    #Get fingerprints in advance
    fingerPrintsAttrs = []
    fingerPrintsRes = {}
    if FingerPrints:
        for ex in data:
            mol = str(ex[smilesName].value)
            try:
                chemMol = rdk.Chem.MolFromSmiles(mol,True)
                if not chemMol:
                    chemMol = rdk.Chem.MolFromSmiles(mol,False)
                fingerPrint = rdk.AllChem.GetMorganFingerprint(chemMol,radius)
                resDict = fingerPrint.GetNonzeroElements()
            except:
                continue
            fingerPrintsRes[mol] = {}
            for ID in resDict:
                count = resDict[ID]
                name = toolkitsDef["rdk"]["tag"]+"FP_"+str(ID)
                if name not in [x.name for x in fingerPrintsAttrs]:
                    fingerPrintsAttrs.append(orange.FloatVariable(name))
                fingerPrintsRes[mol][name] = float(count)
        #Add FP attributes even if there was no reference to it. Models will need it as FP not present, i.e. equal 0.0 !
        for fpDesc in FP_desc:
            name = toolkitsDef["rdk"]["tag"]+fpDesc
            if name not in [str(attr.name) for attr in fingerPrintsAttrs]:
                fingerPrintsAttrs.append(orange.FloatVariable(name))
    #Test attrTypes
    for ex in data:
        try:
             attrObj = []
             molStr = str(ex[smilesName].value)
             chemMol = rdk.Chem.MolFromSmiles(molStr,True)
             if not chemMol:
                chemMol = rdk.Chem.MolFromSmiles(molStr,False)
             mol = rdk.readstring("mol", rdk.Chem.MolToMolBlock(chemMol))
             moldesc = mol.calcdesc(myDescList)
             for desc in myDescList:
		 if type(moldesc[desc]) == str:
                     attrObj.append(orange.StringVariable(toolkitsDef["rdk"]["tag"] + desc))
                 else:
                     attrObj.append(orange.FloatVariable(toolkitsDef["rdk"]["tag"] + desc))

             #Process fingerprints
             if FingerPrints:
                 for desc in [fp for fp in fingerPrintsAttrs if fp.name not in attrObj]:
                     attrObj.append(desc)#orange.FloatVariable(desc.name))
             break
        except:
            continue    


    resData = orange.ExampleTable(orange.Domain([data.domain[smilesName]] + attrObj,0))     
    badCompounds = 0
    for ex in data:
        newEx = orange.Example(resData.domain)   # All attrs: ?, ?, ?, ..., ?
        newEx[smilesName] = ex[smilesName]
        molStr = str(newEx[smilesName].value)
        # OBS - add something keeping count on the number of unused smiles
        try:
             chemMol = rdk.Chem.MolFromSmiles(molStr,True)
             if not chemMol:
                chemMol = rdk.Chem.MolFromSmiles(molStr,False) 
             mol = rdk.readstring("mol", rdk.Chem.MolToMolBlock(chemMol))
             #mol = rdk.readstring("smi", molStr)
             moldesc = mol.calcdesc(myDescList)
             for desc in myDescList:
                 newEx[toolkitsDef["rdk"]["tag"]+desc] = moldesc[desc]
 
             #Process fingerprints
             if FingerPrints:
                 for desc in fingerPrintsAttrs:
                     if desc.name in fingerPrintsRes[molStr]:
                         newEx[desc.name] = fingerPrintsRes[molStr][desc.name]
                     else:
                         newEx[desc.name] = 0.0
             resData.append(newEx)
        except: 
            badCompounds += 1
    print "Compounds in original data:       ",len(data)
    print "Compounds able to calculate descs:",len(resData)
    print "Ignored Compounds:                ",badCompounds

    return resData
 
def getCinfonyDescResults(origData,descList,radius=1):
    """Calculates the cinfony descriptors on origData
       maintains the input variables and class
       Adds the Cinfony descritors 
            Returns a new Dataset"""
    if not origData or not descList: return None
    smilesName = getSMILESAttr(origData)
    if not smilesName: return None
    #Create a new domain saving original smiles and other attributes
    newDomain = orange.Domain([attr for attr in origData.domain if attr is not origData.domain.classVar] + [orange.StringVariable("origSmiles")],origData.domain.classVar)
    data = dataUtilities.DataTable(newDomain, origData)
    # Standardize SMILES
    for ex in data:
        ex["origSmiles"] = ex[smilesName].value
    #TODO: Create a method in dataUtilities to standardize the attribute smilesName in place having the attr origSmiles as ID
    if "AZutilities.extraUtilities" in sys.modules and hasattr(extraUtilities, "StandardizeSMILES"):
         # Call a method for standardizing the SMILES in Data.
         # The method is expected to change the attribute defined as smiAttr in data object
         #                                 +->Data     +-> SMILES attribuite name     +->Compound Name or attribute to act as an ID"
         extraUtilities.StandardizeSMILES(data,      smiAttr = smilesName,           cName="origSmiles") 
    results = []

    # Calculate available descriptors
    res = getObabelDescResult(data,descList)
    if res: results.append(res)
    res = getRdkDescResult(data,descList,radius)
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
    # return None if no results at all 
    if not results:
        return None
    resData = results[0]
    if len(results) > 1:
        for res in results[1:]:
            resData = dataUtilities.horizontalMerge(resData, res, smilesName, smilesName)
    data = dataUtilities.horizontalMerge(data, resData, smilesName, smilesName)
    # Revert the SMILES back to it's original state
    for ex in data:
        ex[smilesName] = ex["origSmiles"]
    #Remove the origSmiles attributes
    data = dataUtilities.DataTable(orange.Domain([attr for attr in data.domain if attr.name != "origSmiles" and attr is not data.domain.classVar],data.domain.classVar),data)
    return data

def getAvailableDescs(descSet = "all"):
    """
    descSet : all, rdk, rdkFP, rdkPhysChem, cdk, webel, obabel
    """
    
    #Get descs from obabel
    if "obabel" in toolkitsEnabled:
        obabelDescs = [toolkitsDef["obabel"]["tag"]+desc for desc in obabel.descs]
    else:
        obabelDescs = []
    #Get descs from RDKit
    if "rdk" in toolkitsEnabled:
        rdkDescs = [toolkitsDef["rdk"]["tag"]+desc for desc in rdk.descs] + [toolkitsDef["rdk"]["tag"]+"FingerPrints"]
        rdkFP = [toolkitsDef["rdk"]["tag"]+"FingerPrints"]
        rdkPhysChem = [toolkitsDef["rdk"]["tag"]+desc for desc in rdk.descs]
    else:
        rdkDescs = [] 
    #Get cdk from CDK
    if "cdk" in toolkitsEnabled:
        cdkDescs = [toolkitsDef["cdk"]["tag"]+desc for desc in cdk.descs]
    else:
        cdkDescs = []

    #Get descs from webel     
    try:
        if "webel" in toolkitsEnabled: 
            webelDescs = [toolkitsDef["webel"]["tag"]+desc for desc in webel.getdescs()]
        else:
            webelDescs = []
    except:
        print "WARNING: Could not get the descriptors from webel!"
        if "webel" in toolkitsEnabled: toolkitsEnabled.remove("webel")
        webelDescs = []

    if descSet == "all":
        cinfonyDesc = obabelDescs + rdkDescs + cdkDescs + webelDescs
    elif descSet == "rdk":
        cinfonyDesc = rdkDescs 
    elif descSet == "rdkFP":
        cinfonyDesc = rdkFP
    elif descSet == "rdkPhysChem":
        cinfonyDesc = rdkPhysChem
    elif descSet == "cdk":
        cinfonyDesc = cdkDescs
    elif descSet == "webel":
        cinfonyDesc = webelDescs
    elif descSet == "obabel":
        cinfonyDesc = obabelDescs
    else:
        print descSet+" Is not an allowed descriptor set specification"
        cinfonyDesc = None

    return cinfonyDesc



