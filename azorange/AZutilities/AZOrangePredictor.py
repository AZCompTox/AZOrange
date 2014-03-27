import os,copy
import string
import math
import time

from pprint import pprint

DescMethodsAvailable = []

import orange
try:
    from AZutilities import ClabUtilities as ClabUtilities
    DescMethodsAvailable.append("clab")
except:
    print "WARNING: ClabUtilities not available!"

try:
    from AZutilities import getSignatures
    DescMethodsAvailable.append("sign")
except:
    print "WARNING: getSignatures not available!"

try:
    from AZutilities import getCinfonyDesc
    DescMethodsAvailable.append("cinf")
except:
    print "WARNING: getCinfonyDesc not available!"

try:
    from AZutilities import getBBRCDesc
    DescMethodsAvailable.append("bbrc")
except:
    print "WARNING: getBBRCDesc not available"

from AZutilities import dataUtilities
from AZutilities import miscUtilities
from AZutilities import descUtilities
from trainingMethods import AZBaseClasses

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Features.FeatDirUtilsRD import findNeighbors
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from rdkit.Chem.Draw import MolDrawing
from collections import defaultdict
from PIL import Image
import ConfigParser
"""
    Module for getting predictions from the AZOrange models defined by a full path.
"""

class AZOrangePredictor:

    highConf = None
    medConf = None
    lowConf = None
    highConfString = "High"
    medConfString = "Medium"
    lowConfString = "Low"
    modelDef = None

    modelDef = None
    # testing levels DEBUG

    significanceThreshold = None # to be used in regression models
    predictionOutcomes =  None   # to be used in Classification: [PosGradComponent, NegGradComponent]


    def NNsearch(_self, _smi):
        pass


    def getDef(self, section, option=None):
        """
                If Option is None, it will return True or False depending on the section being present or not.
        """
        if not self.modelDef: return ""
        if self.modelDef.has_section(section):
            if option is None:
                return True
            if self.modelDef.has_option(section,option):
                return self.modelDef.get(section, option)
        elif option is None:
            return False
        return ""


    def loadDefs(self):
        self.preDefSignatureFile = self.getDataFile(self.modelLocation)

        defFilePath = os.path.join(self.mountPoint,"data/modelDef.ini")
        if not os.path.isfile(defFilePath): return

        self.modelDef=ConfigParser.ConfigParser()
        self.modelDef.read(defFilePath)   

        thresholds = eval(self.getDef("confidence","thresholds"))
        if thresholds:
            self.lowConf = thresholds[0]
            self.highConf = thresholds[1] 
            self.medConf = (0.0 + self.lowConf + self.highConf)/2.0
        else:
            self.lowConf = None
            self.highConf = None
            self.medConf = None

        self.predictionOutcomes =  eval(self.getDef("general","label_def"))
        significanceThreshold = self.getDef("significance","threshold")
        if significanceThreshold:
            self.significanceThreshold = eval(significanceThreshold)
        else:
            self.significanceThreshold = None


    def __init__(self, modelPath):
        self.model = None
        self.mountPoint = None
        self.modelLocation = modelPath
        self.preDefSignatureFile = self.getDataFile(modelPath)

        self.model = AZBaseClasses.modelRead(self.modelLocation) 
        if not self.model:
            print "ERROR: Cannot load model ",modelPath
            return None


    def getClabDescSignList(self, smiles, getMolFile=False):
        # Create an Orange ExampleTable with a smiles attribute
        smilesAttr = orange.EnumVariable("SMILEStoPred", values = [smiles])
        myDomain = orange.Domain([smilesAttr], 0)
        smilesData = dataUtilities.DataTable(myDomain, [[smiles]])
        #    Calculate descriptors defined in the model files
        try:
            descList = self.model.varNames
        except:   # Consensus object different
            attributes = self.model.domain.variables
            descList = []
            for attr in attributes:
                descList.append(attr.name)
        #    Determine Signature and non-Signature descriptor names
        cinfonyDesc, clabDesc, signatureHeight, bbrcDesc, signDesc = descUtilities.getDescTypes(descList)
        #    Signatures
        if "sign" in DescMethodsAvailable and signatureHeight:
            print "Calculating signatures..."
            preCalcData = dataUtilities.DataTable(self.preDefSignatureFile)
            startHeight = 0                # Not used desc ignored in model prediction
            endHeight = signatureHeight
            dataSign,cmpdSignDict, cmpdSignList, sdfStr  = getSignatures.getSignatures(smilesData, startHeight, endHeight, preCalcData, returnAtomID=True)
        else:
            cmpdSignList = [[]]
            sdfStr = ""
        if not getMolFile:
            return (clabDesc,cmpdSignList[0])
        elif not sdfStr:
            return (clabDesc,cmpdSignList[0],"","")
        # create a mol file
        molFile = miscUtilities.generateUniqueFile(desc="NN", ext = "mol")
        file= open(molFile,"w")
        molStr=""
        for line in sdfStr[0]:
            if "$$$$" in line:
                break
            molStr += line
            file.write(line)
        file.close()

        return (clabDesc,cmpdSignList[0],molFile,molStr)


    def getNNAtoms(self, molStr, cAtoms, hight):
        atoms = copy.deepcopy(cAtoms)
        # create an RDKit mol
        mol = Chem.MolFromMolBlock(molStr,True,False)
        if not mol:
            print "Could not create mol for compound "
            return []
        adj = GetAdjacencyMatrix(mol)
        visitedAtoms = []
        for n in range(hight):
          for atom in copy.deepcopy(atoms):
             if atom not in visitedAtoms:
                lNN = findNeighbors(atom,adj)
                visitedAtoms.append(atom)
                for lnn in lNN:
                    if lnn not in atoms:
                        atoms.append(lnn)
        atoms.sort()
        return atoms

        

    def createSignImg(self,smi,signature,atomColor,imgPath, endHeight = None):
        colors = []
        print "Creating signature image..."
        if not signature or not atomColor or not smi:
            print "Missing inputs:",str([smi,signature,atomColor])
            return "","",[], []
        if hasattr(self.model, "specialType") and self.model.specialType == 1:
            # Create an Orange ExampleTable with a smiles attribute
            smilesAttr = orange.EnumVariable("SMILEStoPred", values = [smi])
            myDomain = orange.Domain([smilesAttr], 0)
            smilesData = dataUtilities.DataTable(myDomain, [[smi]])
            preCalcData = None
            startHeight = 0
            dataSign,cmpdSignDict, cmpdSignList, sdfStr  = getSignatures.getSignatures(smilesData, startHeight, endHeight, preCalcData, returnAtomID=True)
            cmpdSignList = cmpdSignList[0]
            CLabDesc = []
            # create a mol file
            tmpFile = miscUtilities.generateUniqueFile(desc="NN", ext = "mol")
            file= open(tmpFile,"w")
            molStr=""
            for line in sdfStr[0]:
                if "$$$$" in line:
                    break
                molStr += line
                file.write(line)
            file.close()
        else: 
            CLabDesc,cmpdSignList, tmpFile, molStr  =  self.getClabDescSignList(smi, getMolFile=True)
        if not cmpdSignList or not tmpFile:
            print "Couldn't get the cmpd list or the mol file"
            return "","",[], []
        # create an RDKit mol
        mol = Chem.MolFromMolFile(tmpFile,True,False)
        if not mol:
            mol = Chem.MolFromMolFile(tmpFile,False,False)
        if not mol:
            print "Could not create mol for: ",smi
            return "","",[], []
        adj = GetAdjacencyMatrix(mol)
        # find the NN
        hights = []
        for i in miscUtilities.Range(0,len(cmpdSignList),mol.GetNumAtoms()):
            hList = cmpdSignList[i:i+mol.GetNumAtoms()]
            if len(hList):
                hights.append(cmpdSignList[i:i+mol.GetNumAtoms()])
       
        atoms = []
        hight = None
        for idx,h in enumerate(hights):
            if signature in h:
                for i,a in enumerate(h):
                    if a == signature:
                        atoms.append(i)
                hight = idx
                break
        if len(atoms) == 0:
            print "ERROR: Could not find the atom for ",signature
            return "signatureNOTfound","",[],[]
        #print "IniAtoms: ",atoms
        visitedAtoms = []
        for n in range(hight):
          for atom in copy.deepcopy(atoms):
             if atom not in visitedAtoms:    
                lNN = findNeighbors(atom,adj)
                visitedAtoms.append(atom)
                for lnn in lNN:
                    if lnn not in atoms: 
                        atoms.append(lnn)
        atoms.sort()
        os.system("rm " + tmpFile)
        #Specify the atom colors
        colors=[atomColor]*len(atoms)

        if not imgPath:
            return "",molStr,atoms,colors 
        try:
                #Draw the image
                MolDrawing.elemDict=defaultdict(lambda : (0,0,0))
                Draw.MolToImageFile(mol,imgPath,size=(300, 300), kekulize=True, wedgeBonds=True, highlightAtoms=atoms)
                #Color the Highlighted atoms with the choosen atomColor.
                # Only using one color
                if atomColor == 'r':
                    rgb = (255,0,0)
                elif atomColor == 'g':
                    rgb = (0,255,0)
                else:
                    rgb = (0,0,255)    #Blue
                    
                img = Image.open(imgPath)
                img = img.convert("RGBA")
                pixdata = img.getdata()
                newData = list()
                for item in pixdata:
                  if item[0] == 255 and item[1] == 0 and item[2] == 0:
                    newData.append(rgb + (255,) )
                  else:
                    newData.append(item)
                img.putdata(newData)
                img.save(imgPath)

                if os.path.isfile(imgPath):
                    return imgPath,molStr,atoms,colors
                else:
                    return "",molStr,atoms,colors
        except:
                return "",molStr,atoms,colors

        
    def getDataFile(self, modelPath):

        files = os.listdir(modelPath)
        infoFile = None
        for file in files:
            if string.find(file, "Saved.tab") != -1 or string.find(file, "ImputeData.tab") != -1 or string.find(file, "trainDomain.tab") != -1 :
                infoFile = file 
        return os.path.join(modelPath, infoFile)


    def getSmilesData(self, smiles):
        # Create an Orange ExampleTable with a smiles attribute
        smilesAttr = orange.StringVariable("SMILEStoPred")  
        myDomain = orange.Domain([smilesAttr], 0)
        self.smilesData = dataUtilities.DataTable(myDomain, [[smiles]])


    def setDescriptors(self, ex):

        self.exToPred = ex


    def getDescriptors(self, smiles):
        self.getSmilesData(smiles)

        # Calculate descriptors defined in the model files
        descList = self.model.varNames
   
        savedSmilesData = dataUtilities.DataTable(self.smilesData)

        #Try 3 time to get All compounds descriptors
        nTry = 3       
        errorDesc = "" 
        while nTry > 0:
           try:
           #if True:
                traceLog = "Model Location:"+str(self.modelLocation)+"\n"
                nBadEx = 0        
                # Determine Signature and non-Signature descriptor names
                cinfonyDesc, clabDesc, signatureHeight, bbrcDesc, signDesc = descUtilities.getDescTypes(descList)
                # Signatures
                if "sign" in DescMethodsAvailable and signatureHeight:
                    traceLog += "Calculating signatures...\n"
                    print "Calculating signatures...."
                    preCalcData = dataUtilities.DataTable(self.preDefSignatureFile)
                    startHeight = 0                # Not used desc ignored in model prediction
                    endHeight = signatureHeight  
                    self.smilesData  = getSignatures.getSignatures(self.smilesData, startHeight, endHeight, preCalcData)

                # C-Lab desc
                if "clab" in DescMethodsAvailable and clabDesc:
                    traceLog += "Calculating C-Lab...\n"
                    print "Calculating C-Lab desc...."
                    self.smilesData = ClabUtilities.appendCLabDesc(clabDesc, self.smilesData)

                # Cinfony
                if cinfonyDesc:
                    traceLog += "Calculating Cinfony...\n"
                    print "Calculating Cinfony desc..."
                    self.smilesData = getCinfonyDesc.getCinfonyDescResults(self.smilesData, cinfonyDesc, radius = 5)

                # bbrcDesc
                if "bbrc" in DescMethodsAvailable and bbrcDesc:
                    traceLog += "Calculating BBRC...\n"
                    print "Calculating BBRC desc..."
                    self.smilesData = getBBRCDesc.getBBRCDescResult(self.smilesData, algo = "FTM", minSupPar = 1, descList = bbrcDesc)

                # Detect if the descripts calaculation or something else went wrong!
                for ex in self.smilesData:
                   if sum([ex[attr].isSpecial() for attr in self.smilesData.domain.attributes]) == len(self.smilesData.domain.attributes):
                        nBadEx +=1
                if nBadEx:
                    traceLog += "WARNING: Desc. Calculation: From the "+str(len(self.smilesData))+" compounds, "+str(nBadEx)+" could not be calculated!\n"
                    print "WARNING: Desc. Calculation: From the "+str(len(self.smilesData))+" compounds, "+str(nBadEx)+" could not be calculated!"
                    print "WARNING:   Tying again..."
                    self.smilesData = dataUtilities.DataTable(savedSmilesData)
                    nTry -= 1
                else:
                    nTry = 0
           #else:
           except Exception, e:
                errorDesc = "Error Calculating Descriptors:;"+traceLog+str(e)+";"
                nTry -= 1
        if errorDesc:
           raise Exception(errorDesc)
        self.exToPred = self.smilesData





    def getClabTasksAndSignatures(self, smiles):
        if "clab" not in DescMethodsAvailable or "sign" not in DescMethodsAvailable:
            return
        self.getSmilesData(smiles) 

        # Signatures
        preCalcData = dataUtilities.DataTable(self.preDefSignatureFile)
        startHeight = 0
        endHeight = 1
        dataSign  = getSignatures.getSignatures(self.smilesData, startHeight, endHeight, preCalcData)

        # C-Lab descriptors
        self.exToPred = ClabUtilities.appendCLabTasks(self.clabTasks, dataSign)

    def getClabTask(self, smiles):
        if "clab" not in DescMethodsAvailable:
            return
        self.getSmilesData(smiles) 

        # C-Lab descriptors
        self.exToPred = ClabUtilities.appendCLabTasks(self.clabTasks, self.smilesData)

    def predict(self):

        try:
            prediction = self.model(self.exToPred[0]).value
        except Exception, e:
            prediction = None
            raise Exception("Could not predict: " + str(e))

        return prediction

    def processSignificance(self, smi, prediction, orderedDesc, res, resultsPath, exWithDesc = None, idx = 0, topN = 1, regMinIsDesired = True):
        """descs* = [(1.3, ["LogP"]), (0.2, ["[So2]", ...]), ...]
           res =  { "signature"     : "",       
                    "imgPath"       : "",      for placing the results 
                    "non-signature" : "",
                    "molStr"        : "",
                    "atoms"         : []
                    "color"         : [(r,g,b),(),...]}
        
           It uses for Classificartion: 
                        self.predictionOutcomes that must define [BADlabel, GOODlabel] in this same order
            and for Regression:
                        self.significanceThreshold for which a GOOD prediction is BELOW the threshold
        
        orderedDesc = { "molStr":''...,                                                              # only on specialType=1
                        "height":2,                                                                  # only on specialType=1
                        "atoms":[1,2,3],                                                             # only on specialType=1
                          'Continuous': {
                          'DOWN':[ [('[F]', -0.008885456983609475), ... ('[F]',-0,0001)],
                                 [('[O3]', -0.007324209285573964)],
                                 [('[C3]([C3][C3])', -0.0047175657931883405)],
                                 [('[C3]', -0.00389763719161594)]],
                           'UP': [[('[Car]([Car][Nar])', 0.009768981717302358)],
                                 [('[HC]([C3])', 0.009135563633559857)]]},
                       'Discrete': {'DOWN': [], 'UP': []}}
 
        """
        atomColor = None
      
        orderedDesc_sign =    {'Continuous': {'DOWN': [], 'UP': []},
                               'Discrete': {'DOWN': [], 'UP': []}}

        orderedDesc_nonSign = {'Continuous': {'DOWN': [], 'UP': []},
                               'Discrete': {'DOWN': [], 'UP': []}}
        

        if hasattr(self.model, "specialType") and self.model.specialType == 1:
            print "This is a special model nr 1: It calculates itself the Significant Signatures"
            endHeight = orderedDesc["height"]
            try:
                molStr = orderedDesc["molStr"]
                atoms = eval(orderedDesc["atoms"])
            except:
                atoms = None
                molStr = None
            if not molStr or type(atoms)!=list or not atoms:
                atoms = None
                molStr = None

            orderedDesc_sign = orderedDesc
        else:
            atoms = None
            endHeight = None
            molStr = None
            cinfonyDesc, clabDesc, signatureHeight, bbrcDesc, signDesc = descUtilities.getDescTypes([attr.name for attr in self.model.domain.attributes])
            for attrType in ['Continuous', 'Discrete']:
                for vector in ['UP','DOWN']:
                    for ord in range(len(orderedDesc[attrType][vector])):
                        signEmpty=True
                        nonSignEmpty=True
                        for attr in orderedDesc[attrType][vector][ord]:
                            if attr[0] in signDesc:
                                if signEmpty:
                                    signEmpty = False
                                    orderedDesc_sign[attrType][vector].append([])
                                orderedDesc_sign[attrType][vector][-1].append(attr)
                            else:
                                if nonSignEmpty:
                                    nonSignEmpty = False
                                    orderedDesc_nonSign[attrType][vector].append([])
                                orderedDesc_nonSign[attrType][vector][-1].append(attr)

        #Process color to use if highlight is used
        outComeIsRev = None
        if self.model.classVar.varType == orange.VarTypes.Discrete:
            if self.predictionOutcomes is None:
                print "WARNING: Cannot process Significance, Missing definition of predictionOutcomes for the EndPoint"
                return
            theGoodPred = str(self.predictionOutcomes[1])
            theBadPred = str(self.predictionOutcomes[0])
            if [str(p) for p in self.model.classVar.values] == self.predictionOutcomes:
                outComeIsRev = False
            elif [str(p) for p in self.model.classVar.values][::-1] == self.predictionOutcomes:
                outComeIsRev = True
            else:
                print "ERROR: User outcome ordered list is not consistens toth model: ",\
                      self.predictionOutcomes, "<-->",self.model.classVar.values

            if prediction == theGoodPred:
                atomColor = 'g'
            else:
                atomColor = 'r'
        else:
            if self.significanceThreshold is None:
                print "WARNING: Cannot process Significance, Missing definition of significanceThreshold for the EndPoint"
                return
            if prediction < self.significanceThreshold:    # It is a GOOD prediction
                atomColor = 'g'
            else:
                atomColor = 'r'
        #Process Signatures
        if exWithDesc: # Precalculated signatures
            # cmpdSignList differ from when it is calc from smiles. However, not used. 
            dataSign, cmpdSignDict, cmpdSignList, sdfStr = self.getSignDataStruct(exWithDesc)
        else:
            # OBS Hard coded for signatures 0 to 1.  
            smilesData = self.getAZOdata(smi)
            if "sign" in DescMethodsAvailable:
                dataSign, cmpdSignDict, cmpdSignList, sdfStr = getSignatures.getSignatures(smilesData, 0, 1, returnAtomID = True, useClabSmiles = False)
            else:
                dataSign, cmpdSignDict, cmpdSignList, sdfStr = None,[{}]*len(smilesData),[[]]*len(smilesData),"" 
        # If signSVM model already returning one sign as the most significant
        if not (hasattr(self.model, "specialType") and self.model.specialType == 1):
            downAbs = 0.0
            rankIdxDown = 0
            elemIdxDown = 0
            if len(orderedDesc_sign["Continuous"]["DOWN"]):
                for rankIdx in range(len(orderedDesc_sign["Continuous"]["DOWN"])):
                    if downAbs != 0.0:
                        break
                    for elemIdx in range(len(orderedDesc_sign["Continuous"]["DOWN"][rankIdx])):
                        # Test that the signature exists in the molecule
                        if orderedDesc_sign["Continuous"]["DOWN"][rankIdx][elemIdx][0] in cmpdSignDict[0].keys():
                            downAbs = abs(orderedDesc_sign["Continuous"]["DOWN"][rankIdx][elemIdx][1])
                            rankIdxDown = rankIdx
                            elemIdxDown = elemIdx
                            break
        else: 
            if len(orderedDesc_sign["Continuous"]["DOWN"]):
                downAbs = abs(orderedDesc_sign["Continuous"]["DOWN"][0][0][1])
            else:
                downAbs = 0.0
                        
        # If signSVM model already returning one sign as the most significant
        if not (hasattr(self.model, "specialType") and self.model.specialType == 1):
            upAbs = 0.0
            rankIdxUp = 0
            elemIdxUp = 0
            if len(orderedDesc_sign["Continuous"]["UP"]):
                for rankIdx in range(len(orderedDesc_sign["Continuous"]["UP"])):
                    if upAbs != 0.0:
                        break
                    for elemIdx in range(len(orderedDesc_sign["Continuous"]["UP"][rankIdx])):
                        # Test that the signature exists in the molecule
                        if orderedDesc_sign["Continuous"]["UP"][rankIdx][elemIdx][0] in cmpdSignDict[0].keys():
                            upAbs = abs(orderedDesc_sign["Continuous"]["UP"][rankIdx][elemIdx][1])
                            rankIdxUp = rankIdx
                            elemIdxUp = elemIdx
                            break
        else: 
            if len(orderedDesc_sign["Continuous"]["UP"]):
                upAbs = abs(orderedDesc_sign["Continuous"]["UP"][0][0][1])
            else:
                upAbs = 0.0

        # Could happen that all derivatives are smaller than epsilon
        if orderedDesc_sign["Continuous"]["UP"] or  orderedDesc_sign["Continuous"]["DOWN"]:
            if topN == 1:   # Preserve syntax for Plato, only one significant signature
                if upAbs > downAbs:
                    if not (hasattr(self.model, "specialType") and self.model.specialType == 1):
                        MSDsign = orderedDesc_sign["Continuous"]["UP"][rankIdxUp][elemIdxUp][0]
                        MSDdv = orderedDesc_sign["Continuous"]["UP"][rankIdxUp][elemIdxUp][1]
                    else:
                        MSDsign = orderedDesc_sign["Continuous"]["UP"][0][0][0]
                        MSDdv = orderedDesc_sign["Continuous"]["UP"][0][0][1]
                elif downAbs > upAbs:
                    if not (hasattr(self.model, "specialType") and self.model.specialType == 1):
                        MSDsign = orderedDesc_sign["Continuous"]["DOWN"][rankIdxDown][elemIdxDown][0]
                        MSDdv = orderedDesc_sign["Continuous"]["DOWN"][rankIdxDown][elemIdxDown][1]
                    else:
                        MSDsign = orderedDesc_sign["Continuous"]["DOWN"][0][0][0]
                        MSDdv = orderedDesc_sign["Continuous"]["DOWN"][0][0][1]
                elif downAbs != 0.0:
                    if not (hasattr(self.model, "specialType") and self.model.specialType == 1):
                        MSDsign = orderedDesc_sign["Continuous"]["DOWN"][rankIdxDown][elemIdxDown][0]
                        MSDdv = orderedDesc_sign["Continuous"]["DOWN"][rankIdxDown][elemIdxDown][1]
                    else:
                        MSDsign = orderedDesc_sign["Continuous"]["DOWN"][0][0][0]
                        MSDdv = orderedDesc_sign["Continuous"]["DOWN"][0][0][1]
                else:
                    MSDsign = None
                    MSDsign = 0
            else:
                MSDsign = []
                MSDdv = []
                for idx in range(topN):
                    if abs(orderedDesc_sign["Continuous"]["DOWN"][0][0][1]) > abs(orderedDesc_sign["Continuous"]["UP"][0][0][1]):
                        MSDsign.append(orderedDesc_sign["Continuous"]["DOWN"][0][0][0])
                        MSDdv.append(orderedDesc_sign["Continuous"]["DOWN"].pop(0)[0][1])
                    else:
                        MSDsign.append(orderedDesc_sign["Continuous"]["UP"][0][0][0])
                        MSDdv.append(orderedDesc_sign["Continuous"]["UP"].pop(0)[0][1])
        else:
            MSDsign = None
            MSDsign = 0

        #Process non-signatures
        if self.model.classVar.varType == orange.VarTypes.Discrete and outComeIsRev:
            UP = "DOWN"
            DOWN = "UP"
        elif self.model.classVar.varType != orange.VarTypes.Discrete and not regMinIsDesired:
            UP = "DOWN"
            DOWN = "UP"
        else:
            UP = "UP"
            DOWN = "DOWN"
        #Process DiscreteAttrs
        MSDnonSign = ""
        nD_DOWN = len(orderedDesc_nonSign["Discrete"][DOWN])
        if nD_DOWN: 
            for n in range(min(topN,nD_DOWN)):
                if topN > 1:
                    MSDnonSign += str(n+1)+": "
                MSDnonSign += string.join(["Change "+x[0] for x in orderedDesc_nonSign["Discrete"][DOWN][n]],'\n')+'\n'

        #Process Continuous attributes 
        for n in range(topN):
                if topN > 1:
                    order = str(n+1)+": "
                else:
                    order = ""
                if len(orderedDesc_nonSign["Continuous"][DOWN]):
                    downAbs = abs(orderedDesc_nonSign["Continuous"][DOWN][0][0][1])
                else:
                    downAbs = 0.0
                if len(orderedDesc_nonSign["Continuous"][UP]):
                    upAbs = abs(orderedDesc_nonSign["Continuous"][UP][0][0][1])
                else:
                    upAbs = 0.0

                if orderedDesc_nonSign["Continuous"][UP] and upAbs >= downAbs:
                    TOPmsd = orderedDesc_nonSign["Continuous"][UP].pop(0)
                    MSDnonSign += order + string.join(["Decrease "+x[0] for x in TOPmsd],'\n')+'\n'
                if orderedDesc_nonSign["Continuous"][DOWN] and downAbs >= upAbs:
                    TOPmsd = orderedDesc_nonSign["Continuous"][DOWN].pop(0)
                    MSDnonSign += order + string.join(["Increase "+x[0] for x in TOPmsd],'\n')+'\n'

        res["non-signature"] = MSDnonSign



        # Most probably Signatures will always be associated with Discrete attributes. Nevertheless, it happens that some are Continuous, and therefore
        #  we will be using signatures reported as Continuous if any
        if not MSDsign:
            res["imgPath"] = ""
            res["signature"] = ""
            res["signarure_deriv_val"] = 0
            return
        if resultsPath and os.path.isdir(resultsPath):
            imgPath = os.path.join(resultsPath,"significance_"+str(idx)+"_"+str(time.time()).replace(".",'')+".png")
        else:
            imgPath = ""
        # Call the method to create the image/mol specifying the color of the hilighted atoms  
        if exWithDesc == None: # Don't use with precalc desc
            if molStr and atoms and endHeight is not None and not imgPath:
                print "Using molStr and atoms from Learner Significant Signature"
                res["imgPath"]=''
                res["molStr"] = molStr
                allAtoms = self.getNNAtoms(molStr, atoms, endHeight)
                res["atoms"] = allAtoms
                res["color"] = [atomColor]*len(allAtoms)
            else:
                res["imgPath"] , res["molStr"], res["atoms"], res["color"] = self.createSignImg(smi,MSDsign,atomColor,imgPath,endHeight)
        #Fix the significant descriptors so that it is a formated string
        res["signature"] = MSDsign
        res["signarure_deriv_val"] = MSDdv 

        return res        


    def getSignDataStruct(self, exToPred):
        
        cmpdSignDict = {}
        cmpdSignList = []
        for attr in exToPred.domain.attributes:
            if string.find(attr.name, "[") == 0:  # Is signature
                if float(exToPred[0][attr].value) > 0.0:  # The molecule has the signature
                    cmpdSignDict[attr.name] = exToPred[0][attr].value
                    cmpdSignList.append(attr.name)
        sdfStr = None
        dataSign = exToPred

        return dataSign, [cmpdSignDict], [cmpdSignList], sdfStr


    def getAZOdata(self, smi):
        """
        Create an orange data set with a smiles attribute
        """
        smilesAttr = orange.StringVariable("SMILES")
        smilesDomain = orange.Domain(smilesAttr, 0)
        smilesData = dataUtilities.DataTable(smilesDomain, [[smi]])
        return smilesData


    def getSDs(self, smi, prediction, exWithDesc = None, resultsPath = "", idx = 0, topN=1, c_step = None, regMinIsDesired = True):
        # descs will  have a list containg the first most significant non-signature descriptor, 
        #   and the first most significant signatures descriptor in the respective order. Ex:
        #               ["LogP","[So2]"]
        #               ["","[So2]"]       No non-signatures found
        #               ["LogP",""]        No signatures found
        
        descsUP = ["",""]
        descsDOWN = ["",""]
        res =  { "signature"          : "",       
                 "signarure_deriv_val": 0 ,
                 "imgPath"            : "",     
                 "non-signature"      : "",
                 "molStr"             : "", 
                 "atoms"              : [],
                 "color"              : []}
        # Calculate the signatures id SMILES
        #CLabDesc,signList = self.getClabDescSignList(smi)
        if hasattr(self.model,'getTopImportantVars') and self.exToPred:
            orderedDesc = self.model.getTopImportantVars(self.exToPred[0],0, absGradient = False, c_step = c_step, getGrad = True)
            #orderedDesc = {'Discrete':   {'DOWN': [['SELMA_GC_type_058', 'SELMA_GC_type_025'], ['SELMA_GC_type_060'], ...]
            #                              'UP':   [['SELMA_GC_type_011'], ['SELMA_GC_type_026'], ...]},
            #               'Continuous': {'DOWN': [['SELMA_GC_type_002'], ['SELMA_GC_type_053'], ...]
            #                              'UP':   [['SELMA_GC_type_006'], ['SELMA_GC_type_024'], ...]} }                    
            # or None
            # or     {'Discrete': [], 'Continuous': []}
            if orderedDesc and "NA" not in orderedDesc:
                res = self.processSignificance(smi, prediction, orderedDesc, res, resultsPath, exWithDesc, idx = idx, topN=topN, regMinIsDesired = regMinIsDesired)
            else:
                print "Model does not have the information needed to compute the Significance"

        return res
        

if __name__ == "__main__":
    #modelPath = "../../tests/source/data/DescModel.model"  # Just RDK descriptors and RDK Fingerprints
    modelPath = "../../tests/source/data/BBRC_RDK_RDKFP.model"
    smi = "CCC" 
    predictor = AZOrangePredictor(modelPath)
    #Needed for classification
    predictor.predictionOutcomes = ["POS", "NEG"]#["NEG","POS"]#["POS", "NEG"]#

    #Needed for Regression
    predictor.significanceThreshold = 5.0
    
    print "========== Running =========="
    print "Calculating descriptors..."
    predictor.getDescriptors(smi)
    print "Predicting...", smi, "..."
    prediction = predictor.predict()
    print prediction
    print "Finding significant descriptors..."
    significance = predictor.getSDs(smi, prediction)
    print "========== Results =========="
    print "Prediction:",prediction
    print "Significant descriptors:",significance
    print

