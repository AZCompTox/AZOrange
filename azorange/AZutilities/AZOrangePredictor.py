import os,copy
import string
import math
import time

DescMethodsAvailable = []

import orange
try:
    from AZutilities import ClabUtilities as ClabUtilities
    DescMethodsAvailable.append("Clab")
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

    def getDef(self, section, option):
        if not self.modelDef: return ""
        if self.modelDef.has_section(section):
            if self.modelDef.has_option(section,option):
                return self.modelDef.get(section, option)

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
        self.significanceColors = eval(self.getDef("significance","colors"))


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
        cinfonyDesc, clabDesc, signatureHeight, bbrcDesc = self.getDescTypes(descList)
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
        

    def createSignImg(self,smi,signature,atomColor,imgPath):
        colors = []
        print "Creating signature image..."
        if not signature or not atomColor or not smi:
            print "Missing inputs:",str([smi,signature,atomColor])
            return "","",[], []
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
                print "highlightAtoms:",atoms
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
        """When only CLab descriptors are present. """

        # Create an Orange ExampleTable with a smiles attribute
        smilesAttr = orange.StringVariable("SMILEStoPred", values = [smiles])
        myDomain = orange.Domain([smilesAttr], 0)
        self.smilesData = dataUtilities.DataTable(myDomain, [[smiles]])


    def getDescTypes(self, descList):
        clabDescList = []
        cinfonyDescList = []
        bbrcDesc = []
        signatureHeight = None
        cinfonyDescs = getCinfonyDesc.getAvailableDescs()
        cinfonyTags = [tk["tag"] for tk in getCinfonyDesc.toolkitsDef.values()]
        for desc in descList: 
            # Cinfony
            if desc in cinfonyDescs or sum([tag in desc for tag in cinfonyTags]):
                cinfonyDescList.append(desc)
            # None signature CLab
            elif "clab" in DescMethodsAvailable and string.find(desc, "[") != 0:
                if desc == "SMILES" or desc == "ID":
                    print "Warning"
                    print "The model was built with SMILES or IDs!!"
                else:
                    clabDescList.append(desc)
            # BBRCDesc
            elif "bbrc" in DescMethodsAvailable and desc[0] == "[" and desc[1] == "#" and "&" in desc:
                bbrcDesc.append(desc)
            # Signature descriptor by exclusion
            elif "sign" in DescMethodsAvailable:  # else!!
                newSignatureHeight = getSignatures.getSignatureHeight(desc) 
                if newSignatureHeight > signatureHeight: signatureHeight = newSignatureHeight
        return cinfonyDescList, clabDescList, signatureHeight, bbrcDesc


    def getDescriptors(self, smiles):
        self.getSmilesData(smiles)

        # Calculate descriptors defined in the model files
        descList = self.model.varNames
   
        savedSmilesData = dataUtilities.DataTable(self.smilesData)

        #Try 3 time to get All compounds descriptors
        nTry = 3        
        while nTry > 0:
                nBadEx = 0        
                # Determine Signature and non-Signature descriptor names
                cinfonyDesc, clabDesc, signatureHeight, bbrcDesc = self.getDescTypes(descList)
                # Signatures
                if "sign" in DescMethodsAvailable and signatureHeight:
                    print "Calculating signatures...."
                    preCalcData = dataUtilities.DataTable(self.preDefSignatureFile)
                    startHeight = 0                # Not used desc ignored in model prediction
                    endHeight = signatureHeight  
                    self.smilesData  = getSignatures.getSignatures(self.smilesData, startHeight, endHeight, preCalcData)

                # C-Lab desc
                if "clab" in DescMethodsAvailable and clabDesc:
                    print "Calculating C-Lab desc...."
                    self.smilesData = ClabUtilities.appendCLabDesc(clabDesc, self.smilesData)

                # Cinfony
                if cinfonyDesc:
                    print "Calculating Cinfony desc..."
                    self.smilesData = getCinfonyDesc.getCinfonyDescResults(self.smilesData, cinfonyDesc, radius = 5)

                # bbrcDesc
                if "bbrc" in DescMethodsAvailable and bbrcDesc:
                    print "Calculating BBRC desc..."
                    self.smilesData = getBBRCDesc.getBBRCDescResult(self.smilesData, algo = "FTM", minSupPar = 1, descList = bbrcDesc)

                # Detect if the descripts calaculation or something else went wrong!
                for ex in self.smilesData:
                   if sum([ex[attr].isSpecial() for attr in self.smilesData.domain.attributes]) == len(self.smilesData.domain.attributes):
                        nBadEx +=1
                if nBadEx:
                    print "WARNING: Desc. Calculation: From the "+str(len(self.smilesData))+" compounds, "+str(nBadEx)+" could not be calculated!"
                    print "WARNING:   Tying again..."
                    self.smilesData = dataUtilities.DataTable(savedSmilesData)
                    nTry -= 1
                else:
                    nTry = 0

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
        except:
            prediction = None

        return prediction

    def processSignificance(self, smi, prediction, descsUP, descsDOWN, res, resultsPath, idx = 0):
        """descs* = ["LogP","[So2]"]
           res =  { "signature"     : "",       
                    "imgPath"       : "",      for placing the results 
                    "non-signature" : "",
                    "molStr"        : "",
                    "atoms"         : []
                    "color"         : [(r,g,b),(),...]}
        """
        atomColor = None
        #Define the rules to choose from descsUP or descsDOWN and the color to highlight
        if self.model.classVar.varType == orange.VarTypes.Discrete:
            if self.predictionOutcomes is None:
                print "WARNING: Cannot process Significance, Missing definition of predictionOutcomes for the EndPoint"
                return
            revPOut = list(self.predictionOutcomes)
            revPOut.reverse()
            if [str(p) for p in self.model.classVar.values] == self.predictionOutcomes:
                UP = descsUP
                DOWN = descsDOWN
            elif [str(p) for p in self.model.classVar.values] == revPOut:
                UP = descsDOWN
                DOWN = descsUP
            else:
                print "ERROR: User outcome ordered list is not consistens toth model: ",\
                      self.predictionOutcomes, "<-->",self.model.classVar.values

            if  prediction == self.predictionOutcomes[0]:
                res["signature"] = UP[1]
                res["non-signature"] = UP[0]
                atomColor = 'r'
            elif prediction == self.predictionOutcomes[1]:
                res["signature"] = DOWN[1]
                res["non-signature"] = DOWN[0]
                atomColor = 'g'
            else:
                print "ERROR: the precicted value (",prediction,") is not a known predictString"
        else:
            if self.significanceThreshold is None:
                print "WARNING: Cannot process Significance, Missing definition of significanceThreshold for the EndPoint"
                return
            if prediction >= self.significanceThreshold:
                res["signature"] = descsUP[1]
                res["non-signature"] = descsUP[0]
                atomColor = 'r'
            else:
                res["signature"] = descsDOWN[1]
                res["non-signature"] = descsDOWN[0]
                atomColor = 'g'
        if not res["signature"]:
            res["imgPath"] = ""
            return

        if resultsPath and os.path.isdir(resultsPath):
            imgPath = os.path.join(resultsPath,"significance_"+str(idx)+"_"+str(time.time()).replace(".",'')+".png")
        else:
            imgPath = ""
        # Call the method to create the image/mol specifying the color of the hilighted atoms    
        res["imgPath"] , res["molStr"], res["atoms"], res["color"] = self.createSignImg(smi,res["signature"],atomColor,imgPath)



    def getSDs(self, smi, prediction, resultsPath = "", idx = 0):
        # descs will  have a list containg the first most significant non-signature descriptor, 
        #   and the first most significant signatures descriptor in the respective order. Ex:
        #               ["LogP","[So2]"]
        #               ["","[So2]"]       No non-signatures found
        #               ["LogP",""]        No signatures found
        
        descsUP = ["",""]
        descsDOWN = ["",""]
        res =  { "signature"     : "",       
                 "imgPath"       : "",     
                 "non-signature" : "",
                 "molStr"        : "", 
                 "atoms"         : [],
                 "color"         : (0,0,0)}
        # Calculate the signatures id SMILES
        CLabDesc,signList = self.getClabDescSignList(smi)
        if hasattr(self.model,'getTopImportantVars') and self.exToPred:
            orderedDesc = self.model.getTopImportantVars(self.exToPred[0],0, absGradient = False)
            if not orderedDesc:
                print "Model does not have the information needed to compute the Significance"
                return res
            # find the first signature and non-signature descriptors with positive Gradient
            for desc in orderedDesc["UP"]:
                if "sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None:
                    if not descsUP[1] and desc in signList:
                        descsUP[1] = desc
                elif not descsUP[0]:
                    descsUP[0] = desc
            # find the first signature and non-signature descriptors with negative Gradient
            for desc in orderedDesc["DOWN"]:
                if "sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None:
                    if not descsDOWN[1] and desc in signList:
                        descsDOWN[1] = desc
                elif not descsDOWN[0]:
                    descsDOWN[0] = desc
            self.processSignificance(smi, prediction, descsUP, descsDOWN, res, resultsPath, idx = idx)

        return res
        

if __name__ == "__main__":
    #modelPath = "../../tests/source/data/DescModel.model"  # Just RDK descriptors and RDK Fingerprints
    modelPath = "../../tests/source/data/BBRC_RDK_RDKFP.model"
    #modelPath = "/home/palmeida/RFmodel"  # with webel descriptors
    smi = "C123C5C(O)C=CC2C(N(C)CC1)Cc(ccc4O)c3c4O5"  # NEG
    #smi = "CCC"  # POS

    predictor = AZOrangePredictor(modelPath)
    #Needed for classification
    predictor.predictionOutcomes = ["POS", "NEG"]

    #Needed for Regression
    #predictor.significanceThreshold = 0.4
    
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

