import os,copy
import string
import math
import time

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
        # Create an Orange ExampleTable with a smiles attribute
        smilesAttr = orange.StringVariable("SMILEStoPred")  
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
        errorDesc = "" 
        while nTry > 0:
           try:
                traceLog = "Model Location:"+str(self.modelLocation)+"\n"
                nBadEx = 0        
                # Determine Signature and non-Signature descriptor names
                cinfonyDesc, clabDesc, signatureHeight, bbrcDesc = self.getDescTypes(descList)
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

    def processSignificance(self, smi, prediction, orderedDesc, res, resultsPath, idx = 0):
        """descs* = ["LogP","[So2]"]
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
        """
        atomColor = None
        predictionIsGood = None
        
        DiscMSD = []
        ContMSD = []
        #Define the rules to choose from UP or DOWN vectors                          and the color to highlight
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
            # Choose from proper vector according Documentation
            #if the  labels at the model are oposite to the self.predictionOutcomes, switch
            if (prediction == theGoodPred and not outComeIsRev) or (prediction == theBadPred and outComeIsRev):
                DiscMSD = orderedDesc["Discrete"]["UP"]
                ContMSD = orderedDesc["Continuous"]["DOWN"]
            elif (prediction == theBadPred and not outComeIsRev) or (prediction == theGoodPred and outComeIsRev):
                DiscMSD = orderedDesc["Discrete"]["DOWN"]
                ContMSD = orderedDesc["Continuous"]["UP"]
            else:
                print "ERROR: the precicted value (",prediction,") is not a known predictString"
            
            if prediction == theGoodPred:
                atomColor = 'g'
            else:
                atomColor = 'r'
        else:
            if self.significanceThreshold is None:
                print "WARNING: Cannot process Significance, Missing definition of significanceThreshold for the EndPoint"
                return
            # Choose from proper vector according Documentation
            if prediction < self.significanceThreshold:    # It is a GOOD prediction
                DiscMSD = orderedDesc["Discrete"]["UP"]
                ContMSD = orderedDesc["Continuous"]["DOWN"]
                atomColor = 'g'
            else:
                DiscMSD = orderedDesc["Discrete"]["DOWN"]
                ContMSD = orderedDesc["Continuous"]["UP"]
                atomColor = 'r'

        res["signature"] = {"Discrete":[], "Continuous":[]}
        res["non-signature"] = {"Discrete":[], "Continuous":[]}


        #Discrete attributes
        # find the first signature and non-signature descriptors
        groupNsign = None
        groupNnon_sign = None
        for n,descList in enumerate(DiscMSD):
            for desc in descList:
                if "sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None:
                    if groupNsign is None:
                        groupNsign = n
                elif groupNnon_sign is None:
                    groupNnon_sign = n
            if groupNnon_sign is not None and groupNsign is not None:
                break
        # Group the Discrete attributes signatures and non-signatures
        if groupNsign is not None:
            for desc in DiscMSD[groupNsign]:
                if "sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None:
                    res["signature"]["Discrete"].append(desc)
        if groupNnon_sign is not None:
            for desc in DiscMSD[groupNnon_sign]:
                if not ("sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None):
                    res["non-signature"]["Discrete"].append(desc)


        # Continuous attributes
        # find the first signature and non-signature descriptors
        groupNsign = None
        groupNnon_sign = None
        for n,descList in enumerate(ContMSD):
            for desc in descList:
                if "sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None:
                    if groupNsign is None:
                        groupNsign = n
                elif groupNnon_sign is None:
                    groupNnon_sign = n
            if groupNnon_sign is not None and groupNsign is not None:
                break
        # Group the Continuous attributes signatures and non-signatures
        if groupNsign is not None:
            for desc in ContMSD[groupNsign]:
                if "sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None:
                    res["signature"]["Continuous"].append(desc)
        if groupNnon_sign is not None:
            for desc in ContMSD[groupNnon_sign]:
                if not ("sign" in DescMethodsAvailable and getSignatures.getSignatureHeight(desc) is not None):
                    res["non-signature"]["Continuous"].append(desc)


        #Fix the significant descriptors so that it is a formated string
        MSDstr = ""
        if len(res["non-signature"]["Continuous"]):
            MSDstr += "Continuous: \n"
            for attr in res["non-signature"]["Continuous"]:
                MSDstr += "   "+attr+"\n"
            MSDstr += "\n"
        if len(res["non-signature"]["Discrete"]):
            MSDstr += "Discrete: \n"
            for attr in res["non-signature"]["Discrete"]:
                MSDstr += "   "+attr+"\n"
            MSDstr += "\n"
        res["non-signature"] = MSDstr



        # Most probably Signatures will always be associated with Discrete attributes. Nevertheless, it happens that some are Continuous, and therefore
        #  we will be using signatures reported as Continuous if any
        if not len(res["signature"]["Discrete"]) and not len(res["signature"]["Continuous"]):
            res["imgPath"] = ""
            res["signature"] = ""
            return
        if len(res["signature"]["Discrete"]):
            signature = res["signature"]["Discrete"][0]
        else:
            signature = res["signature"]["Continuous"][0] 
        if resultsPath and os.path.isdir(resultsPath):
            imgPath = os.path.join(resultsPath,"significance_"+str(idx)+"_"+str(time.time()).replace(".",'')+".png")
        else:
            imgPath = ""
        # Call the method to create the image/mol specifying the color of the hilighted atoms   
        res["imgPath"] , res["molStr"], res["atoms"], res["color"] = self.createSignImg(smi,signature,atomColor,imgPath)
        #Fix the significant descriptors so that it is a formated string
        MSDstr = ""
        if len(res["signature"]["Continuous"]):
            MSDstr += "Continuous: \n"
            for attr in res["signature"]["Continuous"]:
                MSDstr += "   "+attr
                if attr == signature and res["atoms"]:    # Signal with * the signature that is being highlited
                    MSDstr += " (*)"
                MSDstr += "\n"
            MSDstr += "\n"
        if len(res["signature"]["Discrete"]):
            MSDstr += "Discrete: \n"
            for attr in res["signature"]["Discrete"]: 
                MSDstr += "   "+attr
                if attr == signature and res["atoms"]:   # Signal with * the signature that is being highlited
                    MSDstr += " (*)"
                MSDstr += "\n"
            MSDstr += "\n"
        res["signature"] = MSDstr



    def getSDs(self, smi, prediction, resultsPath = "", idx = 0, c_step = None):
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
                 "color"         : []}
        # Calculate the signatures id SMILES
        CLabDesc,signList = self.getClabDescSignList(smi)
        if hasattr(self.model,'getTopImportantVars') and self.exToPred:
            orderedDesc = self.model.getTopImportantVars(self.exToPred[0],0, absGradient = False, c_step = c_step)
            #orderedDesc = {'Discrete':   {'DOWN': [['SELMA_GC_type_058', 'SELMA_GC_type_025'], ['SELMA_GC_type_060'], ...]
            #                              'UP':   [['SELMA_GC_type_011'], ['SELMA_GC_type_026'], ...]},
            #               'Continuous': {'DOWN': [['SELMA_GC_type_002'], ['SELMA_GC_type_053'], ...]
            #                              'UP':   [['SELMA_GC_type_006'], ['SELMA_GC_type_024'], ...]} }                    
            # or None
            # or     {'Discrete': [], 'Continuous': []}

            if orderedDesc and "NA" not in orderedDesc:
                self.processSignificance(smi, prediction, orderedDesc, res, resultsPath, idx = idx)
            else:
                print "Model does not have the information needed to compute the Significance"

        return res
        

if __name__ == "__main__":
    #modelPath = "../../tests/source/data/DescModel.model"  # Just RDK descriptors and RDK Fingerprints
    #modelPath = "../../tests/source/data/BBRC_RDK_RDKFP.model"
    modelPath = "/home/palmeida/Workspace/BBRC/model"
    #modelPath = "/home/palmeida/RFmodel"  # with webel descriptors
    #modelPath = "/home/azorangeLive/liveModels_V2/geneTox/currentModelVersion/liveModel/Model"
    #modelPath = "/home/palmeida/liveModels/NaV1.5/currentModelVersion/liveModel/Model/"
    #modelPath = "/home/palmeida/liveModels/hERG/currentModelVersion/liveModel//Model/C0.model/"
    #modelPath = "/home/palmeida/liveModels/AhRrat/currentModelVersion/liveModel/Model/"
    #smi = "C123C5C(O)C=CC2C(N(C)CC1)Cc(ccc4O)c3c4O5"  # NEG
    #smi = "CCC"  # POS
    smi = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    #testData = orange.ExampleTable("../../tests/source/data/QSAR_10mols.tab")
    #testData = orange.ExampleTable("/home/palmeida/Workspace/BBRC/train.tab")
    #res = []
    #for ex in testData:
    if True:
            #print "-------------------------------------------------------------------"
            #smi = ex["SMILES"].value
            print "SMILES: ",smi
            predictor = AZOrangePredictor(modelPath)
            #Needed for classification
            #predictor.predictionOutcomes = ["NEG","POS"]#["POS", "NEG"]
            predictor.predictionOutcomes = ["1","2"]#["POS", "NEG"]

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
            #res.append(str(significance)+" : "+prediction)
            print

#    for r in res:
#        print r
