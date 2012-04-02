import logging
import os,string,sys
from AZutilities import dataUtilities 
import random
import bbrc
import orange
from cinfony import rdk


logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


class BBRC(object):
    def __init__(self, verbose = 0):
        self.iniDir = os.getcwd()
        self.runDir = "/tmp/BBRC_"+str(random.random())[2:]
        self.SMI = []
        self.CLASS = []
        self.ID = []
        self.verbose = verbose
        self.minsup = 6               # -f   default 2
        self.Backbone = True          # -b   default  True
        self.ChisqActive = True       #      default True
        self.DynamicUpperBound = True #      default True
        self.ChisqSig = None          # default is -1.0 but cannot assign -1.0 !
        # if the temporal dataSet is specifyed, p=0 f=1
        self.data = None
        self.active = None  # Only available for Classification
        # Constructor parameters:
        #fminer_lazar
        #fminer_smarts
        #fminer_pvalues
        #fminer_no_aromatic_wc
        #fminer_silent
        #fminer_nr_hits

        #FMINER_LAZAR : Produce output in linfrag format which can be used as input to Lazar (e.g. export FMINER_LAZAR=1).
        #FMINER_PVALUES : Produce p-values instead of chi-square values (e.g. export FMINER_PVALUES=1).
        #FMINER_NO_AROMATIC_WC: Disallow aromatic wildcard bonds on aliphatic bonds, when aromatic perception was switched off (
        #                       '-a') (e.g. export FMINER_NO_AROMATIC_WC=1).
        #FMINER_SILENT : Redirect STDERR (debug output) of fminer to local file 'fminer_debug.txt'
        #FMINER_NR_HITS : Display (in the occurrence lists) the number of times each fragment occurs in a molecule.

        try:
            #                          lazar    smarts    pvalue    no_aromatic_wc    silent   fminer_nr_hits  
            self.MyFminer = bbrc.Bbrc( True,    True,     False,    False,            0,       True)
        except:
            #FallBack to default constructor. Env vars must be already set: FMINER_LAZAR 1; FMINER_SMARTS 1; FMINER_PVALUES 0
            self.MyFminer = bbrc.Bbrc()

    def __setBBRCOptions(self):
        if not self.MyFminer:
            print "Missing initialization"
            return

        self.MyFminer.SetDynamicUpperBound(self.DynamicUpperBound)
        if self.ChisqSig is not None:
            self.MyFminer.SetChisqSig(self.ChisqSig)
        self.MyFminer.SetBackbone(self.Backbone)
        self.MyFminer.SetChisqActive(self.ChisqActive)
        self.MyFminer.SetConsoleOut(0)
        self.MyFminer.SetAromatic(1)
        self.MyFminer.SetMinfreq(self.minsup)   # same as -f

    def __createBBRCInputs(self):
        if not self.data:
            print "ERROR: Data must be loaded first!"
            return None
        if self.active and self.active not in self.data.domain.classVar.values:
            print "ERROR: '"+str(self.active)+"' is not part of the class values!"
            return None

        smilesName = dataUtilities.getSMILESAttr(self.data)
        print "SMILES attr detected: ",smilesName
        for idx,ex in enumerate(self.data):
            if not self.active:
                activity = 0    #It is unknown 
            elif ex.getclass().value == self.active:
                activity = 1 
            else:
                activity = 0
            ID = idx+1 # ID is the number of coumpound in self.data which is the number os the example (1 based!)
            self.MyFminer.AddCompound(str(ex[smilesName].value), ID)
            self.MyFminer.AddActivity(activity, ID)


    def __runBBRC(self):
        if self.verbose: print "Running BBRC for "+ str(repr(self.MyFminer.GetNoCompounds())) +" compounds..."
        # gather results for every root node in vector instead of immediate output
        lines = []
        for j in range(0, self.MyFminer.GetNoRootNodes()-1):
            result = self.MyFminer.MineRoot(j)
            for i in range(0, result.size()-1):
                 lines.append( result[i].strip() )
                 #print result[i]  # DEBUG
        return lines

    def __parseBBRCoutput(self,res):
        #Parse the results to an orange tab file
        if self.verbose: print "Parsing BBRC results. Please wait..."
        nCompounds = len(self.data)
        allDesc = []
        allIDs = []
        for line in res:
            allDesc.append(line.split("\t")[0].strip())
            allIDs.append([ int(x) for x in line.split("\t")[1][1:-1].strip().split(" ")])

        # Find the Descriptors that are required to be at the output file, but they are not among allDesc
        missingDesc = []
        desAttr = []
        selDesc = [x for x in allDesc]
        newDomainAttrs = [attr for attr in self.data.domain.attributes] + \
                         [orange.FloatVariable(name) for name in selDesc] 
        newDomain = orange.Domain(newDomainAttrs, self.data.domain.classVar)
        if self.verbose: 
            print "Original domain lenght: ",len(self.data.domain)
            print "New domain lenght     : ",len(newDomain)
            print "\n0%"+" "*98+"100%"
            print "|"+"-"*100+"|"
            sys.stdout.write("|")
            sys.stdout.flush()

        newData = dataUtilities.DataTable(newDomain) 
        for idx,ex in enumerate(self.data):
            newEx = orange.Example(newDomain,ex)
            if self.verbose: 
                if nCompounds < 100:
                    sys.stdout.write("=")
                elif idx%(int(nCompounds/100)) == 0:
                    sys.stdout.write("=")
                sys.stdout.flush()

            ID = idx+1   # ID is the number of coumpound in self.data which is the number os the example (1 based!)
            for dIdx,d in enumerate(selDesc):
                if ID in allIDs[dIdx]:
                    newEx[d] = 1.0
                else:
                    newEx[d] = 0.0
            newData.append(newEx)
        if self.verbose: 
            if nCompounds < 100:
                sys.stdout.write("="*(100-nCompounds+1))
            print ""
        return newData

    def getDesc(self, data):
        self.data = data
        self.__setBBRCOptions()
        self.__createBBRCInputs()
        res = self.__runBBRC() 

        if res:
            originalDesc = [attr.name for attr in self.data.domain]
            newData = self.__parseBBRCoutput(res)
            conflictRen = [attr.name for attr in self.data.domain if attr.name not in originalDesc]
            conflictOrig= [attr for attr in originalDesc if attr not in self.data.domain]
            if conflictRen or conflictOrig:
                print "WARNING!!  It has been detected conflictuous descriptors."
                print "           Descriptors that existed in the input data and were also calculated "
                print "             were renamed."
                print "           Descriptors in conflict: ",conflictOrig," -> ",conflictRen
            calculatedDesc = [attr.name for attr in newData.domain if attr.name not in self.data.domain]
            print "Original data attributes: ",len(self.data.domain)
            print "Calculated descriptors: ",len(calculatedDesc)
            print "  ...of which already in the original data (orig. renamed): ",len([attr for attr in calculatedDesc if attr in originalDesc])
            print "  ...of which were added: ",len([attr for attr in calculatedDesc if attr not in self.data.domain])
            print "Output data attributes: ",len(newData.domain)
            return newData
        else:
            print "No output from BBRC!"
            return None
        

#TopLevel interface
def getSMARTSrecalcDesc(data, smarts):
    """ Calculates structural descriptors for test and training data.
                In other words, checks for the substructure occurrence (0/1) in the 
                test or prediction molecules. Uses RDK.
                Expects the test/prediction data and a list of SMARTS strings.
                Returns the data including the new features. 
    """
    smilesName = dataUtilities.getSMILESAttr(data)
    if not smilesName or type(smarts) != list or not len(smarts): 
        print "Please check the input parameters"
        return None
                
    existingAttrs = [attr for attr in smarts if attr in data.domain]
    if existingAttrs:
        print "The input data cannot contain the smarts to be calculated!"
        return None

    newdomain = orange.Domain(data.domain.attributes + \
                              [orange.FloatVariable(attr, numberOfDecimals=1) for attr in smarts],\
                              data.domain.classVar )
    newdata = orange.ExampleTable(newdomain, data)
       
    for ex in newdata:
        smile = str(ex[smilesName].value)
        mol = rdk.Chem.MolFromSmiles(smile)
        if mol is None: 
            continue
        for smrt in smarts:
            patt = rdk.Chem.MolFromSmarts(smrt)
            if mol.HasSubstructMatch(patt):
                ex[smrt] = 1.0
            else:
                ex[smrt] = 0.0
    return newdata


def getStructuralDescResult(dataIN, algo = "FTM", minSupPar = 2, ChisqSig = None, active = None, verbose = 0, descList = []):
    """ delegate to different algorithm methods 
    """
    if not descList:
        descList = []
    outData = None
    if active is not None:
        activeLabel = active
    else:
        if dataIN.domain.classVar:
            activeLabel = dataIN.domain.classVar.values[0]   # For BBRC the active class can be any since it will only use the "count"
        else:
            activeLabel = None
            
    if (algo == "FTM"):              # Using BBRC without class correlation
        BBRCCalc = BBRC(verbose = verbose)
        BBRCCalc.minsup = minSupPar
        BBRCCalc.active = activeLabel
        #Disanling class correlation
        BBRCCalc.DynamicUpperBound = False 
        BBRCCalc.ChisqSig = 0.0
        BBRCCalc.Backbone = False

        outData = BBRCCalc.getDesc(dataIN)
    elif (algo == "BBRC"):
        BBRCCalc = BBRC(verbose = verbose)
        BBRCCalc.minsup = minSupPar
        BBRCCalc.active = activeLabel
        if ChisqSig is not None:       
            if ChisqSig < 0 or ChisqSig > 1:
                print "ERROR: ChisqSig must be defined between 0 and 1"
                return None 
            BBRCCalc.ChisqSig = ChisqSig
        else:
            BBRCCalc.ChisqSig = 0.95
        outData = BBRCCalc.getDesc(dataIN)
    elif (algo == "LAST-PM"):
        outData = getFMinerDescResult(data,minsup,algo)
    else:
        print "Algorithm "+str(algo)+" is unknown!"
    if not outData:
        return None
    newAttrs = [attr.name for attr in outData.domain if attr.name not in dataIN.domain]
    if descList:
        desAttrs = [attr for attr in newAttrs if attr not in descList]
    else:
        desAttrs = []
    print "Structural descriptors requested: "+str(len(descList) or  "ALL")
    print "Structural descriptors returned: "+str(len(newAttrs)-len(desAttrs))
    if desAttrs:
        outData = dataUtilities.attributeDeselectionData(outData, desAttrs)
    unknownAttrs = [attr for attr in descList if attr not in outData.domain]
    print "Attributes not found among the structural descriptors: ",len(unknownAttrs)," (set to 0.0)"
    outData = dataUtilities.attributeAddData(outData, unknownAttrs, orange.FloatVariable, 0.0)
    return outData


if __name__=="__main__":
    #logging.basicConfig(level=logging.DEBUG)
    log.setLevel(logging.DEBUG)
    dataPath = os.path.join(os.environ["AZORANGEHOME"],"tests/source/data/QSAR_10mols.tab")
    #dataPath = "/home/palmeida/RDK_RF.tab"
    dataIN =  dataUtilities.DataTable(dataPath)
    algoPar = "FTM"#"BBRC"
    minSupPar = 6
    dl = None
    outData = getStructuralDescResult(dataIN,algoPar,minSupPar, ChisqSig = None ,verbose = 1, descList = dl)
    if not outData:
        print "Could not get BBRC descriptors!"
    else:
        outData.save("./outBBRC.tab")
        print "OK!"



