import logging
import os,string,sys
from AZutilities import dataUtilities 
import random
import bbrc
import orange

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
        self.active = "POS"  # Only available for Classification
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

        #                          lazar    smarts    pvalue    no_aromatic_wc    silent                     nr_hits
        self.MyFminer = bbrc.Bbrc( True,    True,     False,    False,            not bool(self.verbose),    True)

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
        if self.active not in self.data.domain.classVar.values:
            print "ERROR: '"+str(self.active)+"' is not part of the class values!"
            return None

        smilesName = dataUtilities.getSMILESAttr(self.data)
        for idx,ex in enumerate(self.data):
            if ex.getclass().value == self.active:
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

    def __parseBBRCoutput(self,res, ctrlDescSet):
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
        if ctrlDescSet:# NOT AVAILABLE TODO
            if not os.path.isfile(ctrlDescSet):
                print "Unable to find Control Descriptors Set: ",ctrlDescSet
            else:
                f = open(ctrlDescSet)
                trainDesc = f.readline().strip().split("\t")[2:-1]
                f.close()
                missingDesc = [x for x in trainDesc if x not in allDesc]
                #Attributes to be desellected
                desAttr = [x for x in allDesc if x not in trainDesc]
                selDesc = [x for x in allDesc if x not in desAttr]     

        #TODO  Add attrs and values to the NEWdata!
        newDomainAttrs = [attr for attr in self.data.domain.attributes] + \
                         [orange.FloatVariable(name) for name in selDesc] + \
                         [self.data.domain.classVar]
        newDomain = orange.Domain(newDomainAttrs)
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

    def getDesc(self, data, ctrlDescSet = None):
        if ctrlDescSet and os.path.isfile(ctrlDescSet):  # Not implemented: TODO
            # if the temporal dataSet is specifyed, p=0 minsup=1
            self.ChisqSig = 0.0
            self.minsup = 1
            if self.verbose: 
                print "Because a controll dataset was specifyed, the following parameters were changed:"
                print "  p was set to 0.0"
                print "  f was set to 1"
        self.data = data
        self.__setBBRCOptions()
        self.__createBBRCInputs()
        res = self.__runBBRC() 
        if res:
            return self.__parseBBRCoutput(res,ctrlDescSet)
        else:
            print "No output from BBRC!"
        

#TopLevel interface
def getStructuralDescResult(dataIN, algo, minSupPar, active = None, verbose = 0):
    """ delegate to different algorithm methods 
    """
    if active is not None:
        activeLabel = active
    else:
        activeLabel = dataIN.domain.classVar.values[0]   # For BBRC the active class can be any since it will only use the "count"

    if (algo == "FTM"):              # Using BBRC without class correlation
        BBRCCalc = BBRC(verbose = verbose)
        BBRCCalc.minsup = minSupPar
        BBRCCalc.active = activeLabel
        #Disanling class correlation
        BBRCCalc.DynamicUpperBound = False 
        BBRCCalc.ChisqSig = 0.0
        BBRCCalc.Backbone = False

        return BBRCCalc.getDesc(dataIN)
    elif (algo == "BBRC"):
        BBRCCalc = BBRC(verbose = verbose)
        BBRCCalc.minsup = minSupPar
        BBRCCalc.active = activeLabel
        return BBRCCalc.getDesc(dataIN)
    elif (algo == "LAST-PM"):
        return getFMinerDescResult(data,minsup,algo)
    else:
        print "Algorithm "+algo+" is unknown!"




if __name__=="__main__":
    #logging.basicConfig(level=logging.DEBUG)
    log.setLevel(logging.DEBUG)
    dataIN =  dataUtilities.DataTable("./testSMILES.tab")
    algoPar = "FTM"
    minSupPar = 4 
    outData = getStructuralDescResult(dataIN,algoPar,minSupPar, verbose = 1)
    if not outData:
        print "Could not get BBRC descriptors!"
    else:
        outData.save("./outBBRC.tab")



