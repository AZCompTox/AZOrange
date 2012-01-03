import os,string,sys
from AZutilities import dataUtilities 
import random
import bbrc
import orange


VERSION = "2.0"

class BBRC(object):
    def __init__(self):
        self.iniDir = os.getcwd()
        self.runDir = "/tmp/BBRC_"+str(random.random())[2:]
        self.SMI = []
        self.CLASS = []
        self.ID = []
        self.verbose = 0
        self.minsup = 6               # -f   default 2
        #self.chisqSig = -1.0         # -p   default -1 
        self.Backbone = True          # -b
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

        #                          lazar    smarts    pvalue    no_aromatic_wc    silent    nr_hits
        self.MyFminer = bbrc.Bbrc( True,    True,     False,    False,             False,    True)

    def __setBBRCOptions(self):
        if not self.MyFminer:
            print "Missing initialization"
            return
        self.MyFminer.SetBackbone(self.Backbone)
        self.MyFminer.SetConsoleOut(0)
        self.MyFminer.SetAromatic(1)
        self.MyFminer.SetMinfreq(self.minsup)   # same as -f
        #self.MyFminer.SetChisqSig(self.chisqSig)  # same as -p 


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
        if self.verbose: print "Running version ",VERSION
        print "Running BBRC for ", repr(self.MyFminer.GetNoCompounds())," compounds..."
        # gather results for every root node in vector instead of immediate output
        lines = []
        for j in range(0, self.MyFminer.GetNoRootNodes()-1):
            result = self.MyFminer.MineRoot(j)
            for i in range(0, result.size()-1):
                 lines.append( result[i].strip() )
                 print result[i]  # DEBUG
        return lines

    def __parseBBRCoutput(self,res, ctrlDescSet):
        #Parse the results to an orange tab file
        print "Parsing BBRC results. Please wait..."
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
        print "Original domain lenght: ",len(self.data.domain)
        print "New domain lenght     : ",len(newDomain)
        print "\n0%"+" "*98+"100%"
        print "|"+"-"*100+"|"
        sys.stdout.write("|")
        sys.stdout.flush()

        newData = dataUtilities.DataTable(newDomain) 
        for idx,ex in enumerate(self.data):
            newEx = orange.Example(newDomain,ex)
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
        if nCompounds < 100:
            sys.stdout.write("="*(100-nCompounds+1))
        print ""
        return newData

    def getDesc(self, data, ctrlDescSet = None):
        if ctrlDescSet and os.path.isfile(ctrlDescSet):  # Not implemented: TODO
            # if the temporal dataSet is specifyed, p=0 minsup=1
            self.chisqSig = 0.0
            self.minsup = 1
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
def getStructuralDescResult(dataIN, algo, minSupPar, active = None):
    """ delegate to different algorithm methods 
    """
    if (algo == "FTM"):              # Using BBRC without class correlation
        BBRCCalc = BBRC()
        BBRCCalc.minsup = minSupPar
        BBRCCalc.active = active
        BBRCCalc.Backbone = False    # No class correlation
        return BBRCCalc.getDesc(dataIN)
    elif (algo == "BBRC"):
        BBRCCalc = BBRC()
        BBRCCalc.minsup = minSupPar
        BBRCCalc.active = active
        BBRCCalc.Backbone = True
        return BBRCCalc.getDesc(dataIN)
    elif (algo == "LAST-PM"):
        return getFMinerDescResult(data,minsup,algo)




if __name__=="__main__":
    dataIN =  dataUtilities.DataTable("./testSMILES.tab")
    algoPar = "LAST-PM"
    minSupPar = 4 
    active = "2"
    outData = getStructuralDescResult(dataIN,algoPar,minSupPar,active)
    if not outData:
        print "Could not get BBRC descriptors!"
    else:
        outData.save("./outBBRC.tab")



