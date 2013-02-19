import string
import AZOrangeConfig as AZOC


DescMethodsAvailable = []

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



def getDescTypes(descList):
        clabDescList = []
        cinfonyDescList = []
        bbrcDesc = []
        signDesc = []
        signatureHeight = None
        cinfonyDescs = getCinfonyDesc.getAvailableDescs()
        cinfonyTags = [tk["tag"] for tk in getCinfonyDesc.toolkitsDef.values()]
        for desc in descList:
            if desc in AZOC.SMILESNAMES:
                continue 
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
                signDesc.append(desc)
                if newSignatureHeight > signatureHeight: signatureHeight = newSignatureHeight
        return cinfonyDescList, clabDescList, signatureHeight, bbrcDesc, signDesc


if __name__=="__main__":
    descs = ["[Car]", "[Npl]","Acid","Base","ACDlogD74","[C3]","z1A2A","ACDlogD65"]
    cinfonyDesc, clabDesc, signatureHeight, bbrcDesc, signDesc = getDescTypes(descs)

    print "cinfonyDesc    : ",cinfonyDesc
    print "cclabDesc      : ",clabDesc
    print "signatureHeight: ",signatureHeight
    print "bbrcDesc       : ",bbrcDesc
    print "signDesc       : ",signDesc




