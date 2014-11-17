from AZutilities import dataUtilities
import itertools
import string
import random


def getDescComb(data, nAttr, Ndesc):

    # Create a string with integers 0-nAttr
    idxList = []
    for idx in range(nAttr):
        idxList.append(str(idx))
    descStr = string.join(idxList, "")
    print descStr

    # Create a list of lists containing the names of the variables to select
    attrListList = []
    # Loop over all relevant number of descriptors
    for length in range(1, nAttr+1):

        # Create all combinations of length up to length
        combinations = itertools.combinations(descStr, length)
    
        # Select the attributes with the corresponding indecies in data and create a list of attribute name lists 
        for combo in combinations:
            attrList = []
            for idx in combo:
                attrList.append(data.domain.attributes[int(idx)].name)
            attrListList.append(attrList)
            
    
    # Randomly select only Ndesc combinations
    indecies = []
    for idx in range(Ndesc):
        indecies.append(random.randint(0, len(attrListList)))
    attrList = []
    for idx in indecies:
        attrList.append(attrListList[idx])

    return attrList


def descSelection(data, NdescComb):

    nAttr = len(data.domain.attributes)
    print "Number of attributes ", nAttr
    print "Maximum number of desc combinations ", pow(2, nAttr)
    print "Ndesc must be lower than the max number of desc combinations"
    print NdescComb

    # Randomly sample Ndesc combinations
    attrList = getDescComb(data, nAttr, NdescComb)

    # Rank the accuracy of each descriptor by averaging the accuracy of all models including a descriptor

    # Select all descriptors above median accuracy and repeat the random sampling of desc combinations

    return attrList


if __name__ == "__main__":

    dataFile = "trainDataAllEP.txt"
    data = dataUtilities.DataTable(dataFile)

    attrList = ["IT03423_Seq_BF", "hERG_IW_pIC50", "IT03423_BF", "IT03423_perc101_BF", "Caco2_intrinsic", "ACDlogD74", "Conc_QTc", "IT03713_BF", "IT10850_BF", "IT22015_BF", "IT22016_BF"]
    data = dataUtilities.attributeSelectionData(data, attrList)

    NdescComb = 100  # Number of desc combinations to sample in the first iteration
    attrList = descSelection(data, NdescComb)
    print attrList
