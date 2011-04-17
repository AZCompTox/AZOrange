import orange
import string
import types
import sys
import os
import time
import commands
import random
from opencv import ml
from opencv import cv
import AZOrangeConfig as AZOC
from AZutilities import miscUtilities
version = 9
verbose = 0

#DC: " " "NA" "ERROR" "."  "?" ->  "?"
#DK: "*" "MIS_VAL" "~" "nan"   ->  "~"
#ALWAYS USE UPPER CASE HERE and when using these lists, transform also the test string to upper to make them case insensitive
DClist = ["","?"," ","NA","N/A","ERROR", ".","NAN"]
DC = "?"
DKlist = ["~","*","MIS_VAL"]
DK = "~"


def getIsMissing(ex):
    """Returns True if any value in the example object is missing."""
    isMissing = False
    for value in ex:
        if value.isSpecial(): isMissing = True
    return isMissing

##scPA
def SeedDataSampler(data, nFolds):
    """ Samples the data for being used in Folds: Pseudo-Random Selected based on a Seed
        It assures that the sampling is constant for the same dataset using the same number of folds
        Outputs the respective fold indices, not the actual data.
        Output example for 3 Folds:
                [ [0 0 1 1 0 1 0 0 0]
                  [0 1 0 0 0 0 1 0 1]
                  [1 0 0 0 1 0 0 1 0] ]
        Where (0s) represents the train set
              (1s) represents the test set
           for each fold
        The seed used will be based on data dimensions and nFolds.
        It is assured that examples are used as test exampels in only one testSet.
        In the case when nEx is not divisible by nFolds, the remaining examples will never 
          be part of a testSet, although they will always be included in the trainSets. 
        Usage sample:
            DataIdxs = dataUtilities.SeedDataSampler(self.data, self.nExtFolds)
            for foldN in range(self.nExtFolds):
                trainData = self.data.select(DataIdxs[foldN],negate=1)
                testData = self.data.select(DataIdxs[foldN]) 
    """
    if not data or not nFolds:
        return None
    nEx = len(data)
    nTestEx = int(nEx/nFolds)
    if not nTestEx:
        return None
    usedReg = [0] * nEx
    foldsIdxs = []
    for n in range(nFolds):
        foldsIdxs.append([0] * nEx)

    random.seed(nEx + nFolds + len(data.domain.attributes))
    for idx in range(nFolds):
        for x in range(nTestEx):
            Zidx = int(random.random() * (usedReg.count(0) - 1))
            #Find the index of the Zidx'th Zero  
            realIdx = [idxReg[0] for idxReg in enumerate(usedReg) if idxReg[1]==0][Zidx]
            usedReg[realIdx] = 1
            foldsIdxs[idx][realIdx] = 1

    return foldsIdxs



def getPossibleMetas(data, checkIndividuality = False):
    """Retuns a list of attributes that sould be considered meta but they were not.
        Criteria for considering an attribute as meta:
          1) a) Attributes which have more than 20 different values 
              AND 
             b) (individuality) - Less than half of the values appear more than once 
                                  (Only used if checkIndividuality == True)
          OR
          2) The attribute is of String Type
    """
    maxValuesFactor = 2   #orange:2 factor Ex: 2 -> nValues/2
    valuesCountLimit = 2  #orange:2 Limit for the Number of occurences of a value. after this number of occurences, the value will count for the maxValuesFactor
    possMeta = []
    for idx,attr in enumerate(data.domain.attributes):   #Go through all Attributes
        # Only test the attrs that are strings or have more than 20 values
        #   note: For the string ones there is no way of knowing how many different values exists
        if (attr.varType == orange.VarTypes.String):
            possMeta.append(attr.name)
        elif hasattr(attr,"values") and (len(attr.values) > 20):
            if checkIndividuality:                
                nVals = len(attr.values)
                valuesCount = [0] * nVals
                for ex in data:
                    if ex[idx].isSpecial():
                        continue
                    valuesCount[int(ex[idx])] += 1
                    if sum([i>=valuesCountLimit for i in valuesCount]) == nVals/maxValuesFactor: 
                        break
                if sum([i>=valuesCountLimit for i in valuesCount]) < nVals/maxValuesFactor:
                    possMeta.append(attr.name)
                #else:
                #    print "Warning: Attribute almost-meta: ",attr.name
            else:
                possMeta.append(attr.name)
    return possMeta               

class DataTable(orange.ExampleTable):
    """DataTable(filename[, classLess=False, noMeta=True, removeDuplicatedVars=False] | domain[, examples] | examples)"""
    def __new__(cls, *argTuple, **kwds):
        # Load/create the dataset with option to all variables to be unique and new if not specifyed differently
        if "createNewOn" not in kwds:
            kwds["createNewOn"] = orange.Variable.MakeStatus.OK 

        # Specifying classLess=True will allow to load a txt data whithout a class var.
        classLess = False
        if "classLess" in kwds:
            if kwds["classLess"] == True:
                classLess = True
            kwds.pop("classLess")

        # Specifying noMeta=False will allow to load a txt and auto gess the meta attributes.
        noMeta = True
        if "noMeta" in kwds:
            if kwds["noMeta"] == False:
                noMeta = False
            kwds.pop("noMeta")
        # Specifying removeDuplicatedVars=True will remove the variable that are exacly the same (Name and values).
        removeDuplicatedVars = False
        if "removeDuplicatedVars" in kwds:
            if kwds["removeDuplicatedVars"] == True:
                removeDuplicatedVars = True
            kwds.pop("removeDuplicatedVars")

        # Load the data with the default ExampleTable object
        attributeLoadStatus = None
        metaAttributeLoadStatus = None
        data = orange.ExampleTable(*argTuple,**kwds)
        #Convert any symbold that are not allowed in the unicode format
        for attr in data.domain:
            for sym in AZOC.CVT_SYM:
                attr.name = attr.name.replace(sym,AZOC.CVT_SYM[sym])
        # Fix the variables that have no name to be 'Unnamed'
        UnnamedVars = fixUnnamedNames(data.domain)

        if hasattr(data,"attributeLoadStatus"): attributeLoadStatus = data.attributeLoadStatus
        if hasattr(data,"metaAttributeLoadStatus"): metaAttributeLoadStatus = data.metaAttributeLoadStatus

        # Remove any duplicated Attributes is flaged to do so. This is done by loading the data with createNewOn as default
        if removeDuplicatedVars:
            (remVars, newDomain) = getDomainWithoutDupSameAttr(data.domain)
            #The assignment made this way prevent the data to be duplicated
            if remVars: 
                data = orange.ExampleTable(newDomain,data)
                if hasattr(data,"metaAttributeLoadStatus"): metaAttributeLoadStatus = data.metaAttributeLoadStatus
        else:
            remVars = {}
                

        #When loading a txt file, make all meta attributes regular attributes 
        if noMeta and len(data.domain.getmetas()) > 0 and argTuple and type(argTuple[0]) == str and os.path.splitext(os.path.split(argTuple[0])[-1])[-1]==".txt":
            NewAttrs = [attr for attr in data.domain.attributes] + \
                       [data.domain[attrID] for attrID in data.domain.getmetas()]
            NewDomain = orange.Domain(NewAttrs,data.domain.classVar)
            data = orange.ExampleTable(NewDomain,data)

        # Specifying classLess=True will allow to load a txt data whithout a class var.
        if classLess:
            NewAttrs = [attr for attr in data.domain]
            NewDomain = orange.Domain(NewAttrs,0)
            for attrID in data.domain.getmetas():
                NewDomain.addmeta(attrID,data.domain[attrID])

            data = orange.ExampleTable(NewDomain,data)

        # Fix the domain in case of duplicated name attributes: DUplicated Names will be renamed
        dupVars = fixDuplicatedNames(data.domain)
        # create the self object based on the new domain without any duplicated names
        self = orange.ExampleTable.__new__(cls, data)
        #Assign the DataTable attributes attributeLoadStatus, metaAttributeLoadStatus, duplicatedFixedVars, unnamedFixedVars, removedVars
        if attributeLoadStatus:
            self.attributeLoadStatus = attributeLoadStatus
        if metaAttributeLoadStatus:
            self.metaAttributeLoadStatus = metaAttributeLoadStatus

        self.duplicatedFixedVars = dupVars
        self.unnamedFixedVars = UnnamedVars
        self.removedVars = remVars
        return self

def fixUnnamedNames(domain):
    # For the variable that have no name, assign to then the name: "Unnamed" to be compatible with 
    #   all methods that require the vars to have a name.
    # NOTE: This function DOES MODIFY the input domain. Use it carefully!!
    UnnAttrs = []
    idx = 0
    for col,attr in enumerate(domain):
        if not str(attr.name).strip():
            attr.name = "Unnamed"+(idx and "_"+str(idx) or "")
            idx += 1
            UnnAttrs.append(col)
    return UnnAttrs

def CvMat2List(cvMatrix):
    listMatrix = []
    for row in range(cvMatrix.rows):
        listMatrix.append([0]*cvMatrix.cols)
        for col in range(cvMatrix.cols):
            listMatrix[row][col] = cvMatrix[row,col]
    return listMatrix

def List2CvMat(pyList, type="CV_32FC1"):
    cvType = eval("cv."+type)
    Mat = cv.cvCreateMat(1,len(pyList),cvType)
    if cvType == cv.CV_32FC1:
        for idx,x in enumerate(pyList):
            Mat[idx] = float(x)
    else:
        for idx,x in enumerate(pyList):
            Mat[idx] = int(x)
    return Mat


def Example2CvMat(ex,varNames,thisVer = True,getMissingMask = False):
    """Converts an example to a CvMat by picking the variables from the order specified in varNames
       The output CvMat is a row matrix
       The Conversion is only made for attributes defined in varNames, the class, if present in the ex, will be ignored
       The len of the CvMat returned is the same as ex.domain.attributes
       If the len of varNames is not in accordance with the len of ex.domain.attributes, None will be returned
       The ex passed in should already have been fixed by ExFix when using it for predictions.
          - The ExFix already fixed the order
    """
    if len(ex.domain.attributes) != len(varNames):
        return None
    # Creates a matrix for the attributes
    mat = cv.cvCreateMat(1,len(ex.domain.attributes),cv.CV_32FC1)   #rows, cols, type
    #Fill the cvMat with data
    if thisVer:
        if getMissingMask or '?' in ex:
            missing_mask = cv.cvCreateMat(1,len(varNames),cv.CV_8UC1)   #rows, cols, type
            cv.cvSetZero(missing_mask)
            for idx in xrange(len(varNames)):
                if ex[idx].isSpecial():
                    cv.cvSetReal1D(missing_mask,idx, 1)
                else:
                    cv.cvSetReal1D(mat,idx,float(ex[idx]))
            return (mat,missing_mask)
        else:
            # This loop will assure that the var order is correct according to the Classifier domain when using this for prediction
            for idx in xrange(len(varNames)):
                # This will convert any discrete attribute value to the proper order in the values list
                #   No problem about converting to float since the CvMat is expecting one float 
                #   If the attribute is Disc. it will use it's repective value
                # We HAVE to be SURE that the ex is already according to the classifyer domain when using this for predictions
                #mat[0,idx] = float(ex[idx])
                cv.cvSetReal1D(mat,idx,float(ex[idx]))  # No difference in speed for 1D vectors
                #cv.cvmSet(mat,0,idx,float(ex[idx]))   
            return mat
    else: # for backCompatibility with older models. This is deprecated and will not include the missing mask
        for idx in xrange(len(varNames)):
            value = ex[idx].value
            if miscUtilities.isNumber(value):
                cv.cvSetReal1D(mat,idx, float(value))
            else:
                cv.cvSetReal1D(mat,idx, 0)
        return mat

def ExampleTable2CvMat(data,unfoldDescreteClasses=False):
    """ Converts an orange ExampleTable to an opencv CvMat
        if unfoldDescreteClasses is True, the responses will be a float and unfolded into as many columns as class values 
        It returns a dict with the matrix and the responses in CvMat objects:
            returned["matrix"]         The data without the response variable
            returned["responses"]      The column with the responsed from the response variable
            returned["varTypes"]       the variable types in the CvMat object including the response var if exist
            returned["varNames"]       Ordered variable names used in the CvMat object
            returned["responseName"]   Ordered variable names used in the CvMat object
    """
    CvMatices = {}
    # Create the vars with attribute names in proper order
    CvMatices["varNames"] = [attr.name for attr in data.domain.attributes]
    if data.domain.classVar:
        CvMatices["responseName"] = [data.domain.classVar.name]
    else:
        CvMatices["responseName"] = []

    # Creates a matrix for the attributes
    CvMatices["matrix"] = cv.cvCreateMat(len(data),len(data.domain.attributes),cv.CV_32FC1)   #rows, cols, type
    CvMatices["varTypes"] = cv.cvCreateMat(1,len(data.domain),cv.CV_8UC1) # define all var types + the output var type

    # Fill the cvMat with data from exampleTable
    # This iterations take almoast 100% of the time spent in this procedure... very slow!
    # using cv.cvmSet makes it 2x faster, but this cannot be used in vectors of type CV_8SC1, there may be a bug
        
    discreteAttrs = []
    for idx,attr in enumerate(CvMatices["varNames"]):
        if data.domain[attr].varType == orange.VarTypes.Discrete:
            discreteAttrs.append(True)
            #if data.domain.classVar.varType != orange.VarTypes.Discrete:  #DEBUGSEGFAULT
            #    cv.cvSet1D(CvMatices["varTypes"],idx, ml.CV_VAR_NUMERICAL)        #DEBUGSEGFAULT  In regression problems the RF segfaults when there are attributes defined as discrete
            #else:                                                         #DEBUGSEGFAULT
            cv.cvSet1D(CvMatices["varTypes"],idx, ml.CV_VAR_CATEGORICAL)
        else:
            discreteAttrs.append(False)
            cv.cvSet1D(CvMatices["varTypes"],idx, ml.CV_VAR_NUMERICAL)

        
    #cv.cvSetZero(CvMatices["matrix"])
    # Convert to Numpy:
    #    numPyData[0]        - All examples and attributes
    #    numPyData[0][i]     - example i
    #    numPyData[0][i][j]  - attribute j of example i
    #    numPyData[1]        - Class variable
    #    numPyData[1][i]     - Class value of example i

    matrix = CvMatices["matrix"]
    if data.hasMissingValues():
        numPyData = data.toNumpyMA()
        missing_data_mask = numPyData[0].mask
        #Create the CvMat objects for the masks
        CvMatices["missing_data_mask"] = cv.cvCreateMat(len(data),len(data.domain.attributes),cv.CV_8UC1)   #rows, cols, type
        cv.cvSetZero(CvMatices["missing_data_mask"])
        orngMatrix = numPyData[0] # this will convert discrete attributes to number being each int the order of the respecive value
        lenVars = range(len(CvMatices["varNames"]))
        for idxEx,ex in enumerate(orngMatrix):
            for idx in lenVars:
                if missing_data_mask[idxEx][idx]:
                    #it will not set this particular value since it is undefined. It's however marked in the missing_data_mask
                    cv.cvSet2D(CvMatices["missing_data_mask"],idxEx,idx,1)
                else:
                    cv.cvmSet(matrix,idxEx,idx,ex[idx])
    else:
        numPyData = data.toNumpyMA() 
        CvMatices["missing_data_mask"] = None    
        orngMatrix = numPyData[0] # this will convert discrete attributes to number being each int the order of the respecive value
        lenVars = range(len(CvMatices["varNames"]))
        for idxEx,ex in enumerate(orngMatrix):
            for idx in lenVars:
                cv.cvmSet(matrix,idxEx,idx,ex[idx])

   #Fill the cvMat Response Var
    if data.domain.classVar:
        if data.domain.classVar.varType == orange.VarTypes.Discrete: 
            if  unfoldDescreteClasses:
                CvMatices["responses"]  = cv.cvCreateMat(len(data),len(data.domain.classVar.values),cv.CV_32FC1) 
                CvMatices["varTypes"][0,-1] = ml.CV_VAR_CATEGORICAL  
            else:
                CvMatices["responses"]  = cv.cvCreateMat(len(data),1,cv.CV_32SC1)  # In classification use CV_32SC1, while in regression use: CV_32FC1
                CvMatices["varTypes"][0,-1] = ml.CV_VAR_CATEGORICAL   # The element corresponding to the Class
        else:
            CvMatices["responses"]  = cv.cvCreateMat(len(data),1,cv.CV_32FC1)  # In classification use CV_32SC1, while in regression use: CV_32FC1
            CvMatices["varTypes"][0,-1] = ml.CV_VAR_NUMERICAL   # The element corresponding to the Class
        orngMatrix = numPyData[1]  # this will convert discrete classes to number being each int the order of the respecive value
        if data.hasMissingValues():
            missing_res_mask = numPyData[1].mask
            CvMatices["missing_res_mask"] = cv.cvCreateMat(CvMatices["responses"].rows,CvMatices["responses"].cols,cv.CV_8UC1)   #rows, cols, type
            cv.cvSetZero(CvMatices["missing_res_mask"])
            if unfoldDescreteClasses and data.domain.classVar.varType == orange.VarTypes.Discrete:
                cv.cvSet(CvMatices["responses"],0)
                for idxEx,classV in enumerate(orngMatrix):
                    if missing_res_mask[idxEx]:
                        cv.cvSet(CvMatices["missing_res_mask"][idxEx],1)
                    else:
                        #Here classV is the index of the class value of the actual class
                        cv.cvmSet(CvMatices["responses"],idxEx,int(classV),1)
            else:
                for idxEx,classV in enumerate(orngMatrix):
                    if missing_res_mask[idxEx]:
                        cv.cvSetReal1D(CvMatices["missing_res_mask"],idxEx, 1)
                    else:
                        cv.cvSetReal1D(CvMatices["responses"],idxEx, classV)   #This returns the index of the values within .values or the classs value when class is continuous
        else:
             CvMatices["missing_res_mask"]  = None
             if unfoldDescreteClasses and data.domain.classVar.varType == orange.VarTypes.Discrete:
                 cv.cvSet(CvMatices["responses"],0)
                 for idxEx,classV in enumerate(orngMatrix):
                     #Here classV is the index of the class value of the actual class
                     cv.cvmSet(CvMatices["responses"],idxEx,int(classV),1)
             else:
                 for idxEx,classV in enumerate(orngMatrix):
                     cv.cvSetReal1D(CvMatices["responses"],idxEx, classV)   #This returns the index of the values within .values or the classs value when class is continuous

    else:
        CvMatices["responses"] = None

    return CvMatices


def CvMat2orangeResponse(value, classVar, foldClass = False):
    """ Converts the 'value' returned from a classification/regression made by a openCV learner
        to a Class value of the passed classVar 
        if foldClass is true, the value is unfolded, and needs to be identifyed based on classVar 
           according to the values of classVar. This is used when a classifier does not handle natively discrete class
    """
    discreteClass = classVar.varType == orange.VarTypes.Discrete
    if value == None:   
        return None
    elif foldClass:
        if discreteClass and max(list(value)) != max(list(value)): return classVar('?')
        return discreteClass and classVar(list(value).index(max(list(value)))) or classVar(value[0])
    elif value != value:   # This will be false if value is "nan"
        return classVar('?')
    else:
        return discreteClass and classVar(int( max(min(value,len(classVar.values)), 0) ) )  or classVar(value)


def getQuickDataSize(dataPath):
    """ Fast procedure to get the size of data. Returns a Dict:
        returned["N_EX"]             - Number of examples
        returned["N_ATTR"]           - Number of Attributes
        returned["N_CLASS"]          - Number of classes:   1 if present, 0 if Classless
        returned["maybeOrangeTab"]   - Flag of Possible well formated orange tab file (True or False)
        returned["discreteClass"]    - Flag indicating the type of class: 1:discrete, 0:continuous, -1: unknown
        returned["metaAndIgnore"]    - Number of variables not take in account because they were marked for ignore, or they were meta attributes, -1: unknown
    """
    NLines = 0
    NVars = 0
    NClass = 0
    PossTab = False
    discreteClass = -1

    if not os.path.isfile(dataPath):
        return {"N_EX":NLines, "N_ATTR":NVars-NClass, "N_CLASS":NClass, "maybeOrangeTab":PossTab, "discreteClass":discreteClass}
    try:
        NLines = int(commands.getstatusoutput("wc -l "+dataPath)[1].strip().split()[0])
        headers = commands.getstatusoutput("head -n 3 "+dataPath)[1].split("\n")
        FirstLine = [x.strip() for x in headers[0].split("\t")]
        SecondLine = [x.strip() for x in headers[1].split("\t")]
        ThirdLine = [x.strip() for x in headers[2].split("\t")]

        NVars = len(ThirdLine)
        className = FirstLine[-1]
        if "class" in ThirdLine:
            NClass = 1
            classIdx = ThirdLine.index("class")
            className = FirstLine[classIdx]
        elif "c" in ThirdLine:
            NClass = 1
            classIdx = ThirdLine.index("c")
            className = FirstLine[classIdx]
        else:
            NClass = 0

        #If it is an orange tab file, header is 3 lines
        NLines -= 3
        PossTab = True
        #If afteralll it is not an orange tab file, add 2 lines (the first must be the var names!)
        metaAndIgnore = 0
        values=["", "i","ignore","c","class","m","meta"]
        for attr in ThirdLine:
            if attr not in values:
                if "-dc" not in attr:
                    NLines += 2
                    PossTab = False
                    break
            else:
                if attr in ["i","ignore","m","meta"]:
                    metaAndIgnore += 1
        if PossTab and NClass == 1:
            NVars -= metaAndIgnore
            if SecondLine[classIdx].lower() in ["c","continuous"]:
                discreteClass = 0
            else:
                discreteClass = 1
        else:
            discreteClass = -1
    except:
        return {"N_EX":NLines, "N_ATTR":NVars-NClass, "N_CLASS":NClass, "maybeOrangeTab":False, "discreteClass":discreteClass,"className":"","metaAndIgnore":-1}
    return {"N_EX":NLines, "N_ATTR":NVars-NClass, "N_CLASS":NClass, "maybeOrangeTab":PossTab, "discreteClass":discreteClass,"className":className,"metaAndIgnore":metaAndIgnore}



def changeClassName(data,newName):
    """ Changes the class var name of 'data' to 'newName'
        if the newName already exists in the 'data' attributes, a numbered suffix will be added Ex: Activity_1
        It creates a new domain based on the data's domain, changes the class name and copy the data 
        in 'data' to a new ExampleTable of the new domain.
        - Returns the new ExampleTable or None if something is missing
        - The input data (domain) is not modified!
    """
    if not newName or not data or not data.domain.classVar:
        return None
    newDomain = data.domain

    #create the new class attribute
    newClassName = newName
    suffix = 1
    while newClassName in [x.name for x in data.domain.attributes]:
        newClassName = newName + "_" + str(suffix)
        suffix = suffix + 1
    # Create the new Domain with the class Var being a NEW attribute
    newDomain = orange.Domain(data.domain.attributes, data.domain.classVar.clone())
    # Create the new data with new Domain
    newData = DataTable(newDomain)
        
    #Fill the new example table by converting the examples to the new domain
    ExampleFix = ExFix(newDomain,None,True)
    for ex in data:
        newData.append(ExampleFix.fixExample(ex))

    # Rename the class var name in the new domain
    newData.domain.classVar.name = newClassName

    return newData


def addDiscreteClass(data, className, value, values=[]):
    """
    Adds a categorical response variable to a data set
    Makes the input 'data' classless, adds a new discrete attribute with name 'className' and makes this attribute the class.
    The value 'value' is set to the new attribute 'className' of all examples in 'data'
    The name of the new class attribute may be the same input 'className', or if the name already existed in the domain,
     it uses the 'className' followed by "_" and a number Ex: "Activity_1"
    This method returns a new data with one more attribute in the domain and the class being the new added attribute
    Optionally, you can specify the values of the new class attribute
    returns None if some error occured
    """
    if not data or not className:
        return None
    if type(className) not in  types.StringTypes:
        return None
    else:
        # Make the class name a String Type
        className = str(className)
    if values==None or type(values) != types.ListType:
        return None
    else:
        for idx,v in enumerate(values):
            if type(v) not in  types.StringTypes:
                return None
            else:
                values[idx] = str(v)
    
    #check if the value is already in values. If not, add it
    if value.upper() in DClist: value = DC
    if value.upper() in DKlist: value = DK
    if value.upper() != DC and value.upper() != DK and value not in values:
        values.append(str(value))
    #create the new class attribute
    if className not in data.domain: 
        newClassName = className
        newClass = orange.EnumVariable(newClassName, values=values)
    else:
        newClass = None
        suffix = 1
        while not newClass:
            newClassName = className + "_" + str(suffix)
            if newClassName in data.domain:
                suffix = suffix + 1
            else:
                newClass = orange.EnumVariable(newClassName, values=values)
    # Create the new Domain
    newDomain = orange.Domain(data.domain,newClass)
    # Replace original existing data with same data but with new domain
    newData = DataTable(newDomain,data)

    # Assign the 'value' to all examples
    for ex in newData:
        ex[newClassName] = str(value)

    return newData



def concatenate(datasets,useFirstAsLeader=False,mergeDomains=True):
    """
        Concatenates the ExampleTables present in the 'datasets' list
        if useFirstAsLeader is True, the domain of the first dataset in the datasets list (the primary Table) 
           will be used and none of the new attributes in the other datasets will be used. The flag 'mergeDomains' 
           will be ignored and the values of the first dataset attributes that does not exist in the other datasets will be '?' 
        if mergeDomains is True, the attributes of all domains will be merged and if False, only the attributes present 
           on all datasets domains  will be used (useFirstAsLeader must be False in order to use mergeDomains)

        Returns a list with first value being the concatenated data and the second a Dict with concatenate results
        (string names of attributes): 
                    [   
                        ConcatenatedExampleTable,
                        {
                           "New1"                 : AttributesInD1NotPresentInD2,
                           "New2"                 : AttributesInD1NotPresentInD1,
                           "Equal"                : AttributesOnBothTablesDomainsAndCompatible,
                           "Conflict"             : DiscreteAttributesOnBothTableDomainsButDiffValues,
                           "Diff"                 : AttributesOnBothTableDomainsButNotCompatible, 
                           "NewDiffAttrNames"     : AttributeNamesGivenToDiffD2Attributes,
                           "classVar"             : AttributeDefinedAsClass
                        }
                    ]
    """

    def localConcatenation(D1,D2,allDatasets,useD1Leader,merge):
        # If we are using D1 as Leader, merge must be True in order to include all attributes from D1 
        # even if the attribute does not exist in D2
        if useD1Leader: merge=True
        attributesNamesD1 = [ x.name for x in D1.domain]
        attributesNamesD2 = [ x.name for x in D2.domain]
        commonAttr=[]           # Attributes with common names
        newAttrD1=[]            # Attributes not in the original D1 domain
        newAttrD2=[]            # Attributes not in the original D2 domain
        conflictAttr=[]         # Discrete Attributes in both Domains but different values 
        diffAttr=[]             # Attributes in both Domains but not compatible at all
        diffAttrLeaderD1=[]     # those attr of diffAttr that have as Leader the D1 Domain
        diffAttrLeaderD2=[]     # those attr of diffAttr that have as Leader the D2 Domain
        equalAttr=[]            # Attributes on both domains representing a unique attribute
        classVar = ""

        # Place attribute names in respective lists        
        for attr in attributesNamesD2:
            if attr not in attributesNamesD1:
                if merge: newAttrD2.append(attr)
            else:
                commonAttr.append(attr)
        if merge:
            for attr in attributesNamesD1:
                if attr not in attributesNamesD2:
                    newAttrD1.append(attr)
        for attr in commonAttr:
            #print attr
            if (D1.domain[attr].varType == D2.domain[attr].varType == orange.VarTypes.String) or (D1.domain[attr].varType == D2.domain[attr].varType == orange.VarTypes.Continuous) or (D1.domain[attr]==D2.domain[attr]):
                equalAttr.append(attr)
            else:
                if D1.domain[attr].varType==orange.VarTypes.Discrete and D2.domain[attr].varType==orange.VarTypes.Discrete:
                    conflictAttr.append(attr) 
                else:
                    diffAttr.append(attr)
        #  Create a list of attributes for the new domain according to the repective created lists and passed flags
        NewDomainAttr=[]
        for attr in equalAttr:
            NewDomainAttr.append(D1.domain[attr])
        for attr in newAttrD1:
            NewDomainAttr.append(D1.domain[attr])
        if not useD1Leader:
            for attr in newAttrD2:
                NewDomainAttr.append(D2.domain[attr])

        # Decide from the incompatible attributes what will be the Leader, if not defined already by useD1Leader
        for attr in diffAttr:
            if useD1Leader: # There is already a Leader!
                diffAttrLeaderD1.append(attr)
                NewDomainAttr.append(D1.domain[attr])
            else:           #Find the Leader among all datasets present in the concatenation
                # create a list of indices of datasets where the attribute is present:
                idxAttrList = [idx for idx,testData in enumerate(allDatasets) if (testData and hasattr(testData,"domain") and (attr in testData.domain))]
                totalExamples=0
                OKCvtToCont = 0
                for idx in idxAttrList:
                    for ex in allDatasets[idx]:
                        totalExamples += 1
                        OKCvtToCont += miscUtilities.isNumber(ex[attr].value) and True #his way assure only 1 is added at a time
                #If >=50% can be securely converted to Continuous, use the continuous attribute
                if OKCvtToCont >= (totalExamples/2):
                    if D1.domain[attr].varType==orange.VarTypes.Continuous:
                        diffAttrLeaderD1.append(attr)
                        NewDomainAttr.append(D1.domain[attr])
                    else:
                        diffAttrLeaderD2.append(attr)
                        NewDomainAttr.append(D2.domain[attr])
                else:
                    if D1.domain[attr].varType==orange.VarTypes.Discrete:
                        diffAttrLeaderD1.append(attr)
                        NewDomainAttr.append(D1.domain[attr])
                    else:
                        diffAttrLeaderD2.append(attr)
                        NewDomainAttr.append(D2.domain[attr])
                
        # For the Conflict Attributes only need to join the values from both domains for the specific attributes
        for attr in conflictAttr:
            allValues = D1.domain[attr].values.native()
            for value in D2.domain[attr].values.native():
                if value not in allValues:
                    allValues.append(value)
            NewDomainAttr.append(orange.EnumVariable(str(attr), values = allValues)) 
        # Find the attribute to use as class (priority to D1 class)
        if D1.domain.classVar:
            classVar = D1.domain.classVar.name
        elif D2.domain.classVar and not useD1Leader:
            classVar = D2.domain.classVar.name
        else:
            classVar = ""
        # Create the new domain with chosen class attribute
        if not useD1Leader:
            #Sort attributes by name
            NewDomainAttr.sort(lambda x,y: cmp(x.name,y.name))
        else:
            #Assure that if useD1Leader, then the Domain must be in same order as D1
            tmp=[None]*len(NewDomainAttr)
            for attr in NewDomainAttr:
                tmp[D1.domain.index(attr.name)] = attr
            NewDomainAttr = tmp 
            
        classVarAttr = [x for x in NewDomainAttr if x.name==classVar]
        
        if not classVarAttr:
            classVar = ""
            newDomain = orange.Domain(NewDomainAttr,0)
        else:
            NewDomainAttr.remove(classVarAttr[0])
            newDomain = orange.Domain(NewDomainAttr+classVarAttr)
            classVar = newDomain.classVar.name
                
        # Create an empty table with new domain
        newTable = DataTable(newDomain)       
        # Insert the examples from DataTable1
        for ex in D1:
            newEx=orange.Example(newDomain) # ->   [?, ?, ?, ?, ?, ... , ?]
            for attr in newAttrD1:
                newEx[attr]=ex[attr].value
            for attr in equalAttr:
                newEx[attr]=ex[attr].value
            for attr in conflictAttr:
                newEx[attr]=ex[attr].value
            for attr in diffAttr:
                if attr in diffAttrLeaderD1: # D1 is the Leader
                    newEx[attr]=ex[attr].value
                else:   # D2 is the Leder
                    if D2.domain[attr].varType==orange.VarTypes.Discrete: # The Leader is Discrete
                        if str(ex[attr].value) in D2.domain[attr].values:
                            newEx[attr] = str(ex[attr].value)
                        else:
                            newEx[attr] = '?'
                    else:                                                  # The Leader is Continuous
                        if miscUtilities.isNumber(ex[attr].value):
                           newEx[attr] = float(ex[attr].value)
                        else:
                            newEx[attr] = '?'
            newTable.append(newEx)
        # Insert the examples from DataTable2
        for ex in D2:
            newEx=orange.Example(newDomain) # ->   [?, ?, ?, ?, ?, ... , ?]
            if not useD1Leader:
                for attr in newAttrD2:
                    newEx[attr]=ex[attr].value
            for attr in equalAttr:
                newEx[attr]=ex[attr].value

            for attr in conflictAttr:
                newEx[attr]=ex[attr].value

            for attr in diffAttr:
                if attr in diffAttrLeaderD2: # D2 is the Leader
                    newEx[attr]=ex[attr].value
                else:   # D1 is the Leder
                    if D1.domain[attr].varType==orange.VarTypes.Discrete: # The Leader is Discrete
                        if str(ex[attr].value) in D1.domain[attr].values:
                            newEx[attr] = str(ex[attr].value)
                        else:
                            newEx[attr] = '?'
                    else:                                                 #The Leader is Continuous
                        if miscUtilities.isNumber(ex[attr].value):
                           newEx[attr] = float(ex[attr].value)
                        else:
                            newEx[attr] = '?'
            newTable.append(newEx)
        return (newTable, diffAttr)

    firstDataIdx = None
    # Use the first valid dataset to start the concatenation
    for idx,data in enumerate(datasets):
        if data:
            newTable = datasets[idx]
            firstDataIdx = idx
            break
    if firstDataIdx == None:
        if verbose >0: print "ERROR: None of the datasets provided are valid!"
        return (None, {})
    elif firstDataIdx != 0:
        if verbose >0: print "WARNING: The first dataset provided is not valid, the data with index",firstDataIdx,"was used as first dataset."
    allDiffAttr = []
    for idx,data in enumerate(datasets[firstDataIdx+1:]):
        if data and hasattr(data,"domain"):
            newTable, diffAttr = localConcatenation(newTable, data, datasets, useFirstAsLeader, mergeDomains)
            allDiffAttr += diffAttr
        else:
            if verbose >0: print "WARNING: The dataset of index",idx+firstDataIdx+1,"is not valid and was ignored."
    diffList = {}
    for attr in allDiffAttr:
        diffList[attr] = allDiffAttr.count(attr)

    return (newTable, diffList)




def loadSMI(filePath, domain = None, **args):
    """Loads an SMI space separated file with 2 columns and no header.
       Returns an ExampleTable with 2 attributes: SMILES (first column on SMI) and MOLNAME (second column of SMI)
       It ignores args and domain input variables, just for compatibility propose
    """
    if not os.path.isfile(filePath):
            if verbose >0: print "ERROR: Could not open file for reading: ", filePath
            return None
    file = open(filePath)
    lines = file.readlines()

    # Create New Domain (classless)
    domain = orange.Domain([orange.StringVariable("SMILES") , orange.StringVariable("MOLNAME")] , 0)
    # Create an empty example table with domain "domain"
    examples = DataTable(domain)

    # Fill the ExampleTable with the values from file
    for idx,line in enumerate(lines):
        sLine = line.split()
        if len(sLine) >= 2:
            examples.append([sLine[0],sLine[1]])
        elif len(sLine) == 1:
            examples.append([sLine[0],'?'])
    examples.attributeLoadStatus = [orange.Variable.MakeStatus.OK] * 2
    return examples

###############################################################################################################
#                  Set of functions to find/fix duplicated names in a specific domain                         #

def findDuplicatedNames(domain):
    """If the domain contains duplicated attribute names, it returns a dict with 
       duplicated names and the number of attributes for each.
       -You would probably want to use the function fixDuplicatedNames(domain) which will call this one!
    """
    domainNames = [attr.name for attr in domain]
    dup = {}
    for idx,attr in enumerate(domainNames):
        if domainNames.count(attr) > 1:
            if attr in dup:
                dup[attr].append(idx)
            else:
                dup[attr] = [idx]
    return dup

def fixDuplicatedNames(domain,reverse = False):
    """ fix the names of domain specified in the dupVectro (returned by findDuplicatedNames) by adding a sequencial 
        numbered suffix to the attributes with same name ex: Activity_1
        The reverse flag indicates that the rename will start from the class towards the first Attribute in the domain
        NOTE: This function DOES MODIFY the input domain. Use it carefully!!
              It tries to rename first the attributes and leave the class name unchanged
    """
    dupVector = findDuplicatedNames(domain)
    for attr in dupVector:
        dupVector[attr].sort(reverse = reverse)
        for idx,n in enumerate(dupVector[attr][:-1]):
            domain[n].name = attr + "_"+str(idx+1)
    return dupVector

#                                                                                                             #
###############################################################################################################


def getDataWithoutUnusedValues(data, checkClass = True, removeOneValued=False):
    """Create a new domain if, for the discrete attrobutes, there are unused values. The new
       domain attributes will not contain the unused values anymore.
       if the domain was changed, create a new datatable with the new domain and the old data.
       if nothing changes (there were no unused values) then return the same data
    """
    #create a new domain with new attributes although this method will add "R_" to the attribute name of changed attributes
    newDomain =orange.RemoveUnusedValues(data, 0, removeOneValued, checkClass)
    changedAttributes = []
    # only do anything if the resultant domain did in fact change
    if newDomain != data.domain:
         #On the new Domain rename the changed attributes back to their original names
         for newAttr in newDomain:
             if newAttr.name not in data.domain and newAttr.name[0:2]=="R_" and newAttr.name[2:] in data.domain:
                 newDomain[newAttr].name = newAttr.name[2:]
                 if verbose > 0: 
                     print "Changed values of attribute "+str(newAttr)
                     print "   Before: "+str(data.domain[newAttr.name].values)
                     print "   Now:    "+str(newAttr.values)
         data = DataTable(newDomain, data)
    return data

  
def getDomainWithoutDupSameAttr(domain):
    """ For the attributes in domain that are the same, remove them from domain leaving just one of them
        This returns a dict with attribute names removed and how many od them were actually removed and 
        the respective new Domain in a tuple: (<namesDict>, <newDomain>)
        NOTE: Two attribute are considered the same if orange loaded them as being the SAME i.e. they had
              the same name and represented the same variable, and so, a change to the variable domain would
              be applied to both. It doesn't mean that they had the same values!
    """
    metaIds = domain.getmetas().keys()
    metaIds.sort()
    classVar = domain.classVar
    metaAttrs = [domain[i] for i in metaIds]
    allAttrs = [attr for attr in domain] + [domain[i] for i in metaIds]
    sameDuplicates = []
    for idxA,attrA in enumerate(allAttrs):
        for idxB in range(idxA+1,len(allAttrs)):
            if (allAttrs[idxA] is allAttrs[idxB]):
                if idxB not in sameDuplicates:
                    sameDuplicates.append(idxB)
    ##sameDuplicates.sort(reverse = True)
    RemovedVars = [allAttrs[var].name for var in sameDuplicates]
    if not RemovedVars:
        return ({},domain)
    allAttrs = [attr for idx,attr in enumerate(allAttrs) if idx not in sameDuplicates]
    newDomainAttr = [attr for attr in allAttrs if attr not in metaAttrs]
    newDomainMetas = [attr for attr in metaAttrs if attr in allAttrs]
    if classVar and classVar in allAttrs:
        newDomainAttr.remove(classVar)
        classVar = [classVar]
        newDomain = orange.Domain(newDomainAttr + classVar)
    else:
        classVar=[]
        newDomain = orange.Domain(newDomainAttr,0)

    if newDomainMetas:
        for idx,attr in enumerate(newDomainMetas):
            newDomain.addmeta(-(idx-1),attr)
    rmList = {}
    for attr in RemovedVars:
        rmList[attr] = RemovedVars.count(attr)
    return (rmList, newDomain)
 

def getCopyWithoutMeta(data):
    """
    Complete removal of all meta information from an Example or ExampleTable object 
    This method does not modify the data/example passed in the argument but it has to duplicate the dataset
    """
    attributes=data.domain[0:len(data.domain)]
    attrDomain=orange.Domain(attributes)
    if type(data)==orange.Example:
        newdata = orange.Example(attrDomain, data)
    else:
        newdata = DataTable(attrDomain, data)
    return newdata


def ConvertEx2Domain(domain,example,examplesFixedLog=None,fixClass = False):
    print "The method ConvertEx2Domain is no longer available. Please use the class ExFix"
    return None


class ExFix(object):
    """
    Class for fixing the examples to a particular domain
    """
    def __init__(self,domain = None,examplesFixedLog=None,fixClass = False):
            self.domain = domain    # Domain to match each example to be fixed
            self.exDomain = None    # Domain of the example used to test the compatibility. All exampleas will have to mach this, or else 
                                    #    new compatibility tests need to be performed
            self.examplesFixedLog = examplesFixedLog
            self.fixClass = fixClass
            self.attrNamesIdx = {}
            self.fixFuns = {}#  Ex:  {"Attr1":labmda}    #Varibles that needs to be fixed for the particular exDomain
            self.ready = False
            #Vars for setting the increment counter of fixed examples
            self.fixOrder = False
            self.fixTypes = False
            self.fixCount = False
        
    def __cvtDiffDisc(self,attr):   # Just used for discrete attributes
            if attr.value in self.domain[attr.variable.name].values:
                return attr.value
            else:
                self.__composeVarLog(attr, isUnknown=True)
                self.fixTypes=True
                return '?'

    def __cvt2Disc(self,attr):  # Converts any attr disc or cont to discrete
            value = str(attr.value)
            self.fixTypes=True
            if value in self.domain[attr.variable.name].values:
                self.__composeVarLog(attr)
                return value
            else:
                self.__composeVarLog(attr,isUnknown=True)
                return '?'
    def __cvt2Str(self,attr):  # Converts any attr to String
            self.fixTypes=True
            self.__composeVarLog(attr)
            return str(attr.value)

    def __cvt2Cont(self,attr):  # Converts any attr to continuous
            self.fixTypes=True
            if miscUtilities.isNumber(attr.value):
                self.__composeVarLog(attr)
                return float(attr.value)
            else:
                self.__composeVarLog(attr, isUnknown=True)

                return '?'

    # The lambda function is used instead this one
    #def __cvt2Ukn(self,attr): 
    #        return '?'

    # The lambda function is used instead this one
    #def __cvtFullCompatible(self,attr):
    #        return attr.value
        
    def set_domain(self, domain = None):
        if not domain or not self.domain or domain is not self.domain:
            self.ready = False
            self.domain = domain

    def set_examplesFixedLog(self, examplesFixedLog = None):
        self.examplesFixedLog = examplesFixedLog

    def fixExample(self,example):
            if self.ready and self.exDomain is example.domain: 
                if not self.fixFuns:
                    return example
                else:
                    newEx=orange.Example(self.domain) #newEx=[?,?,?,...,?]
                    self.fixTypes = False
                    for idx,attr in enumerate(self.domain.attributes):
                        # Here is the bottle-Neck: Reading or writing to ExampleTable using "ATTRName" 
                        #    is aboout 4x slower that using the attribute itself or the index
                        newEx[idx] = self.fixFuns[attr.name](example[self.attrNamesIdx[attr.name]])
                    if self.fixClass and newEx.domain.classVar:
                        newEx[newEx.domain.classVar] = self.fixFuns[newEx.domain.classVar.name](example[self.attrNamesIdx[newEx.domain.classVar.name]])
                    if self.fixTypes:
                        self.__composeNExFixedLog()
                    return newEx
            # This is only reached the first time the fixExample is called in order to fill the fixFuns
            #   OR if the domain needs new evaluation
            self.ready = False
            if not self.domain or not example:
                return None
            self.exDomain = example.domain
            # Using this dictionary with the correspondence of attr names and idx makes the fix much faster
            self.attrNamesIdx = {}
            for idx,attr in enumerate(self.exDomain):
                self.attrNamesIdx[attr.name] = idx
            if self.domain is self.exDomain:
                self.ready = True
                self.fixFuns = {}
                return example
            if self.fill_fixFuns():
                self.ready = True
                return self.fixExample(example)
            else:
                return None
            
    def __composeVarLog(self, attr, isUnknown=False):
            # Only compose the log if a handler was provided in self.examplesFixedLog
            if (self.examplesFixedLog!=None) and (type(self.examplesFixedLog) == types.DictType):
                varBefore = str(self.exDomain[attr.variable.name]).replace("'" + str(attr.variable.name)  + "'" , "").strip()
                varAfter = str(self.domain[attr.variable.name]).replace("'" + str(attr.variable.name)  + "'","").strip()
                if varBefore==varAfter and self.domain[attr.variable.name].varType == orange.VarTypes.Discrete:
                    changeTxt = "EnumVar - Attribute Values were different"
                else:
                    changeTxt = varBefore +" to " + varAfter
                if ("Vars needing type fix" in self.examplesFixedLog):
                    if attr.variable.name not in self.examplesFixedLog["Vars needing type fix"]:
                        self.examplesFixedLog["Vars needing type fix"][attr.variable.name] = changeTxt
                else:
                        self.examplesFixedLog["Vars needing type fix"] = {attr.variable.name : changeTxt}
                if isUnknown and \
                "(some impossible conversions" not in self.examplesFixedLog["Vars needing type fix"][attr.variable.name]:
                    self.examplesFixedLog["Vars needing type fix"][attr.variable.name] += " (some impossible conversions. It was set to '?' for some examples.)"
                    if verbose >0: print "The attr '"+attr.variable.name+"' ["+str(self.exDomain.index(attr.variable.name))+"] of the following example was set to '?'.\n   "

    def __composeNExFixedLog(self):
            if (self.examplesFixedLog!=None) and (type(self.examplesFixedLog) == types.DictType):
                #Compose the log that will be used for warning messages
                if self.fixTypes:
                    if "Fixed Types of variables" in self.examplesFixedLog:
                        self.examplesFixedLog["Fixed Types of variables"]+=1
                    else:
                        self.examplesFixedLog["Fixed Types of variables"]=1
                        if verbose >0: print "WARNING: Fixed Types of variables"
                # For Now just using the fixType
                #if self.fixCount:
                #    if "Fixed Number of variables" in self.examplesFixedLog:
                #        self.examplesFixedLog["Fixed Number of variables"]+=1
                #    else:
                #        self.examplesFixedLog["Fixed Number of variables"]=1
                #        if verbose >0: print "WARNING: Fixed Number of variables"
                #if self.fixOrder:
                #    if "Fixed Order of variables" in self.examplesFixedLog:
                #        self.examplesFixedLog["Fixed Order of variables"]+=1
                #    else:
                #        self.examplesFixedLog["Fixed Order of variables"]=1
                #        if verbose >0: print "WARNING: Fixed Order of variables"


    def fill_fixFuns(self):
            """
            returns False if something is missing, or not compatible at all
            returns True if the self.fixFuns was set properly
            """
            # newEx=orange.Example(self.domain) #newEx=[?,?,?,...,?]
            for attr in self.domain:
                if not self.fixClass and self.domain.classVar and (attr == self.domain.classVar): # In this case the classVar will be set to '?'
                    self.fixFuns[attr.name] = lambda x:'?' # self.__cvt2Ukn   
                    continue  # In this case the classVar will be set to '?' by default
                elif attr.name not in self.exDomain:   # When there is at least one attribute missing, they are imcompatible!
                    if self.examplesFixedLog != None:
                        if "Missing Attributes" in self.examplesFixedLog:
                            if attr.name in self.examplesFixedLog["Missing Attributes"]:
                                self.examplesFixedLog["Missing Attributes"][attr.name] += 1
                            else:
                                self.examplesFixedLog["Missing Attributes"][attr.name] = 1
                                if verbose >0: print "WARNING: The attribute " + attr.name +" is missing in the example"
                        else:
                            self.examplesFixedLog["Missing Attributes"] = {attr.name:1}
                            if verbose >0: print "WARNING: The attribute " + attr.name +" is missing in the example"
                    return False #if there was already one attribuite missing, don't bother fixing anythin!
                else: 
                    if attr is self.exDomain[attr.name]:
                        self.fixFuns[attr.name] = lambda x:x.value #self.__cvtFullCompatible
                    elif attr.varType == self.exDomain[attr.name].varType:
                        if attr.varType == orange.VarTypes.Discrete:
                            ## it depends on var value   using __cvtDiffDisc()
                            self.fixFuns[attr.name] = self.__cvtDiffDisc
                        else:
                            self.fixFuns[attr.name] =  lambda x:x.value# self.__cvtFullCompatible
                    elif self.domain[attr].varType == orange.VarTypes.Discrete:
                        ## it depends on var value   using  __cvt2Disc()
                        self.fixFuns[attr.name] = self.__cvt2Disc
                    elif self.domain[attr].varType == orange.VarTypes.String:
                        self.fixFuns[attr.name] = self.__cvt2Str
                    else:  # Thr domain[attr].varType == orange.VarTypes.Continuous
                        ## it depends on var value   using  __cvt2Cont()
                        self.fixFuns[attr.name] = self.__cvt2Cont
            return True
        ##ecPA

def rmAllMeta(data):
    """
    Complete removal of all meta information from an example or ExampleTable object 
    The argument passed to this method is modified, and nothing is returned
    """
    metaAttr = data.domain.getmetas()
    for meta in metaAttr:
##scPA
        if type(data)==orange.Example:
            data.removemeta(meta)
        else:
##ecPA
            data.removeMetaAttribute(meta)
        data.domain.removemeta(meta)

def list2string(inList):
    """Converts a list of numbers to a list of char numbers. """
    stringList = []
    for elem in inList:
        stringList.append(str(elem))
    return stringList


def getScaledVar(attr, domstat, factor = 2):
    """Define a transformation for an attribute which subtracts the average and 
       scales with factor*std. """
    attr_c = orange.FloatVariable(attr.name)
    attr_c.getValueFrom = orange.ClassifierFromVar(whichVar = attr)
    transformer = orange.NormalizeContinuous()
    attr_c.getValueFrom.transformer = transformer
    transformer.average = domstat[attr].avg
    transformer.span = domstat[attr].dev*factor
    return attr_c


def getScaledVarStat(attr, statData, factor = 2):
    """Define a transformation for an attribute which subtracts the average and 
       scales with factor*std. statData should be an exampleTable with the mean in the first row and the
       std in the second. """
    attr_c = orange.FloatVariable(attr.name)
    attr_c.getValueFrom = orange.ClassifierFromVar(whichVar = attr)
    transformer = orange.NormalizeContinuous()
    attr_c.getValueFrom.transformer = transformer
    transformer.average = statData[0][attr]
    transformer.span = statData[1][attr]*factor
    return attr_c


def scaleContVarStat(data, statData, factor = 2):
    """Scales all continous attributes (not the pred/class variable) in a data set.
       The greater the factor, the smaller the output values. 
       The data is scaled with the mean and std of the second data set statData."""

    # Partition the attributes into discrete and continious attributes
    discAttrList = []
    contAttrList = []
    for attribute in data.domain.attributes:
        if attribute.varType != orange.VarTypes.Discrete:
            contAttrList.append(attribute)
        else: discAttrList.append(attribute)

    # Get the continious part of the data set
    newContDomain = orange.Domain(contAttrList, data.domain.classVar)
    contData = DataTable(newContDomain, data)

    # Scale the continious attributes
    #domstat = orange.DomainBasicAttrStat(contData)
    newattrs = []
    for attr in contData.domain.attributes:
        newAttr = getScaledVarStat(attr, statData, factor)
        newattrs.append(newAttr)

    # Create a new data set by merging the scaled continous variables and the discrete ones
    newDomain = orange.Domain(newattrs + discAttrList, data.domain.classVar)
    scaledData = DataTable(newDomain, data)

    return scaledData


def scaleContVar(data, factor = 2):
    """Scales all continous attributes (not the pred/class variable) in a data set.
       The greater the factor, the smaller the output values. 
       The data is scaled with the mean and std of the data set itself."""

    # Partition the attributes into discrete and continious attributes
    discAttrList = []
    contAttrList = []
    for attribute in data.domain.attributes:
        if attribute.varType != orange.VarTypes.Discrete:
            contAttrList.append(attribute)
        else: discAttrList.append(attribute)

    # Get the continious part of the data set
    newContDomain = orange.Domain(contAttrList, data.domain.classVar)
    contData = DataTable(newContDomain, data)

    # Scale the continious attributes
    domstat = orange.DomainBasicAttrStat(contData)
    newattrs = []
    for attr in contData.domain.attributes:
        newAttr = getScaledVar(attr, domstat, factor)
        newattrs.append(newAttr)

    # Create a new data set by merging the scaled continous variables and the discrete ones
    newDomain = orange.Domain(newattrs + discAttrList, data.domain.classVar)
    scaledData = DataTable(newDomain, data)

    return scaledData

def horizontalMerge(dataAin, dataBin, varAin, varBin):
    """
    dataAin and dataBin are orange ExampleTable objects and varA and varB are strings.
    Adding the attributes of dataBin to the dataAin data set requiering identical values of
    the varA and varB attributes for merged examples.
 
    The merged dataset will always contain all examples from dataAin which attributes will be 
    the ones of dataAin extended (added) with the ones from dataBin of the example where the varBin
    had the same value as the value of varAin. 
    For the examples where there is no match for the value of varAin in the varBin of dataBin, the 
    extended attributes will be set to "?"

    Class variable of merged dataset:
        -The class of merged data will be the dataAin or dataBin class or Classless, tested in this order.
    NOTE: The class (if exists) of the dataset which was not used as class of merged Data, will appear in the
          merged Data as an attribute!
          Meta Attributes will also appear as regular attributes

    When the resulting merged dataset would have attributes with same name (happens often with the class var),
    they are renamed (the oned from dataB) by adding a numbered suffix Ex: Activity  ->  Activity_1     

    DataAin is kind of 'the Leader' and should be the one containig the class and the dataBin sould be classless 
    """
    def prepareDataSets(data1,data2,attr1Name,attr2Name):
      """This will return 2 datasets in a tuple (dataA,dataB) assuring that both varA and varB are fully compatible and
       is safe to proceed for the merge action"""
      attr1=data1.domain[attr1Name]
      attr2=data2.domain[attr2Name]
      if (attr1.varType == attr2.varType == orange.VarTypes.Continuous) or (attr1 == attr2):
        #Nothing to be done
        dataA = data1
        dataB = data2
      else:
        if attr1.varType==orange.VarTypes.Discrete and attr2.varType==orange.VarTypes.Discrete:
                #may have different values 
                allValues = attr1.values.native()
                for value in attr2.values.native():
                    if value not in allValues:
                        allValues.append(value)
                NewVar1 = orange.EnumVariable(attr1.name, values = allValues)
                NewVar2 = orange.EnumVariable(attr2.name, values = allValues)
        elif attr1.varType == orange.VarTypes.String or attr2.varType==orange.VarTypes.String:
                # They are not compatible, and we need to use string type since at least one of it it's this type.
                NewVar1 = orange.StringVariable(name=attr1.name)
                NewVar2 = orange.StringVariable(name=attr2.name)
        else:
                # They are not compatible, need to decide what type to take
                totalExamples=0
                OKCvtToCont = 0
                if attr1.varType==orange.VarTypes.Discrete:
                    dataD=data1
                    dataC=data2
                    attrD = attr1
                    attrC = attr2
                else:
                    dataC=data1
                    dataD=data2
                    attrD = attr2
                    attrC = attr1
                for ex in dataD:
                        totalExamples += 1
                        if ex[attrD].value in ['?','~'] or miscUtilities.isNumber(ex[attrD].value):
                            OKCvtToCont += 1
                #If >=50% can be securely converted to Continuous, use the continuous attribute
                if OKCvtToCont >= (totalExamples/2):
                    NewVar1 = orange.FloatVariable(attr1.name)
                    NewVar2 = orange.FloatVariable(attr2.name)
                else:
                    allValues = attrD.values.native()
                    for ex in dataC:
                        if str(ex[attrC].value) not in allValues:
                            allValues.append(str(ex[attrC].value))

                    NewVar1 = orange.EnumVariable(attr1.name, values = allValues)
                    NewVar2 = orange.EnumVariable(attr2.name, values = allValues)

        attrsA = [attr for attr in data1.domain]
        attrAidx = data1.domain.index(NewVar1.name)
        attrsA[attrAidx] = NewVar1

        attrsB = [attr for attr in data2.domain]
        attrBidx = data2.domain.index(NewVar2.name)
        attrsB[attrBidx] = NewVar2

        newDomainA = orange.Domain(attrsA,0)
        newDomainA.addmetas(data1.domain.getmetas())
        newDomainB = orange.Domain(attrsB,0)
        newDomainB.addmetas(data2.domain.getmetas())
        dataA = DataTable(newDomainA)
        dataB = DataTable(newDomainB)

        fixer = ExFix(dataA.domain,fixClass = True)
        for ex in data1:
            dataA.append(fixer.fixExample(ex))
            for meta in data1.domain.getmetas().values(): 
                dataA[-1][meta.name] = ex[meta.name]
        fixer.set_domain(dataB.domain)
        for ex in data2:
            dataB.append(fixer.fixExample(ex))
            for meta in data2.domain.getmetas().values():
                dataB[-1][meta.name] = ex[meta.name]

      return (dataA,dataB)


    if not dataAin or not dataBin or not varAin or not varBin:
        return None
    if varAin not in dataAin.domain or varBin not in dataBin.domain:
        return None

    #Decide What var will be the class
    if dataAin.domain.classVar:
        classVarName = dataAin.domain.classVar.name
        classVarIsA = True
    elif dataBin.domain.classVar:
        classVarName = dataBin.domain.classVar.name
        classVarIsA = False
    else:
        classVarName = None

    #Convert all meta-attributes to regular attributes
    classA = dataAin.domain.classVar
    attrs = [attr for attr in dataAin.domain if attr is not classA] + dataAin.domain.getmetas().values()
    dataA_CvtMeta = DataTable(orange.Domain(attrs,classA),dataAin)
    classB = dataBin.domain.classVar
    attrs = [attr for attr in dataBin.domain if attr is not classB] + dataBin.domain.getmetas().values()
    dataB_CvtMeta = DataTable(orange.Domain(attrs,classB),dataBin) 

    dataSets = prepareDataSets(dataA_CvtMeta,dataB_CvtMeta,varAin,varBin)    

    varListA = []
    varListB = []
    # dataA
    dc = orange.Domain([x.clone() for x in dataSets[0].domain], None)
    for i, a in enumerate(dc):
        a.getValueFrom = lambda ex,f,i=i: ex[i]
    dc.addmetas(dataSets[0].domain.getmetas())
    dataA = DataTable(dc, dataSets[0])
    varListA = dataA.domain.variables.native() + dataA.domain.getmetas().values()    
    # dataB
    dc = orange.Domain([x.clone() for x in dataSets[1].domain], None)
    for i, a in enumerate(dc):
        a.getValueFrom = lambda ex,f,i=i: ex[i]
    dc.addmetas(dataSets[1].domain.getmetas())
    dataB = DataTable(dc, dataSets[1])
    varListB = dataB.domain.variables.native() + dataB.domain.getmetas().values()               
    
    # chosen vars
    varA = [var for var in varListA if var.name == varAin][0]
    varB = [var for var in varListB if var.name == varBin][0]
    #Merge the dataA and dataB
    val2idxDictA = {}
    for eIdx, e in enumerate(dataA):
        val2idxDictA[e[varA].native()] = eIdx
    val2idxDictB = {}
    for eIdx, e in enumerate(dataB):
        val2idxDictB[e[varB].native()] = eIdx
    if val2idxDictA.has_key("?"):
        val2idxDictA.pop("?")
    if val2idxDictA.has_key("~"):
        val2idxDictA.pop("~")
    if val2idxDictA.has_key(""):
        val2idxDictA.pop("")
    if val2idxDictB.has_key("?"):
        val2idxDictB.pop("?")
    if val2idxDictB.has_key("~"):
        val2idxDictB.pop("~")
    if val2idxDictB.has_key(""):
        val2idxDictB.pop("")
    # example table names
    nameA = dataA.name
    if not nameA: nameA = "Examples A"
    nameB = dataB.name
    if not nameB: nameB = "Examples B"
    # create example B with all values unknown
    exBDK = orange.Example(dataB[0])
    for var in varListB:
        exBDK[var] = "?"
    # build example table to append to the right of A
    vlBreduced = list(varListB)
    vlBreduced.remove(varB)
    domBreduced = orange.Domain(vlBreduced, None)
    etBreduced = DataTable(domBreduced)
    for e in dataA:
        dataBidx = val2idxDictB.get(e[varA].native(), None)
        if dataBidx <> None:
            etBreduced.append(dataB[dataBidx])
        else:
            etBreduced.append(orange.Example(domBreduced, exBDK))
     
    # Rename the classVar to use so that it does not get renamed by fixDuplicatedNames
    ClassTag = "@CLASS@"
    if classVarName:
        if not classVarIsA and classVarName in etBreduced.domain:
            etBreduced.domain[classVarName].name = classVarName + ClassTag
        elif not classVarIsA and classVarName in dataA.domain:
            dataA.domain[classVarName].name = classVarName + ClassTag
        elif classVarIsA and classVarName in dataA.domain:
            dataA.domain[classVarName].name = classVarName + ClassTag
        elif classVarIsA and classVarName in etBreduced.domain:
            etBreduced.domain[classVarName].name = classVarName + ClassTag
        else:
            raise Exception("ERROR: Problems finding the class var in horizontalMerge.")
      

    #Using the ExampleTable, duplicated names will be allowed
    etAB = orange.ExampleTable([dataA, etBreduced]) 
    # Fix any existing duplicated names in the output dataset 
    # Here we want to rename the dataB variables and leave the dataA variables with original names 
    # The classVar will not be renamed since it is already tagged. 
    
    fixDuplicatedNames(etAB.domain, reverse = True)

    # Set the correct class var
    varList = etAB.domain.variables.native() + etAB.domain.getmetas().values()
    if classVarName:
        classVar = etAB.domain[classVarName+ClassTag]
        outDomain = orange.Domain([attr for attr in varList if attr.name != classVarName+ClassTag],classVar)
    else:
        outDomain = orange.Domain(varList,None)
    etAB = orange.ExampleTable(outDomain,etAB)
    if classVarName:
        # Rename the vars with same name as future class to XXXX_N
        idx = 1
        for attr in etAB.domain:
            if attr.name == classVarName:
                attr.name = classVarName+"_"+str(idx)
                idx += 1
        # Rename back the class Var
        etAB.domain[classVarName+ClassTag].name = classVarName
    etAB.name = "%(nameA)s (merged with %(nameB)s)" % vars()
    return etAB



def attributeDeselectionData(data, attributeList):
    """Deselects the attributes passed in attributeList from data"""
    newDomainList = []
    for idx in range(len(data.domain.attributes)):
        if string.strip(data.domain.attributes[idx].name) not in attributeList:
            newDomainList.append(data.domain[idx])
    newDomain = orange.Domain(newDomainList, data.domain.classVar)
    #Handle the Meta Attributes
    for attrID in data.domain.getmetas():
        if string.strip(data.domain[attrID].name) not in attributeList:
                newDomain.addmeta(attrID,data.domain[attrID])
    scaledData = data.select(newDomain)
    return scaledData

def attributeSelectionData(data, attributeList):
    """The returned data set will have the attributes in the same order as attribute list."""

    newDomainList = []
    for idx in range(len(data.domain.attributes)):
        if string.strip(data.domain.attributes[idx].name) in attributeList:
                newDomainList.append(data.domain[idx])
    newDomain = orange.Domain(newDomainList, data.domain.classVar)
    #Handle the Meta Attributes
    for attrID in data.domain.getmetas():
        if string.strip(data.domain[attrID].name) in attributeList:
            newDomain.addmeta(attrID,data.domain[attrID])

    scaledData = data.select(newDomain)
    return scaledData

def attributeDeselectionExample(example, attributeList):

    newDomainList = []
    for idx in range(len(example.domain.attributes)):
        if example.domain.attributes[idx].name not in attributeList:
            newDomainList.append(example.domain[idx])
    newDomain = orange.Domain(newDomainList, example.domain.classVar)
    newExample = orange.Example(newDomain, example)
    return newExample


def addDummyClassExample(example, classVariable):

    newDomain = orange.Domain(example.domain, classVariable)
    newExample = orange.Example(newDomain, example)
    return newExample
    

def getNewDomain(example):
    """Takes an orange example instance and returns the corresponding domain modified
       by changing all discrete attributes to numbers.  """
    domainList = []
    for attribute in example.domain.attributes:
        if attribute.varType == orange.VarTypes.Discrete:
            valueList = list2string(range(len(attribute.values)))
            valueDict = dict(zip(attribute.values,valueList))
            domainList.append(orange.EnumVariable(attribute.name, values = valueList))
        else:
            domainList.append(attribute)
    # Append the prediction variable without transformin the values. Prediction values are taken care of seperatly
    # because they need to be back transformed. 
    try: 
        if example.domain.classVar: domainList.append(example.domain.classVar)
    except: pass
    newDomain = orange.Domain(domainList)
    return newDomain

def getNewExample(example, newDomain):
    """ Create the new exampleTable data set where the values of the enumVar are changed in
        accordance with newDomain"""
    exampleList = []
    for attribute in example.domain.attributes:
        if attribute.varType == orange.VarTypes.Discrete:
            valueList = list2string(range(len(attribute.values)))
            valueDict = dict(zip(attribute.values,valueList))
            exampleList.append(valueDict[str(example[attribute])])
        else:
            exampleList.append(example[attribute])
    try: exampleList.append(example[example.domain.classVar])
    except: pass
    newExample = orange.Example(newDomain, exampleList)
    return newExample


def changeEnumVarData(data):
    """Takes an exampleTable instance and converts all its enumVariables (orange.VarType.Discrete)
       to consecutive numbers. This is useful for training methods which can not handle
       character values. Removes any meta data from the returned  exampleTable instance."""
    # Define the new domain. The values of all enumVar are set to an ordered number list.
    # Keep outside of the loop over examples or the returned exampleTable will be
    # confused by differing domain definitions.
    example = data[0]
    newDomain = getNewDomain(example)

    # Create the new exampleTable data set where the values of the enumVar are changed in
    # accordance with newDomain
    scaledData = DataTable(newDomain)
    for example in data:
        newExample = getNewExample(example, newDomain)
        scaledData.append(newExample)
    return scaledData


def changeEnumVarExample(example):
    """Takes an orange example instance and converts all its enumVariables (orange.VarType.Discrete)
       to consecutive numbers. This is useful for training methods which can not handle
       character values. Removes any meta data from the returned example instance."""
    # Define the new domain. The values of all enumVar are set to an ordered number list.
    newDomain = getNewDomain(example)

    newExample = getNewExample(example, newDomain)
    return newExample



##scPA
class scalizer(object):
    def __new__(cls,**kwds):
        self = object.__new__(cls)
        self.__dict__.update(kwds)
        return self

    def __init__(self, **kwds):
        """
        data - Data to be scaled
        nMin - attributes scaling lower limit
        nMax - attributes scaling upper limit
        nClassMin - Class scaling lower limit
        nClassMax - Class scaling upper limit
            if classMin or classMax is None, no scaling will be performed in class 
        scaledData - variable to place the entire scaling data. If is None, only scaling parameters will be calculated
        Returns scaling parameters
            if error occurred, None will be returned
        """
        self.callback = None
        self.numberOfDecimals = 7
        self.verbose = 0
        self.nMin = -1.0
        self.nMax = 1.0
        self.scaleClass = False
        self.nClassMin = None
        self.nClassMax = None
        self.dataAnalyzed = False
        self.file = None
        self.data = None
        self.domain = None
        self.dataDomain = None

        self.varNames = []
        self.varNameIdx = {}
        self.minimums = []
        self.maximums = []
        self.values = []
        # Append arguments to the __dict__ member variable
        self.__dict__.update(kwds)
        if self.file != None:
            if not self.loadScalingValues(self.file):
                return
        else:
            if self.data != None:
                self.analyzeData(self.data)


    def loadScalingValues(self,file):
        """returns True if success, or False if errors found"""
        self.dataAnalyzed = False
        if not os.path.isdir(file):
            if self.verbose >0: print "ERROR: Could not find data files"
            return False
        domainFile = os.path.join(file,"scalizerDomain.tab")
        dataDomainFile = os.path.join(file,"scalizerDataDomain.tab")
        valuesFile = os.path.join(file,"scalizerValues.txt")

        if not os.path.isfile(domainFile) or not os.path.isfile(dataDomainFile) or not os.path.isfile(valuesFile) :
            if self.verbose >0: print "ERROR: Could not open files for reading"
            return False
        try:
            scalingDataF = open(valuesFile,"r")
            lines = scalingDataF.readlines()
            scalingDataF.close 
            self.varNames = eval(lines[0][:-2])
            self.minimums = eval(lines[1][:-2])
            self.maximums = eval(lines[2][:-2])
            self.values = eval(lines[3][:-2])
            otherParams = eval(lines[4][:-2])
            self.nMin = otherParams[0]
            self.nMax = otherParams[1]
            self.scaleClass = otherParams[2]
            self.nClassMin = otherParams[3]
            self.nClassMax = otherParams[4]
            self.domain = DataTable(domainFile,createNewOn=orange.Variable.MakeStatus.OK).domain
            self.dataDomain = DataTable(dataDomainFile,createNewOn=orange.Variable.MakeStatus.OK).domain
            self.varNameIdx={}
            for idx,attr in enumerate(self.domain):
                self.varNameIdx[attr.name] = idx
                attr.numberOfDecimals = self.numberOfDecimals

        except:
            if self.verbose >0: print "ERROR: Could not load the scaling data present on ",file
            return False
        
        #Check loaded variables
        if not len(self.varNames)==len(self.minimums)==len(self.maximums)==len(self.values):
            if self.verbose >0:print "ERROR: scaling data on File was not with correct lenght"
            return False
        self.dataAnalyzed = True
        return True

    def convertClass(self, value):
        """ Converts the class returning the value if continuous-class or the order number of the discrete value !
            The return value is a number
            The input argument 'value' must be a number""" 
        if not self.scaleClass:
            return value
        if not self.dataAnalyzed:
            if self.verbose >0: print  "ERROR: Data was not analyzed yet"
            return False
        else:
            idxVar=self.varNameIdx[self.dataDomain.classVar.name]
      
        res = self.minimums[idxVar]+((value-self.nClassMin)/(self.nClassMax-self.nClassMin))*\
                                        (self.maximums[idxVar]-self.minimums[idxVar])
        return res
        


    def saveScalingValues(self,file):
        """returns True if success, or False if errors found"""
        if not self.dataAnalyzed:
            if self.verbose >0: print  "ERROR: Data was not analyzed yet"
            return False
        try:
            if os.path.isdir(file):
                os.system("rm -rf "+file)
            os.system("mkdir -p "+file)
            DataTable(self.domain).save(os.path.join(file,"scalizerDomain.tab"))
            DataTable(self.dataDomain).save(os.path.join(file,"scalizerDataDomain.tab"))
            scalingDataF = open(os.path.join(file,"scalizerValues.txt"),"w")
        except:
            if self.verbose >0: print "ERROR: Could not open file for writing"
            return False
        scalingDataF.write(str(self.varNames)+"\r\n")        
        scalingDataF.write(str(self.minimums)+"\r\n")        
        scalingDataF.write(str(self.maximums)+"\r\n")        
        scalingDataF.write(str(self.values)+"\r\n")        
        scalingDataF.write(str([self.nMin,self.nMax,self.scaleClass,self.nClassMin, self.nClassMax])+"\r\n")        
        scalingDataF.close()
        return True

    def analyzeData(self,data):
        """returns True if success, or False if errors found"""
        if not data:
            if self.verbose >0: print "No data to Analyze"
            return False
        if data.hasMissingValues():
            if self.verbose >0: print "The data contains missing values, please impute the data first"
            return False
        self.dataAnalyzed = False
        self.domain = None
        self.dataDomain = None
        self.varNameIdx = {}
        dist = orange.DomainDistributions(data)
        for idx,attr in enumerate(data.domain):
            self.varNames.append(attr.name)
            self.varNameIdx[attr.name] = idx
            if attr.varType == orange.VarTypes.Discrete:
                values = []
                for idx in range(len(data.domain[attr].values)):
                    values.append(str(data.domain[attr](idx).value))
                self.minimums.append(0.0)
                self.maximums.append(float(len(values)-1))
                self.values.append(values)
            else:
                self.minimums.append(float(min(dist[attr].keys())))
                self.maximums.append(float(max(dist[attr].keys())))
                self.values.append([])

        # create empty dataset with all attributes being of type float
        if not data.domain.classVar:
            self.domain = orange.Domain([orange.FloatVariable('%s'%x) for x in self.varNames],False)
        else:
            self.domain = orange.Domain([orange.FloatVariable('%s'%x) for x in self.varNames])
        for attr in self.domain:
            attr.numberOfDecimals = self.numberOfDecimals
        self.dataDomain = data.domain
    
        self.dataAnalyzed = True
        return True

    def getScalingValues(self):
        if not self.dataAnalyzed:
            if self.verbose >0: print  "Warning: Data was not analyzed yet"
            return {}
        scalingValues = {}
        scalingValues["Var Names"] = self.varNames
        scalingValues["Minimums"] = self.minimums
        scalingValues["Maximums"] = self.maximums
        scalingValues["Discrete values"] = self.values
        scalingValues["Scaled vars Min"] = self.nMin
        scalingValues["Scaled vars Max"] = self.nMax
        scalingValues["Scale Class"] = self.scaleClass
        scalingValues["Scaled Class Min"] = self.nClassMin
        scalingValues["Scaled Class Max"] = self.nClassMax
        return scalingValues

    def scaleEx(self,ex, skipCheck=False):
        """returns an example with the scalizer domain (all vars continuous including the class) and all vars scaled  
           or False if errors found"""
        if skipCheck:
            return self.fastscaleEx(ex)
        #Make sure these are floats
        self.nMin = float(self.nMin)
        self.nMax = float(self.nMax)
        if not self.dataAnalyzed:
            if self.verbose >0: print  "ERROR: Data was not analyzed yet"
            return False

        if self.nMax<=self.nMin:
            if self.verbose >0: print "Attributes scaling limits are not correct!"
            return False

        if self.scaleClass:
            if self.nClassMin == None or self.nClassMax == None:
                if self.verbose >0: print "Warning: Class scaling limits were not defined. The same attributes limits will be used"
                self.nClassMax = self.nMax
                self.nClassMin = self.nMin
            elif self.nClassMax <= self.nClassMin:
                if self.verbose >0: print "Class scaling limits are not correct!"
                return False
            self.nClassMin = float(self.nClassMin)
            self.nClassMax = float(self.nClassMax)

        # create empty example [?,?,?,...,?] of the correct domain
        scaledEx = orange.Example(self.domain)
        
        #Scale Attributes  -  Scale function is according to libSVM code
        for attr in ex.domain.attributes:
            if not self.varNameIdx.has_key(attr.name):
                if self.verbose >0: print "Attribute ",attr," was not found in varNames local variable."
                return False
            else:
                idx = self.varNameIdx[attr.name]
            if attr.varType == orange.VarTypes.Continuous:
                if self.maximums[idx] == self.minimums[idx]:   #in case of single-valued attribute
                    scaledEx[idx]=(self.nMax + self.nMin)/2
                else:
                    scaledEx[idx]=self.nMin+(self.nMax-self.nMin)*\
                            (float(ex[attr].value)-self.minimums[idx])/\
                            (self.maximums[idx]-self.minimums[idx])
            else:  # is Discrete
                value_idx = int(ex[attr])
                if self.maximums[idx] == self.minimums[idx]:    #in case of single-valued attribute
                    scaledEx[idx]=(self.nMax + self.nMin)/2
                else:
                    scaledEx[idx]=self.nMin+(self.nMax-self.nMin)*\
                            (value_idx-self.minimums[idx])/\
                            (self.maximums[idx]-self.minimums[idx])
        #Scale class
        classVar = ex.domain.classVar
        if classVar!=None:
            if ex[classVar].value=="?":
                scaledEx[classVar.name]="?"
            else:
                if not self.varNameIdx.has_key(classVar.name):
                    if self.verbose >0: print "Class ",classVar," was not found in varNames local variable!!"
                    return False
                else:
                    idx = self.varNameIdx[classVar.name]
                if classVar.varType == orange.VarTypes.Continuous: 
                    if self.scaleClass:
                        if self.maximums[idx] == self.minimums[idx]:    #in case of single-valued attribute
                            scaledEx[idx]=(self.nClassMax + self.nClassMin)/2
                        else:
                            scaledEx[idx]=self.nClassMin+(self.nClassMax-self.nClassMin)*\
                                (float(ex[classVar].value)-self.minimums[idx])/\
                                (self.maximums[idx]-self.minimums[idx])
                    else: #Cont. but we dont want to scale it
                        scaledEx[idx] = ex[classVar].value
                else:  # is Discrete
                    if self.scaleClass:
                        value_idx = int(ex[classVar])
                        if self.maximums[idx] == self.minimums[idx]:    #in case of single-valued attribute
                            scaledEx[idx]=(self.nClassMax + self.nClassMin)/2
                        else:
                            scaledEx[idx]=self.nClassMin+(self.nClassMax-self.nClassMin)*\
                                    (value_idx-self.minimums[idx])/\
                                    (self.maximums[idx]-self.minimums[idx])
                    else:  #Disc. but don't want to scale
                        scaledEx[idx] = int(ex[classVar])

       # print ex.getclass()," -> ",scaledEx.getclass()
        return scaledEx
      
    def fastscaleEx(self,ex):
        """
        This method skips the check of self.nClassMin, self.nClassMax, nMin and nMax type and ranges.
        It assumes data was already analysed, and the example is in the correct domain format and order, for ex, using fixEx
        
        """
        # create empty example [?,?,?,...,?] of the correct domain
        scaledEx = orange.Example(self.domain)

        #Scale Attributes  -  Scale function is according to libSVM code
        for idx,attr in enumerate(ex.domain.attributes):
            if attr.varType == orange.VarTypes.Continuous:
                if self.maximums[idx] == self.minimums[idx]:   #in case of single-valued attribute
                    scaledEx[idx]=(self.nMax + self.nMin)/2
                else:
                    scaledEx[idx]=self.nMin+(self.nMax-self.nMin)*\
                            (float(ex[attr].value)-self.minimums[idx])/\
                            (self.maximums[idx]-self.minimums[idx])
            else:  # is Discrete
                value_idx = int(ex[attr])
                if self.maximums[idx] == self.minimums[idx]:    #in case of single-valued attribute
                    scaledEx[idx]=(self.nMax + self.nMin)/2
                else:
                    scaledEx[idx]=self.nMin+(self.nMax-self.nMin)*\
                            (value_idx-self.minimums[idx])/\
                            (self.maximums[idx]-self.minimums[idx])

        #Scale class
        classVar = ex.domain.classVar
        if classVar!=None:
            if ex[classVar].value=="?":
                scaledEx[classVar.name]="?"
            else:
                idx = self.varNameIdx[classVar.name]
                if classVar.varType == orange.VarTypes.Continuous:
                    if self.scaleClass:
                        if self.maximums[idx] == self.minimums[idx]:    #in case of single-valued attribute
                            scaledEx[idx]=(self.nClassMax + self.nClassMin)/2
                        else:
                            scaledEx[idx]=self.nClassMin+(self.nClassMax-self.nClassMin)*\
                                (float(ex[classVar].value)-self.minimums[idx])/\
                                (self.maximums[idx]-self.minimums[idx])
                    else: #Cont. but we dont want to scale it
                        scaledEx[idx] = ex[classVar].value
                else:  # is Discrete
                    if self.scaleClass:
                        value_idx = int(ex[classVar])
                        if self.maximums[idx] == self.minimums[idx]:    #in case of single-valued attribute
                            scaledEx[idx]=(self.nClassMax + self.nClassMin)/2
                        else:
                            scaledEx[idx]=self.nClassMin+(self.nClassMax-self.nClassMin)*\
                                    (value_idx-self.minimums[idx])/\
                                    (self.maximums[idx]-self.minimums[idx])
                    else:  #Disc. but don't want to scale
                        scaledEx[idx] = int(ex[classVar])

       # print ex.getclass()," -> ",scaledEx.getclass()
        return scaledEx
 
    def __setattr__(self,name,value):
        try:
            if name in ("nClassMin","nClassMax"):
                if value == None:
                    self.__dict__[name] = None
                else:
                    self.__dict__[name] = float(value)
            elif name in ("nMin","nMax"):                     #Float
                self.__dict__[name] = float(value)
            elif name in ("scaleClass"):                    #Bool
                self.__dict__[name] = bool(value)
            else:
                self.__dict__[name] = value                 #Just assign the value
        except: 
            if self.verbose >0: print "ERROR: Invalid value (",value,") for the  attribute ", name
            return False
        return True

    def scaleAndContinuizeData(self,data):
        """returns a dataset with new domain if success, or None if errors found"""

        if self.nMax<=self.nMin:
            if self.verbose >0: print "Attributes scaling limits are not correct!"
            return None

        if self.scaleClass:
            if self.nClassMin == None or self.nClassMax == None:
                if self.verbose >0: print "Warning: Class scaling limits were not defined. The same attributes limits will be used"
                self.nClassMax = self.nMax
                self.nClassMin = self.nMin
            elif self.nClassMax <= self.nClassMin:
                if self.verbose >0: print "Class scaling limits are not correct!"
                return None

        if not data:
            if self.verbose >0: print "No data to scale"
            return None

        if not self.dataAnalyzed:
            self.analyzeData(data)

        if not self.dataAnalyzed: 
            if self.verbose >0: print  "ERROR: Data could not be analyzed"
            return None

        # Reset the decimals and Clear any old data
        #self.dataDomain = data.domain # DEBUG ??  Remove this??  already done when calling self.analyzeData
        for attr in self.domain:
            attr.numberOfDecimals = self.numberOfDecimals
        scaledData = DataTable(self.domain)


        if self.verbose>1: print "   scaling Examples ..."
        for ex in data:
            scaledEx = scalizer.scaleEx(self,ex)
            if scaledEx:
                scaledData.append(scaledEx)  # The scaleEx will convert the example to the scalizer domain and also make the scale of the respective vars.
            else:       
                if self.verbose >0: print "ERROR: Example ",ex," could not be scaled! It was ignored and will not be present in the scaled data"
            if self.callback!=None:
                self.callback()
        return scaledData        






class scalizerX(scalizer): 
    """ This is a version of the scalizer which returns the class in an orange format
        This is the version correponding to the first(old) scalizer. Avoid using this, use the scalizer class instead"""
    def __new__(cls,**kwds):
        self = scalizer.__new__(cls, **kwds)
        return self

    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        scalizer.__init__(self,**kwds) 
        #self.Xdomain = None
        if self.file != None:
            XdomainFile = os.path.join(self.file,"Xdomain.tab")
            if not os.path.isfile(XdomainFile):
                if self.verbose >0: print "ERROR: Could not open XdominFile.tab file for reading"
                self.Xdomain = None
            else:       
                self.Xdomain = DataTable(XdomainFile).domain

    def saveScalingValues(self,file):
        scalizer.saveScalingValues(self,file)
        DataTable(self.Xdomain).save(os.path.join(file,"Xdomain.tab"))


    def scaleAndContinuizeData(self,data):
        """returns a dataset with new domain if success, or None if errors found"""
        baseScaledData = scalizer.scaleAndContinuizeData(self,data)
        if (baseScaledData == None):
            print "ERROR: data was not scaled by the scalizer class!"
            return None
        # create empty dataset with all attributes being of type float
        if data.domain.classVar==None:
            scaledData=DataTable(orange.Domain([orange.FloatVariable('%s'%x.name) for x in data.domain.attributes],False))
        else:
            if data.domain.classVar.varType==orange.VarTypes.Continuous:
                scaledData=DataTable(orange.Domain([orange.FloatVariable('%s'%x.name) for x in data.domain.attributes]+\
                                    [orange.FloatVariable(data.domain.classVar.name)]))
                scaledData.domain.classVar.numberOfDecimals = data.domain.classVar.numberOfDecimals
            elif hasattr(data.domain.classVar,"values"):
                scaledData=DataTable(orange.Domain([orange.FloatVariable('%s'%x.name) for x in data.domain.attributes]+\
                                    [orange.EnumVariable(data.domain.classVar.name)]))
                scaledData.domain.classVar.values = [str(x) for x in range(len(self.values[-1]))]
            else:
                if self.verbose >0: print "ERROR: Class var is of invalid type"
                return none
        self.Xdomain = scaledData.domain
        for attr in scaledData.domain.attributes:
            if data.domain[attr.name].varType==orange.VarTypes.Continuous:
                attr.numberOfDecimals = data.domain[attr.name].numberOfDecimals
            else:
                attr.numberOfDecimals = self.numberOfDecimals

        classVar = data.domain.classVar
        #The conversion neede is only for the class attribute which will have to be converted to cont. or disc.
        #We know that scaledData and baseScaledData may only differ in the type of class variable
        for ex in baseScaledData:
            if scaledData.domain.classVar.varType == orange.VarTypes.Discrete:
                scaledData.append(orange.Example(scaledData.domain,[attr.value for attr in ex.native()[:-1]]+["?"]))
                #convert the class back to the original value since in this version the discrete classes are not to be scaled.
                value = str(int(scalizer.convertClass(self,ex.getclass().value)))
                if value in scaledData.domain.classVar.values:
                    scaledData[-1].setclass(value)
            else:
                scaledData.append(orange.Example(scaledData.domain,[attr.value for attr in ex.native()]))
            if self.callback!=None:
                self.callback()
        return scaledData        


    def convertClass(self,classDomain, value):
        """ Converts the variable 'value' to the variable of type 'classDomain' after converting back any scale eventually made
            The 'value' parameter must be an orange.Value object
            The return value is also an orange.Value object of type 'classDomain'"""
        if not self.dataAnalyzed:
            if self.verbose >0: print  "ERROR: Data was not analyzed yet"
            return False
        if not hasattr(value,"value") or type(value)!=orange.Value:
            if self.verbose >0: print "ERROR: ",value," is not in correct format!"
            return False
        
        res = orange.Value(classDomain)

        if res.variable.name not in self.varNames:
            if self.verbose >0: print "ERROR: Variable ",res.variable.name," unknown by the scalizer"
            return False
        else:
            idxVar=self.varNameIdx[res.variable.name]

        if classDomain.varType == orange.VarTypes.Continuous:
            if self.scaleClass:
                res.value = scalizer.convertClass(self,value.value)
            else:
                res.value = value.value
        elif hasattr(classDomain,"values"):
            for val in res.variable.values:
                if val not in self.values[idxVar]:
                    if self.verbose >0: print "ERROR: Values ",res.variable.values," are not present in scalizer"
                    return False
            #class was not scaled since it is discrete
            rValue = int(round(float(value.value)))

            if rValue >= (len(self.values)-1):
                res.value = self.values[idxVar][len(self.values[idxVar])-1]
            elif rValue <=0:
                res.value = self.values[idxVar][0]
            else:
                res.value = self.values[idxVar][rValue]
        else:
            if self.verbose >0: print "ERROR: class type is not Continuous nor Discrete"
            return False

        return res
        
    def scaleEx(self,ex):
        """Scales the example 'ex' to an example will all variables Cont. and the class of same type of 
           original class used when analyzing"""
        scaledEx = scalizer.scaleEx(self,ex) # This is always  convertion to a continuous variable.
        #The conversion needed is only for the class attribute which will have to be converted to cont. or disc.
        #We know that scaledEx and self.Xdomain may only differ in the type of class variable. And scaledEx class is always Cont.
        if self.Xdomain.classVar.varType == orange.VarTypes.Discrete:
            fixedEx = orange.Example(self.Xdomain,[attr.value for attr in scaledEx.native()[:-1]]+["?"])
            value = str(ex.getclass().value)
            if value in self.values[-1]:
                fixedEx.setclass(self.values[-1].index(value))
        else:
            fixedEx = orange.Example(self.Xdomain,[attr.value for attr in scaledEx.native()])
        return fixedEx


##ecPA          

def getResponseType(dataFile):
    """
    Return the string 'Regression' or 'Classification' depending on the type of response in dataFile. 
    """
    data = DataTable(dataFile)
    if data.domain.classVar.varType == orange.VarTypes.Discrete:
        return "Classification"
    else:
        return "Regression"

      
def getApproxMemReq(filePath):
    """
    Estimate the required memory (in MB) using the number of elements in the data matrix, N.
    mem = N*8(float representation)*2(Orange data object)*5(Copies of the data matrix during param opt. 
    This number needs to be refined.)
    filePath is the location of an Orange data set on disk
    """
    data = DataTable(filePath)

    nEx = len(data)
    nAttr = len(data.domain.attributes)

    memReq = int(nEx*nAttr*8*2*5*0.000001)
    if memReq < 150:
        memReq = 150

    return memReq     

def rmClassVar(data):
    """
    Make the class varible an attribute and the domain classless
    It also removes any meta-attributes
    """
    newDomain = orange.Domain(data.domain, 0)    
    newData = DataTable(newDomain, data)

    return newData

