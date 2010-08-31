import unittest
import os
import time

import orange
from trainingMethods import AZorngPLS
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC
import AZorngTestUtil
import orngImpute


class dataUtilitiesTest(unittest.TestCase):

    def setUp(self):
        """Sets up the test """
        unusedValuesDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_UnusedValues_Train.tab")
        multiClassDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/iris.tab")
        missValsDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/testMissVals.tab")
        testDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_SmallTest.tab")
        badVarTypePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarType.tab")#The var is discrete but should be continuous
        badVarNamePath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarName.tab")
        badVarOrderPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarOrder.tab")
        self.badVarOrderValuesPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarOrderValues.tab")
        badVarCountPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarCount.tab")
        wMetaPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_W_metas_SmallTest.tab")
        badVarType2Path = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarType2.tab")#there are 2 vars discrete  but should be continuous
        badVarTypeSOMEPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_BadVarTypeSOME.tab")#there are 2 vars discrete  but should be continuous and the order is not correct

        # Read in the data
        self.unusedValuesData = dataUtilities.DataTable(unusedValuesDataPath,createNewOn=orange.Variable.MakeStatus.OK)
        self.multiClassData = dataUtilities.DataTable(multiClassDataPath,createNewOn=orange.Variable.MakeStatus.OK)
        self.missValsData = dataUtilities.DataTable(missValsDataPath,createNewOn=orange.Variable.MakeStatus.OK)
        self.testData = dataUtilities.DataTable(testDataPath,createNewOn=orange.Variable.MakeStatus.OK)
        self.badVarTypeData = dataUtilities.DataTable(badVarTypePath,createNewOn=orange.Variable.MakeStatus.OK)
        self.badVarType2Data = dataUtilities.DataTable(badVarType2Path,createNewOn=orange.Variable.MakeStatus.OK)
        self.badVarTypeSOMEData = dataUtilities.DataTable(badVarTypeSOMEPath,createNewOn=orange.Variable.MakeStatus.OK)
        self.badVarNameData = dataUtilities.DataTable(badVarNamePath,createNewOn=orange.Variable.MakeStatus.OK)
        self.badVarOrderData = dataUtilities.DataTable(badVarOrderPath,createNewOn=orange.Variable.MakeStatus.OK)
        # DO NOT LOAD badVarOrderValuesPath at this time!!!!
        self.badVarCountData = dataUtilities.DataTable(badVarCountPath,createNewOn=orange.Variable.MakeStatus.OK)   #One less example 
        self.wMetaData = dataUtilities.DataTable(wMetaPath,createNewOn=orange.Variable.MakeStatus.OK)

    def test_changeClassName(self):
        """Test the changeClassName method"""
        self.assertEqual(self.testData.domain.classVar.name,"Activity")
        newData = dataUtilities.changeClassName(self.testData,"NewClassName")

        #Test if the class values remain the same
        for idx,ex in enumerate(self.testData):
            self.assertEqual(str(self.testData[idx]["Activity"].value) , str(newData[idx]["NewClassName"].value))

        #Assure that the original data was not changed
        self.assertEqual(self.testData.domain.classVar.name,"Activity")
        #Check if the name is correct
        self.assertEqual(newData.domain.classVar.name,"NewClassName")

        #Change the class name to an existing name in attributes
        newData = dataUtilities.changeClassName(self.testData,"Level")
        self.assertEqual(newData.domain.classVar.name,"Level_1")

        # the class name will not be tested for name duplication
        newData = dataUtilities.changeClassName(newData,"Level")
        self.assertEqual(newData.domain.classVar.name,"Level_1")
        newData = dataUtilities.changeClassName(self.testData,"Activity")
        self.assertEqual(newData.domain.classVar.name,"Activity")

        # If something is missing, None will be returned
        newData = dataUtilities.changeClassName(None,"NAME")
        self.assertEqual(newData,None)
        newData = dataUtilities.changeClassName(self.testData,"")
        self.assertEqual(newData,None)
        
        # Double check if the original data was not changed
        self.assertEqual(self.testData.domain.classVar.name,"Activity")



    def test_getDataWithoutUnusedValues(self):
        """Test the getDataWithoutUnusedValues method"""
        # on Data self.unusedValuesData the attrs RING_HERG and Activity have unused Values!
        learner = AZorngPLS.PLSLearner()

        valuesSel = str(self.unusedValuesData.domain["Attr3"].values)
        valuesActivity = str(self.unusedValuesData.domain["Activity"].values)
        
        #Check that there are in fact unused values: "OTHER" and <7, 8>
        self.assert_(valuesActivity == "<POS, NEG, OTHER>")
        self.assert_(valuesSel == "<1, 2, 3, 4, 5, 6, 7, 8, 9>")
        newData = dataUtilities.getDataWithoutUnusedValues(self.unusedValuesData,True)
        #Check on new dataset that the unused values were removed from the domain
        self.assert_(str(newData.domain["Attr3"].values) == "<1, 2, 3, 4, 5>")
        self.assert_(str(newData.domain["Activity"].values) == "<POS, NEG>")

        #Check that the original domain still exists and has not been modified
        self.assert_(str(self.unusedValuesData.domain["Attr3"].values) == "<1, 2, 3, 4, 5, 6, 7, 8, 9>")
        self.assert_(str(self.unusedValuesData.domain["Activity"].values) == "<POS, NEG, OTHER>")
 



    def test_AddDiscreteClass(self):
        """Tes the add discrete class"""
        newData = dataUtilities.addDiscreteClass(self.testData,"Class_Test","Val_1",["Val_2","Val_3"])
        self.assertEqual(len(self.testData.domain)+1,len(newData.domain),"Domain len not correct")
        self.assertEqual(newData.domain.classVar.name,"Class_Test","Wrong Class var")
        self.assertEqual(str(newData.domain.classVar.values),"<Val_2, Val_3, Val_1>","Wrong values for the class var")
        self.assertEqual(len(self.testData),len(newData),"Missing Examples")

    def test_miss_values_keywords(self):
        """Test the C layer handle of miss values keywords"""
        #dataUtilities.DClist = ["","?"," ","NA","ERROR", "."]
        #dataUtilities.DC = "?"
        #dataUtilities.DKlist = ["~","*","MIS_VAL"]
        #dataUtilities.DK = "~"

        self.assertEqual(len(self.missValsData),150,"Missing Examples")
        self.assertEqual(str(self.missValsData.domain),"[sepal length, sepal width, petal length, petal width, iris]")
        self.assertEqual(str(self.missValsData[0]),"[?, 3.5, 1.4, '?', 'Iris-setosa']")
        self.assertEqual(len(self.missValsData),150,"Missing Examples")
        for i in range(0,7):
            self.assertEqual(self.missValsData[i]["sepal length"],'?',\
            "Attr of example "+str(i)+" Should be a missing value!")
            self.assertEqual(self.missValsData[i]["petal width"],'?',\
            "Attr of example "+str(i)+" Should be a missing value!")
        for i in range(7,11):
            self.assertEqual(self.missValsData[i]["sepal length"],'~',\
            "Attr of example "+str(i)+" Should be a missing value!")
            self.assertEqual(self.missValsData[i]["petal width"],'~',\
            "Attr of example "+str(i)+" Should be a missing value!")


    def test_merge(self):
        """ Test the mergeData functionality"""
        data1 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D1.tab"))
        data2 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D2.tab"))

        mergeAttrName1 = "nameA"
        mergeAttrName2 = "nameB"


        res = dataUtilities.horizontalMerge(data1, data2, mergeAttrName1, mergeAttrName2)
        testD = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D1p2.tab"))
        for idx,ex in enumerate(testD):
            for attr in res[idx].domain:
                self.assertEqual(res[idx][attr.name],ex[attr.name],"D1p2: Example " + str(idx) + " Attr "+attr.name+" value do not match:"+str(res[idx][attr.name])+" != "+str(ex[attr.name]))

        res = dataUtilities.horizontalMerge(data2, data1, mergeAttrName2, mergeAttrName1)
        testD = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D2p1.tab"))
        for idx,ex in enumerate(testD):
            for attr in res[idx].domain:
                self.assertEqual(res[idx][attr.name],ex[attr.name],"D2p1: Examples values do not match")

    def test_merge_classHandle(self):
        """ Test the handle of class var in merge"""
        data1 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D1.tab"))
        data2 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D2.tab"))

        mergeAttrName1 = "nameA"
        mergeAttrName2 = "nameB"
        
        self.assert_(data1.domain.classVar.name != data2.domain.classVar.name,"D1 and D2 should have different classes in order to perform the test")         
        #Test Both with Class -> Class should be the same as data1
        res = dataUtilities.horizontalMerge(data1, data2, mergeAttrName1, mergeAttrName2)
        self.assertEqual(res.domain.classVar.name,data1.domain.classVar.name,"Wrong Class var. Got: "+str(res.domain.classVar.name)+" Expected: "+str(data1.domain.classVar.name))

        #Test data2 classless -> Class should be the same as data1
        data2NoClass = orange.ExampleTable(orange.Domain(data2.domain,None),data2)
        self.assertEqual(data2NoClass.domain.classVar,None,"Data should not have class")
        res = dataUtilities.horizontalMerge(data1, data2NoClass, mergeAttrName1, mergeAttrName2)
        self.assertEqual(res.domain.classVar.name,data1.domain.classVar.name,"Wrong Class var")
        #Test now with different order od datasets
        res = dataUtilities.horizontalMerge(data2NoClass, data1, mergeAttrName2, mergeAttrName1)
        self.assertEqual(res.domain.classVar.name,data1.domain.classVar.name,"Wrong Class var")

        # Test Both classless -> res should also be classless
        data1NoClass = orange.ExampleTable(orange.Domain(data1.domain,None),data1)
        self.assertEqual(data1NoClass.domain.classVar,None,"Data should not have class")
        res = dataUtilities.horizontalMerge(data1NoClass, data2NoClass, mergeAttrName1, mergeAttrName2)
        self.assertEqual(res.domain.classVar, None,"Wrong Class var")

        # Test when both have same class name if they don't get exchanged
        data1 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D1.tab"))
        data2 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D11.tab"))
        self.assertEqual(data1.domain.classVar.name , data2.domain.classVar.name,"D1 and D2 should have same classes in order to perform the test")
        res = dataUtilities.horizontalMerge(data1, data2, "nameA", "nameA")
        self.assertEqual(res.domain.classVar.name,"Activity")
        self.assertEqual(res[3]["nameA"],"N4")
        self.assertEqual(res[3][res.domain.classVar],"POS")
        self.assertEqual(res[3]["Activity_1"],"NEG")
        #Test now with different order od datasets
        res = dataUtilities.horizontalMerge(data2, data1, "nameA", "nameA")
        self.assertEqual(res.domain.classVar.name,"Activity")
        self.assertEqual(res[3]["nameA"],"N4")
        self.assertEqual(res[3][res.domain.classVar],"NEG")
        self.assertEqual(res[3]["Activity_1"],"POS")

        

    def test_merge_stringVar(self):
        """ Test the mergeData functionality based on string var"""
        data1 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/SEIMLS10.tab"))
        data2 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/SEIMLS10_2.tab"))

        mergeAttrName1 = "SEIMLS"
        mergeAttrName2 = "SEIMLS"


        res = dataUtilities.horizontalMerge(data1, data2, mergeAttrName1, mergeAttrName2)
        testD = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/MergedString.tab"))
        for idx,ex in enumerate(testD):
            for attr in res[idx].domain:
                if res[idx][attr.name].isSpecial():
                    self.assertEqual(res[idx][attr.name],ex[attr.name],"MergedStrings: Example " + str(idx) + " Attr "+attr.name+" value do not match:"+str(res[idx][attr.name])+" != "+str(ex[attr.name]))
                elif res[idx][attr.name].varType == orange.VarTypes.String or ex[attr.name].varType == orange.VarTypes.String:
                    self.assertEqual(res[idx][attr.name],ex[attr.name],"MergedString: Examples values do not match")
                else:
                    self.assertEqual(float(res[idx][attr.name]),float(ex[attr.name]),"MergedString: Examples values do not match")
                    self.assertEqual(int(res[idx][attr.name]),int(ex[attr.name]),"MergedString: Examples values do not match")

    def test_merge_With_Metas(self):
        """Test the merge with data having meta Attributes"""        
        data1 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_W_metas_SmallTest.tab"))
        data2 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_No_metas_SmallTest.tab"))
        
        mergeAttrName1 = "Attr3"
        mergeAttrName2 = "Attr3"
        res = dataUtilities.horizontalMerge(data1, data2, mergeAttrName1, mergeAttrName2)
        self.assertEqual(res[0]["Comments"],"ok","Got: "+str(res[0]["Comments"]))
        res = dataUtilities.horizontalMerge(data2, data1, mergeAttrName2, mergeAttrName1)
        self.assertEqual(res[5]["Comments"],"notok","Got: "+str(res[5]["Comments"]))

        #Teste merge based on Meta Attrs
        data2 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/BinClass_W_metas_SmallTest.tab"))        
        mergeAttrName1 = "Comments"
        mergeAttrName2 = "Comments"
        res = dataUtilities.horizontalMerge(data2, data1, mergeAttrName2, mergeAttrName1)
        self.assertEqual(res[0]["Comments"],"ok","Got: "+str(res[0]["Comments"]))
        self.assertEqual(res[5]["Comments"],"notok","Got: "+str(res[5]["Comments"]))




    def test_merge_diff_Attr(self):
        """ Test the mergeData functionality"""
        data1 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D11.tab"))
        data2 = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D22.tab"))

        mergeAttrName1 = "Attr99D"
        mergeAttrName2 = "Attr9C"


        res = dataUtilities.horizontalMerge(data1, data2, mergeAttrName1, mergeAttrName2)
        testD = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mergeTestsData/D1p22.tab"))
        for idx,ex in enumerate(testD):
            for attr in res[idx].domain:
                if res[idx][attr.name].isSpecial():
                    self.assertEqual(res[idx][attr.name],ex[attr.name],"D1p22:Example " + str(idx) + " Attr "+attr.name+" value do not match:"+str(res[idx][attr.name])+" != "+str(ex[attr.name]))
                else:
                    self.assertEqual(float(res[idx][attr.name]),float(ex[attr.name]),"D1p22: Examples values do not match")
                    self.assertEqual(int(res[idx][attr.name]),int(ex[attr.name]),"D1p22: Examples values do not match")


    def test_concatenate(self):
        """ Test the concatenate functionality"""
        #Avtivity = NEG; Attrs: 1,2,3,6,4; attr3=Discrete
        file1=os.path.join(AZOC.AZORANGEHOME,"tests/source/data/concatenateTestsData/D1.tab")
        #Activity = POS and NEG; Attrs: 1,7,3,4,5; attr3=Continuous should be discrete; 
        #Attr1=Discrete should be continuous; Attr4=Discrete should be continuous
        file2=os.path.join(AZOC.AZORANGEHOME,"tests/source/data/concatenateTestsData/D2.tab")
        #Activity = 1 and 0 Coninuous should be discrete; Attrs: 8,2,3,6,4; attr3=Discrete; 
        file3=os.path.join(AZOC.AZORANGEHOME,"tests/source/data/concatenateTestsData/D3.tab")

        D1=dataUtilities.DataTable(file1)
        D2=dataUtilities.DataTable(file2)
        D3=dataUtilities.DataTable(file3)

        res=dataUtilities.concatenate([D1,D2,D3],mergeDomains=True,useFirstAsLeader=False)
        self.assertEqual(len(res[0].domain),9,"C1: Wrong number of attributes")
        self.assertEqual(len(res[0]),30,"C1: Wrong number of examples")
        testD=dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/concatenateTestsData/concatenated1.tab"))
        for idx,ex in enumerate(testD):
            for attr in res[0][idx].domain:
                self.assertEqual(res[0][idx][attr.name],ex[attr.name],"C1: Examples values do not match")

        res=dataUtilities.concatenate([D2,D1,D3],mergeDomains=False,useFirstAsLeader=False)
        self.assertEqual(len(res[0].domain),3,"C2: Wrong number of attributes")
        self.assertEqual(len(res[0]),30,"C2: Wrong number of examples")
        testD=dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/concatenateTestsData/concatenated2.tab"))
        for idx,ex in enumerate(testD):
            for attr in res[0][idx].domain:
                self.assertEqual(res[0][idx][attr.name],ex[attr.name],"C2: Examples values do not match")

        res=dataUtilities.concatenate([D1,D2,D3],mergeDomains=True,useFirstAsLeader=True)
        self.assertEqual(len(res[0].domain),6,"C3: Wrong number of attributes")
        self.assertEqual(len(res[0]),30,"C3: Wrong number of examples")
        testD=dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/concatenateTestsData/concatenated3.tab"))
        for idx,ex in enumerate(testD):
            for attr in res[0][idx].domain:
                self.assertEqual(res[0][idx][attr.name],ex[attr.name],"C3: Examples values do not match")

    def test_Duplicates_and_Same_Attr(self):
        """Test Duplicates and Same attributes fix"""
        data = orange.ExampleTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/DupAndSameVars.tab"))
        origDataDLen = len(data.domain)
        self.assertEqual(origDataDLen,16,"Wrong original domain len: "+str(len(data.domain)))
        # Actual Duplicates
        self.assertEqual("{'MNOLAME': [1, 4, 5, 8, 14], 'SEIMLS_1': [0, 7, 13]}",\
            str(dataUtilities.findDuplicatedNames(data.domain)),"Wrong Duplicates on original Dataset:"+str(dataUtilities.findDuplicatedNames(data.domain)))
        # 1. Remove SAMES
        sames=dataUtilities.getDomainWithoutDupSameAttr(data.domain) 
        data=orange.ExampleTable(sames[1],data)
        self.assertEqual(eval("{'MNOLAME': 4, 'SEIMLS_1': 1}"),sames[0],"Error removing: "+str(sames))
        # 2. Fix Dup Names
        dups = dataUtilities.fixDuplicatedNames(data.domain)
        self.assertEqual(eval("{'SEIMLS_1': [0, 9]}"),dups,"Error fixing Duplicate names: "+str(dups))
        self.assertEqual(eval("{}"),dataUtilities.findDuplicatedNames(data.domain),"Duplicate names not fixed")
        self.assertEqual(len(data.domain),11,"Wrong fixed domain len: "+str(len(data.domain)))
        self.assertEqual(str(data.domain),"[SEIMLS_1_1, MNOLAME, number, ACDC, window, nonSel, winter, bio, specif, SEIMLS_1, Activity]","Wrong attribute names: "+str(data.domain))


    def test_SMI(self):
        """Test loader os SMI files"""
        SMIdata = dataUtilities.loadSMI(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/sample.smi"))
        self.assertEqual(len(SMIdata),2500)
        self.assertEqual(str(SMIdata.domain),"[SMILES, MOLNAME]")

        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirSMItest"+str(time.time()))
        os.mkdir(scratchdir)
        savedSMI=os.path.join(scratchdir,"savedSMI.tab")
        SMIdata.save(savedSMI)
        testD=dataUtilities.DataTable(savedSMI)
        for idx,ex in enumerate(testD):
            for attr in testD.domain:
                self.assertEqual(SMIdata[idx][attr.name],ex[attr],"Examples values do not match")

        if os.path.isfile(savedSMI):
            os.system("rm -rf "+scratchdir)


    def test_scalizer(self):
        """Test the scalizer """
        scaler = dataUtilities.scalizer(data=self.testData)     # The new scalizer 
        scalerX = dataUtilities.scalizerX(data=self.testData)   # derived from scalizer which mimics the old scalizer version for old models
        scaledData = scaler.scaleAndContinuizeData(self.testData)
        scaledDataX = scalerX.scaleAndContinuizeData(self.testData)
        self.assertEqual("[-0.2727273, 0.0000000, 0.0000000, -0.2000000, -1.0000000, -0.3183614, -1.0000000, 1.0000000, -1.0000000, -0.3333333, 1.0000000]",str(scaledData[0]))
        self.assertEqual("[-0.272727, 0.000000, 0.000000, -0.200000, -1.000000, -0.318361, -1.0000000, 1.0000000, -1.0000000, -0.3333333, '1']",str(scaledDataX[0]))
        self.assertEqual("[-1.0000000, 0.0000000, 0.0000000, -0.2000000, -1.0000000, 1.0000000, -1.0000000, 1.0000000, -0.5000000, 0.3333333, 1.0000000]",str(scaledData[6]))
        self.assertEqual("[-1.000000, 0.000000, 0.000000, -0.200000, -1.000000, 1.000000, -1.0000000, 1.0000000, -0.5000000, 0.3333333, '1']",str(scaledDataX[6]))

        self.assertEqual(scalerX.convertClass(self.testData.domain.classVar,scaledData[0].getclass()),self.testData[0].getclass())
        self.assertEqual(scalerX.convertClass(self.testData.domain.classVar,scaledData[5].getclass()),self.testData[5].getclass())
        XscaledEx = scalerX.scaleEx(self.testData[0]) 
        self.assertEqual(str(XscaledEx),"[-0.272727, 0.000000, 0.000000, -0.200000, -1.000000, -0.318361, -1.0000000, 1.0000000, -1.0000000, -0.3333333, '1']")
        self.assertEqual(scalerX.convertClass(self.testData.domain[self.testData.domain.classVar],XscaledEx[XscaledEx.domain.classVar]).value,self.testData[0][self.testData.domain.classVar])

        XscaledEx = scalerX.scaleEx(self.testData[1])
        self.assertEqual(str(XscaledEx),"[-1.000000, 0.000000, 0.000000, 0.200000, -1.000000, -0.676275, 1.0000000, -1.0000000, -0.5000000, -1.0000000, '1']")
        #self.assertEqual(scalerX.convertClass(self.testData.domain[self.testData.domain.classVar],XscaledEx[XscaledEx.domain.classVar]).value,self.testData[1][self.testData.domain.classVar])

        scaler2=dataUtilities.scalizer()
        scaler2.analyzeData(self.testData)
        new=scaler2.scaleAndContinuizeData(self.testData)

        scaler2.scaleClass=True
        new2=scaler2.scaleAndContinuizeData(self.testData)
        self.assertEqual("[-0.2727273, 0.0000000, 0.0000000, -0.2000000, -1.0000000, -0.3183614, -1.0000000, 1.0000000, -1.0000000, -0.3333333, 1.0000000]",str(new2[0]))
        self.assertEqual("[-1.0000000, 0.0000000, 0.0000000, -0.6000000, -1.0000000, -0.0414842, 1.0000000, -1.0000000, -1.0000000, -1.0000000, 0.0000000]",str(new[5]))

        scaler2.nMax=23
        scaler2.nMin=-3
        scaler2.nClassMax=2
        scaler2.nClassMin=-5
        new3=scaler2.scaleAndContinuizeData(self.testData)
        self.assertEqual("[6.4545455, 10.0000000, 10.0000000, 7.4000001, -3.0000000, 5.8613019, -3.0000000, 23.0000000, -3.0000000, 5.6666665, 2.0000000]",str(new3[0]))
        self.assertEqual("[-3.0000000, 10.0000000, 10.0000000, 2.2000000, -3.0000000, 9.4607048, 23.0000000, -3.0000000, -3.0000000, -3.0000000, -5.0000000]",str(new3[5]))

        self.assertEqual(scaler2.convertClass(new3[0].getclass()),int(self.testData[0].getclass()))     
        self.assertEqual(scaler2.convertClass(new3[5].getclass()),int(self.testData[5].getclass()))   

        
        self.assertEqual("['Measure', '[Br]([C])', '[N]([N])', '[O]([C])', '[C]([C][F])', 'Level', 'DiscAttr1', 'DiscAttr2', 'Attr3', 'YetOther', 'Activity']",str(scaler2.varNames))
        self.assertEqual("[[], [], [], [], [], [], ['Red', 'Green', 'Blue'], ['YES', 'NO'], ['1', '2', '3', '4', '5'], ['B', 'A', 'C', '1'], ['POS', 'NEG']]",str(scaler2.values))
        self.assertEqual("[11.0, 0.0, 0.0, 5.0, 1.0, 5.6697301864624023, 2.0, 1.0, 4.0, 3.0, 1.0]",str(scaler2.maximums))
        self.assertEqual("[0.0, 0.0, 0.0, 0.0, 0.0, 0.035923998802900314, 0.0, 0.0, 0.0, 0.0, 0.0]",str(scaler2.minimums))

        ex = scaler2.scaleEx(self.testData[3])
        self.assertEqual(str(new3[3]),str(ex))


        #Test save/load scaling values file
        scratchdir = os.path.join(AZOC.SCRATCHDIR, "scratchdirScalizerTest"+str(time.time()))
        os.mkdir(scratchdir)
        scalingFile =os.path.join(scratchdir,"scaling.att")

        scaler2.saveScalingValues(scalingFile)
        scaledData2 = scaler2.scaleAndContinuizeData(self.testData)
        scaler3=dataUtilities.scalizer(file=scalingFile)
        scaledData3 = scaler3.scaleAndContinuizeData(self.testData)


        self.assertEqual(str(scaledData2[0]),str(scaledData3[0]))
        self.assertEqual(str(scaledData2[2]),str(scaledData3[2]))
        self.assertEqual(str(scaledData2[4]),str(scaledData3[4]))

        scaler4=dataUtilities.scalizer(file=scalingFile)
        self.assertEqual(str(scaledData2[2]),str(scaler4.scaleEx(self.testData[2])))

        #Test with continuous class
        data = dataUtilities.DataTable(os.path.join(AZOC.AZORANGEHOME,"tests/source/data/Reg_No_metas_FullNumeric_Train.tab"))
        scaler4 = dataUtilities.scalizer(data=data)
        scalerX = dataUtilities.scalizer(data=data, scaleClass = True)
        scaledDataX = scalerX.scaleAndContinuizeData(data)
        scaledData4 = scaler4.scaleAndContinuizeData(data)
        self.assertEqual(str(scaledDataX[0]),"[-1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -0.7500000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, 1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -0.5000000, -1.0000000, -0.7500000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, 0.0000000, -0.7500000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, 1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, 1.0000000, -1.0000000, -0.7647059, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -0.3333333, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, 0.0000000, -1.0000000, -1.0000000, -0.6000000, 0.0000000, -1.0000000, -1.0000000, -0.3333333, -1.0000000, -1.0000000, -1.0000000, -0.3333333, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -1.0000000, -0.4086332]")
        self.assertEqual(round(scalerX.convertClass(scaledDataX[0].getclass().value),5),round(data[0].getclass().value,5))
        self.assertEqual(scaler4.convertClass(scaledData4[0].getclass().value),data[0].getclass())
        self.assertEqual(scaler4.convertClass(scaledData4[6].getclass().value),data[6].getclass())

        scaler4.scaleClass = True
        newData = scaler4.scaleAndContinuizeData(data)
        self.assertEqual(scaler4.convertClass(newData[0].getclass()),data[0].getclass())
        self.assertEqual(round(scaler4.convertClass(newData[6].getclass()),6),round(data[6].getclass(),6))  #0.470254
       
        multiclassScaler = dataUtilities.scalizer(data=self.multiClassData)
        scaledData = multiclassScaler.scaleAndContinuizeData(self.multiClassData)

        self.assert_(str(scaledData[0]) == str(multiclassScaler.scaleEx(self.multiClassData[0])) == "[-0.5555557, 0.2500000, -0.8644068, -0.9166667, 0.0000000]")
        self.assert_(str(scaledData[50]) == str(multiclassScaler.scaleEx(self.multiClassData[50])) == "[0.4999999, 0.0000000, 0.2542372, 0.0833333, 1.0000000]")
        self.assert_(str(scaledData[-1]) == str(multiclassScaler.scaleEx(self.multiClassData[-1])) == "[-0.1111111, -0.1666667, 0.3898304, 0.4166666, 2.0000000]")
        #"Classes os scaledData converted back to original"
        self.assertEqual(multiclassScaler.convertClass(scaledData[0].getclass().value),0.0)
        self.assertEqual(multiclassScaler.convertClass(scaledData[50].getclass().value),1.0)
        self.assertEqual(multiclassScaler.convertClass(scaledData[-1].getclass().value),2.0)

#self.assert_(str() == str() == "")
        #"===========OLD SCALIZER============="
        #"scaledData using scaleAndContinuizeData"
        multiclassScalerX = dataUtilities.scalizerX(data=self.multiClassData)
        scaledData = multiclassScalerX.scaleAndContinuizeData(self.multiClassData)
        self.assert_(str(scaledData[0]) == str(multiclassScalerX.scaleEx(self.multiClassData[0])) == "[-0.6, 0.2, -0.9, -0.9, '0']")
        self.assert_(str(scaledData[50]) == str(multiclassScalerX.scaleEx(self.multiClassData[50])) == "[0.5, 0.0, 0.3, 0.1, '1']")
        self.assert_(str(scaledData[-1]) == str(multiclassScalerX.scaleEx(self.multiClassData[-1])) == "[-0.1, -0.2, 0.4, 0.4, '2']")
        # "Classes os scaledData converted back to original"
        self.assertEqual(multiclassScalerX.convertClass(self.multiClassData.domain.classVar,scaledData[0].getclass()).value,"Iris-setosa")
        self.assertEqual(multiclassScalerX.convertClass(self.multiClassData.domain.classVar,scaledData[50].getclass()).value,"Iris-versicolor")
        self.assertEqual(multiclassScalerX.convertClass(self.multiClassData.domain.classVar,scaledData[-1].getclass()).value,"Iris-virginica")


        # "test MultiClass with class Scale"
        data = dataUtilities.DataTable(orange.Domain(self.multiClassData.domain.variables.native()))
        for ex in self.multiClassData:
            data.append(ex)
        multiclassScaler = dataUtilities.scalizer(scaleClass = True, data=data)
        scaledData = multiclassScaler.scaleAndContinuizeData(data)
        # "===========NEW SCALIZER============="
        # "scaledData using scaleAndContinuizeData"
        self.assert_(str(scaledData[0]) == str(multiclassScaler.scaleEx(self.multiClassData[0])) == "[-0.5555557, 0.2500000, -0.8644068, -0.9166667, -1.0000000]")
        self.assert_(str(scaledData[50]) == str(multiclassScaler.scaleEx(self.multiClassData[50])) == "[0.4999999, 0.0000000, 0.2542372, 0.0833333, 0.0000000]")
        self.assert_(str(scaledData[-1]) == str(multiclassScaler.scaleEx(self.multiClassData[-1])) == "[-0.1111111, -0.1666667, 0.3898304, 0.4166666, 1.0000000]")
        # "Classes os scaledData converted back to original"
        self.assertEqual(multiclassScaler.convertClass(scaledData[0].getclass().value),0.0)
        self.assertEqual(multiclassScaler.convertClass(scaledData[50].getclass().value),1.0)
        self.assertEqual(multiclassScaler.convertClass(scaledData[-1].getclass().value),2.0)

        multiclassScalerX = dataUtilities.scalizerX(verbose = 2, scaleClass = True, data=data)
        scaledData = multiclassScalerX.scaleAndContinuizeData(data)
        #"===========OLD SCALIZER============="
        #"scaledData using scaleAndContinuizeData"
        self.assert_(str(scaledData[0]) == str(multiclassScalerX.scaleEx(self.multiClassData[0])) == "[-0.6, 0.2, -0.9, -0.9, '0']")
        self.assert_(str(scaledData[50]) == str(multiclassScalerX.scaleEx(self.multiClassData[50])) == "[0.5, 0.0, 0.3, 0.1, '1']")
        self.assert_(str(scaledData[-1]) == str(multiclassScalerX.scaleEx(self.multiClassData[-1])) == "[-0.1, -0.2, 0.4, 0.4, '2']")
        # "Classes os scaledData converted back to original"
        self.assertEqual(multiclassScalerX.convertClass(data.domain.classVar,scaledData[0].getclass()).value,"Iris-setosa")
        self.assertEqual(multiclassScalerX.convertClass(data.domain.classVar,scaledData[50].getclass()).value,"Iris-versicolor")
        self.assertEqual(multiclassScalerX.convertClass(data.domain.classVar,scaledData[-1].getclass()).value,"Iris-virginica")
        

        
        os.system("rm -rf "+ scratchdir)

 
    def test_ExFix(self):
        """Test the fix of an example according to a specific domain (with order check)"""
        fixLog={}
        ex1Str=str(self.badVarTypeData[12])
        
        ExampleFix = dataUtilities.ExFix(self.testData.domain,fixLog,True)
        exFixed1 = ExampleFix.fixExample(self.badVarTypeData[12])
        self.assert_(ex1Str.replace("'","")==str(exFixed1).replace("'",""),
                        "Example fixed changed relative to the original!\n"+ex1Str+"\n"+str(exFixed1)+"\n")
        self.assert_(exFixed1.domain==self.testData.domain,"Fixed Domain 1 does not mathch to the requested domain")

        exFixed2 = ExampleFix.fixExample(self.badVarOrderData[0])
        self.assert_(exFixed2.domain==self.testData.domain,"Fixed Domain 2 does not mathch to the requested domain")
        self.assert_(str(self.badVarOrderData[0])=="[4.000000, 0.000000, 0.000000, 2.000000, 0.000000, 1.956034, 'Red', 'NO', 'A', '1', 'NEG']",
                        "Original Example is not as expectedl!")
        self.assert_(str(exFixed2)=="[4.000000, 0.000000, 0.000000, 2.000000, 0.000000, 1.956034, 'Red', 'NO', '1', 'A', 'NEG']",
                        "Example fixed is not as expectedl!")

        fixLog.clear()
        for ex in self.badVarOrderData[3:]:
                self.assert_(ExampleFix.fixExample(ex)!=None)

        fixLog.clear()
        for ex in self.badVarTypeData[3:]:
                self.assert_(ExampleFix.fixExample(ex)!=None)
        self.assertEqual(fixLog,{'Fixed Types of variables': 27, 'Vars needing type fix': {'[Br]([C])': 'EnumVariable to FloatVariable'}})

        fixLog.clear()
        for ex in self.badVarType2Data:
                self.assert_(ExampleFix.fixExample(ex)!=None)
        self.assertEqual(fixLog,{'Fixed Types of variables': 30, 'Vars needing type fix': {'[O]([C])': 'EnumVariable to FloatVariable', '[Br]([C])': 'EnumVariable to FloatVariable', 'Level': 'EnumVariable to FloatVariable'}})

        fixLog.clear()
        ExampleFix = dataUtilities.ExFix(self.testData.domain,fixLog)
        for ex in self.badVarTypeSOMEData:
                self.assert_(ExampleFix.fixExample(ex)!=None)
        self.assertEqual(fixLog,{'Fixed Types of variables': 30, 'Vars needing type fix': {'[N]([N])': 'EnumVariable to FloatVariable', 'Measure': 'EnumVariable to FloatVariable'}})

        #test fix of example with missing value

        # Measure - Continuous attribute
        ex3=self.badVarTypeData[10]
        ex3["Measure"]="?"
        ex3Str=str(ex3)
        self.assertEqual(ex3Str[1],"?")
        ExampleFix = dataUtilities.ExFix(self.testData.domain,None,True)
        exFixed3 = ExampleFix.fixExample(ex3)
        self.assert_(ex3Str.replace("'","")==str(exFixed3).replace("'",""),
                        "Example fixed changed relative to the original!\n"+ex3Str+"\n"+str(exFixed3)+"\n")
        self.assert_(exFixed3.domain==self.testData.domain,"Fixed Domain 3 does not mathch to the requested domain")


        #DiscAttr1 -  Discrete attribute
        ex4=self.badVarTypeData[10]
        ex4["DiscAttr1"]="?"
        ex4Str=str(ex4)
        self.assertEqual(ex4Str[57],"?")
        exFixed4 = ExampleFix.fixExample(ex4)
        self.assert_(ex4Str.replace("'","")==str(exFixed4).replace("'",""),
                        "Example fixed changed relative to the original!\n"+ex4Str+"\n"+str(exFixed4)+"\n")
        self.assert_(exFixed4.domain==self.testData.domain,"Fixed Domain 4 does not mathch to the requested domain")

        
        #  [Br]([C]) - the attribute being changed type
        ex5=self.badVarTypeData[10]
        ex5["[Br]([C])"]="?"
        ex5Str=str(ex5)
        self.assertEqual(ex5Str[5],"?")
        exFixed5 = ExampleFix.fixExample(ex5)
        self.assert_(ex5Str.replace("'","")==str(exFixed5).replace("'",""),
                        "Example fixed changed relative to the original!\n"+ex5Str+"\n"+str(exFixed5)+"\n")
        self.assert_(exFixed5.domain==self.testData.domain,"Fixed Domain 5 does not mathch to the requested domain")

        
        #Test incompatible examples
        #print "\nTesting incompatible examples! The vars will be tet to '?'" 
        fixLog.clear()
        ExampleFix = dataUtilities.ExFix(self.testData.domain,fixLog,True)
        self.assert_(ExampleFix.fixExample(self.badVarNameData[0])==None,
                        "VarName: This example is not compatible, It should't be able to convert it!")
        self.assertEqual(fixLog,{'Missing Attributes': {'[O]([C])': 1}})
        ExampleFix.set_examplesFixedLog(None)
        self.assert_(ExampleFix.fixExample(self.badVarTypeData[0])[1]=='?',
                        "VarType: This example is not compatible, It should't be able to convert it!")
        ExampleFix.set_examplesFixedLog(fixLog)
        self.assert_(ExampleFix.fixExample(self.badVarCountData[0])==None,
                        "VarCount: This example is not compatible, It should't be able to convert it!")
        self.assertEqual(fixLog,{'Missing Attributes': {'YetOther': 1, '[O]([C])': 1}})

        #### Teting now real case where fixing example is needed
        from trainingMethods import AZorngPLS
        classifier=AZorngPLS.PLSLearner(self.testData)
        
        # Test different but compatible vartypes
        self.assert_(classifier(self.badVarTypeData[12])=="NEG","VarType: Prediction was not done correcly")
        #test fixed  Different data order  
        self.assert_(classifier(self.badVarOrderData[0])=="NEG","VarOrder(0): Prediction was not done correcly")
        self.assert_(classifier(self.badVarOrderData[12])=="NEG","VarOrder(12): Prediction was not done correcly")
        # Test the same badVarOrderDataSet but also with different order of values in 2 discrete attributes: 
        #    Activity[POS NEG] -> [NEG POS] and Attr3 [2 1 3 5 4] -> [1 2 3 4 5]
        #Predictions must be the same as before since just the order of values are chamged
        badVarOrderValuesData = dataUtilities.DataTable(self.badVarOrderValuesPath,createNewOn=orange.Variable.MakeStatus.OK)
        self.assert_(str(self.badVarOrderData.domain["Activity"].values)=="<POS, NEG>")
        self.assert_(str(badVarOrderValuesData.domain["Activity"].values)=="<NEG, POS>")
        self.assert_(str(self.badVarOrderData.domain["Attr3"].values)=="<1, 2, 3, 4, 5>")
        self.assert_(str(badVarOrderValuesData.domain["Attr3"].values)=="<2, 1, 3, 5, 4>")
        self.assert_(classifier(badVarOrderValuesData[12]).value=="NEG","VarOrderValues (12): Prediction was not done correcly")
        self.assert_(classifier(badVarOrderValuesData[1]).value=="POS","VarOrderValues (0): Prediction was not done correcly")
        example = badVarOrderValuesData[1]
        self.assert_(example["Activity"].value == self.badVarOrderData[1]["Activity"].value)
        self.assert_(int(example["Activity"]) != int(self.badVarOrderData[1]["Activity"]))
        self.assert_(example["Attr3"].value ==  self.badVarOrderData[1]["Attr3"].value)
        self.assert_(int(example["Attr3"]) !=  int(self.badVarOrderData[1]["Attr3"]))

        #print "\nTesting incompatible datasets: Although they have miss values '?' they should be predicted!"
        # Test incompatible different varNames
        self.assert_(classifier(self.badVarNameData[0])==None,"VarName: This prediction should NOT be possible")
        self.assertEqual(classifier.examplesFixedLog,{'Missing Attributes': {'[O]([C])': 1}, 'Fixed Types of variables': 1, 'Vars needing type fix': {'[Br]([C])': 'EnumVariable to FloatVariable'}})
        # Test different and incompatible vartypes
        self.assert_(classifier(self.badVarTypeData[0])=='NEG',"VarType: This prediction should be possible")
        # Test incompatible different varCount
        self.assert_(classifier(self.badVarCountData[0])==None,"VarCount: This prediction should NOT be possible")
        self.assertEqual(classifier.examplesFixedLog,{'Missing Attributes': {'YetOther': 1, '[O]([C])': 1}, 'Fixed Types of variables': 2, 'Vars needing type fix': {'[Br]([C])': "EnumVariable to FloatVariable (some impossible conversions. It was set to '?' for some examples.)"}})

        # Test different but compatible var number
        classifier=AZorngPLS.PLSLearner(self.badVarCountData)
        self.assert_(classifier(self.testData[0])=="POS","VarOrder: Prediction was not done correcly")


        """Test the convertion of an example to a specific domain if possible(no order check)
        also tests the return of None when the conversion is not possible"""
        domain = self.testData.domain
        origDomain = str(domain)
        exOK = self.badVarTypeData[10]  #This example is compatible with domain 
        exOK2 =  self.badVarOrderData[1]
        exOK3 =  self.wMetaData[1]
         
        exKO = self.badVarTypeData[1]   #This example is NOT compatible with domain 
        exKO2 = self.badVarNameData[1]        
        exKO3 = self.badVarCountData[1]
        
        #Check the conversion if different varTypes
        origExOK = str(exOK)
        ExampleFix = dataUtilities.ExFix(domain,None,True)
        newExOK = ExampleFix.fixExample(exOK)
        self.assert_(newExOK.domain==domain,"Domain was not converted properly")
        self.assert_(origDomain==str(domain),"The original domain was changed. It shouldn't")
        self.assert_(origExOK == str(exOK), "The original example was changed. It should't")
        self.assert_(len(newExOK)==len(exOK),"Some values of the original example are not in the converted example")
        for var in newExOK.domain:
            self.assert_(var.name in exOK.domain,"The original example domain did not have the varName"+var.name+" !!?? How come!!??")
            self.assert_(var.name in newExOK.domain,"The new example domain does not have the varName"+var.name)
            self.assert_(str(newExOK[var.name])==str(exOK[var.name]),"The value of var "+var.name+"was changed during the conversion "+str(exOK[var.name])+"->"+str(newExOK[var.name]))

        #Check the convertion of different varOrder
        origExOK2 = str(exOK2)        
        newExOK2 = ExampleFix.fixExample(exOK2)
        self.assert_(newExOK2.domain==domain,"Domain was not converted properly")
        self.assert_(origDomain==str(domain),"The original domain was changed. It shouldn't")
        self.assert_(origExOK2 == str(exOK2), "The original example was changed. It should't")
        self.assert_(len(newExOK2)==len(exOK2),"Some values of the original example are not in the converted example")
        for var in newExOK2.domain:
            self.assert_(var.name in exOK2.domain,"The original example domain did not have the varName"+var.name+" !!?? How come!!??")
            self.assert_(var.name in newExOK2.domain,"The new example domain does not have the varName"+var.name)
            self.assert_(str(newExOK2[var.name])==str(exOK2[var.name]),"The value of var "+var.name+"was changed during the conversion "+str(exOK2[var.name]) +"->"+str(newExOK2[var.name]))
        
        #check the convertion of different varMetas
        origExOK3 = str(exOK3)      
        newExOK3 = ExampleFix.fixExample(exOK3)
        self.assert_(newExOK3.domain==domain,"Domain was not converted properly")
        self.assert_(origDomain==str(domain),"The original domain was changed. It shouldn't")
        self.assert_(origExOK3 == str(exOK3), "The original example was changed. It should't")
        self.assert_(len(newExOK3)==len(exOK3),"Some values of the original example are not in the converted example")
        for var in newExOK3.domain:
            self.assert_(var.name in exOK3.domain,"The original example domain did not have the varName"+var.name+" !!?? How come!!??")
            self.assert_(var.name in newExOK3.domain,"The new example domain does not have the varName"+var.name)            
            self.assert_(str(newExOK3[var.name])==str(exOK3[var.name]),"The value of var "+var.name+"was changed during the conversion "+str(exOK3[var.name]) +"->"+str(newExOK3[var.name]))

        #check the convertion of example to a smaller domain where the example has at least all the domain variables
        domain2 = self.badVarCountData.domain   
        ExampleFix = dataUtilities.ExFix(domain2,None,True)
        origDomain2 = str(domain2)
        origExOK4 = str(exOK)
        newExOK4 = ExampleFix.fixExample(exOK)
        self.assert_(newExOK4.domain==domain2,"Domain was not converted properly")
        self.assert_(origDomain2==str(domain2),"The original domain was changed. It shouldn't")
        self.assert_(origExOK4 == str(exOK), "The original example was changed. It should't")
        for var in newExOK4.domain:
            self.assert_(var.name in exOK.domain,"The original example domain did not have the varName"+var.name+" !!?? How come!!??")
            self.assert_(var.name in newExOK4.domain,"The new example domain does not have the varName"+var.name)
            self.assert_(str(newExOK4[var.name])==str(exOK[var.name]),"The value of var "+var.name+"was changed during the conversion "+str(exOK[var.name]) +"->"+str(newExOK4[var.name]))


        #print "\nTesting impossible conversions which will be set to '?'"
        #Check the return of None to impossible conversions:
        ExampleFix = dataUtilities.ExFix(domain,None,True)
        self.assert_(ExampleFix.fixExample(exKO)[1]=='?',"BadVarType: It shouldn't be able to do this convertion!")
        self.assert_(ExampleFix.fixExample(exKO2)==None,"BadVarName: It shouldn't be able to do this convertion!")
        self.assert_(ExampleFix.fixExample(exKO3)==None,"BadVarCount: It shouldn't be able to do this convertion!")
    def test_rmAllExMeta(self):
        """Test the remove of all meta attributes in one single example"""
        ex = self.wMetaData[0] 
        nVars = len(ex.domain)                 #Store the number and the attributes+classes, they must be the same at the end
        attributes = str(ex)[0:str(ex).find(', {')]
        self.assert_(str(ex).find('{') or str(ex).find('}') >=0,"The initial example had no meta group at all")
        self.assert_(len(ex.getmetas())>=1,"The initial example had no metas at all")
        dataUtilities.rmAllMeta(ex)
        self.assert_(len(ex.domain) == nVars,"Different number of attributes after removing metas")       #Number of attributes did not change
        self.assert_(str(ex)==attributes,"Attributes were changed during removal of metas")    #Attributes remain the same
        self.assert_(len(ex.getmetas()) == 0,"Metas were not removed")            #No meta found
        self.assert_(str(ex).find('{') or str(ex).find('}') ==-1,"Meta group was not removed") #no meta group found

    def test_rmAllDataMeta(self):
        """Test the remove of all meta attributes from an entire dataset"""
        nVars = len(self.wMetaData.domain)                 #Store the number and the attributes, they must be the same at the end
        domain = str(self.wMetaData.domain)[0:str(self.wMetaData.domain).find(', {')]
        self.assert_(str(self.wMetaData.domain).find('{') or str(self.wMetaData.domain).find('}') >=0,"The initial data had no meta group at all")
        self.assert_(len(self.wMetaData.domain.getmetas())>=1,"The initial data had no metas at all")

        dataUtilities.rmAllMeta(self.wMetaData)

        self.assert_(len(self.wMetaData.domain) == nVars,"Different number of attributes after removing metas")       #Number of attributes did not change
        self.assert_(str(self.wMetaData.domain)==domain,"Attributes were changed during removal of metas")    #Attributes remain the same
        self.assert_(len(self.wMetaData.domain.getmetas()) == 0,"Metas were not removed")            #No meta found
        self.assert_(str(self.wMetaData.domain).find('{') or str(self.wMetaData.domain).find('}') ==-1,"Meta group was not removed") #no meta group found
        
    def test_getExCopyWithoutMeta(self):
        """Test the return of a copy of one example but without meta attributes"""
        ex = self.wMetaData[0]
        nVars = len(ex.domain)                 #Store the number and the attributes, they must be the same at the end
        origAttributes = str(ex)
        attributes = str(ex)[0:str(ex).find(', {')]
        self.assert_(str(ex).find('{') or str(ex).find('}') >=0,"The initial example had no meta group at all")
        self.assert_(len(ex.getmetas())>=1,"The initial example had no metas at all")
        
        exNoMeta = dataUtilities.getCopyWithoutMeta(ex)

        #The initial example shall not be modified by this method
        self.assert_(len(ex.domain) == nVars,"nVars: The initial example was modified, it shouldn't") 
        self.assert_(str(ex)==origAttributes,"Attributes: The initial example was modified, it shouldn't")
        self.assert_(str(ex).find('{') or str(ex).find('}') >=0,"MetaGroup: The initial example was modified, it shouldn't")
        self.assert_(len(ex.getmetas())>=1,"Metas: The initial example was modified, it shouldn't")
        #The new example should have no metas now
        self.assert_(len(exNoMeta.domain) == nVars,"Different number of attributes after removing metas")       #Number of attributes did not change
        self.assert_(str(exNoMeta)==attributes,"Attributes were changed during removal of metas")    #Attributes remain the same
        self.assert_(len(exNoMeta.getmetas()) == 0,"Metas were not removed")            #No meta found
        self.assert_(str(exNoMeta).find('{') or str(exNoMeta).find('}') ==-1,"Meta group was not removed") #no meta group found


    def test_getDataCopyWithoutMeta(self):
        """Test the return of a copy of a dataset but without meta attributes"""
        nVars = len(self.wMetaData.domain)                 #Store the number and the attributes, they must be the same at the end
        origAttributes = str(self.wMetaData.domain)
        attributes = str(self.wMetaData.domain)[0:str(self.wMetaData.domain).find(', {')]
        self.assert_(str(self.wMetaData.domain).find('{') or str(self.wMetaData.domain).find('}') >=0,"The initial data had no meta group at all")
        self.assert_(len(self.wMetaData.domain.getmetas())>=1,"The initial data had no metas at all")
        
        noMetaData = dataUtilities.getCopyWithoutMeta(self.wMetaData)

        #The initial example shall not be modified by this method
        self.assert_(len(self.wMetaData.domain) == nVars,"nVars: The initial data was modified, it shouldn't")
        self.assert_(str(self.wMetaData.domain)==origAttributes,"Attributes: The initial data domain was modified, it shouldn't")
        self.assert_(str(self.wMetaData.domain).find('{') or str(self.wMetaData.domain).find('}') >=0,"MetaGroup: The initial data domain was modified, it shouldn't")
        self.assert_(len(self.wMetaData.domain.getmetas())>=1,"Metas: The initial data domain was modified, it shouldn't")
        #The new example should have no metas now
        self.assert_(len(noMetaData.domain) == nVars,"Different number of attributes after removing metas")       #Number of attributes did not change
        self.assert_(str(noMetaData.domain)==attributes,"Attributes were changed during removal of metas")    #Attributes remain the same
        self.assert_(len(noMetaData.domain.getmetas()) == 0,"Metas were not removed")            #No meta found
        self.assert_(str(noMetaData.domain).find('{') or str(noMetaData.domain).find('}') ==-1,"Meta group was not removed") #no meta group found


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(dataUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

