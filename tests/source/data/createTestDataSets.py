import string 
from AZutilities import dataUtilities
import random
import orange
import os

dataFile = "baseData.txt"

fileH = open(dataFile,"r")
lines = fileH.readlines()
fileH.close()

dataFile = "baseDataTmp.txt"
fileH = open(dataFile,"w")
for idx,line in enumerate(lines):
    if idx == 0:
        line = line.replace("[C]([C]=[C])","Measure").replace("activity","Activity")
        line = line.replace("[C](=[C][N][N])\t[C](=[C][N][O])\t[C](=[C][N][S])\t[C](=[C][O])","DiscAttr1\tDiscAttr2\tAttr3\tYetOther")
    if idx == 1:
        line = line.replace("C1=CC(=CC=C1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]\t5959\t2\t0\t0\t0\t2\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0",\
                            "C1=CC(=CC=C1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]\t5959\t2\t0\t0\t0\t2\t0.0\t0.0\t0.0\t0.0\tRed\tYES\tA\tB")
    fileH.write(line)
fileH.close()

data = dataUtilities.DataTable(dataFile,noMeta = False)
dataFile = "baseData.tab"

#[Br]([C])       [C](=[C])       [C](=[C][Br][I])        [C](=[C][Cl])   DiscAttr1       DiscAttr2       Attr3           YetOther        
# continuous      continuous      continuous             continuous      Red Green Blue  YES NO          1 2 3 4 5       A B C 1 2      
#     0                1                2                      3             4                5             6                  7

data.domain["DiscAttr1"].values.append("Green")
data.domain["DiscAttr1"].values.append("Blue")
data.domain["DiscAttr2"].values.append("NO")
for val in ["1","2","3","4","5"]:
    data.domain["Attr3"].values.append(val)
for val in ["A","B","C","1","2"]:
    data.domain["YetOther"].values.append(val)
data.domain["Activity"].values.append("POS")
data.domain["Activity"].values.append("NEG")
data.domain["[O](=[C])"].name = "Level"
data.save(dataFile)
data = dataUtilities.DataTable(dataFile)

NEX = len(data)

for idx,ex in enumerate(data):
  if ex["Activity"].value == "1":
        ex["Activity"] = "POS"
  else:
        ex["Activity"] = "NEG"
  if idx%2 == 0:
    ex[4] = "Red"
  elif idx%3 == 0:
    ex[4] = "Green"
  else:
    ex[4] = "Blue"

for idx,ex in enumerate(data):
    if idx%2:
        ex[5] = "YES"
    else:    
        ex[5] = "NO"



idx = 1
for ex in data:
    ex[6] = str(idx)
    idx = idx+1
    if idx == 6:
       idx = 1

idx = 0
values = ["A","B","C","1","2"]
for ex in data:
    ex[7] = values[idx]
    idx = idx+1
    if idx == 4:
       idx = 0

data = dataUtilities.getDataWithoutUnusedValues(data,True)

random.seed(2)
for ex in data:
    ex["Level"] = ex["Level"] + random.random()



varsS = ["Measure","[Br]([C])","[N]([N])","[O]([C])","[C]([C][F])","Level","DiscAttr1","DiscAttr2","Attr3","YetOther","Activity"]
domain = orange.Domain([data.domain[attr] for attr in varsS])
domain.addmetas(data.domain.getmetas())
data = dataUtilities.DataTable(domain,data)
#LargeDataset
#data.save("")

#============ Classification ====================
# With Meta
fileName = "BinClass"+"_W_metas_"
data.save(fileName + "Train.tab")
os.system("head -n "+str(int(NEX*0.8))+" "+fileName + "Train.tab > "+fileName + "Test.tab")

# No Meta
domain = orange.Domain(data.domain.attributes,data.domain.classVar)
data2 = dataUtilities.DataTable(domain,data)
fileName = "BinClass"+"_No_metas_"
data2.save(fileName + "Train.tab")
os.system("head -n "+str(int(NEX*0.8))+" "+fileName + "Train.tab > "+fileName + "Test.tab")

data2[0][6] = '?'
data2[4][0] = '?'
data2.save(fileName + "Train_missing.tab")

#========BinClass_No_metas_UnusedValues_Train.tab=======
fileName = "BinClass_No_metas_Train.tab"
data3 = dataUtilities.DataTable(fileName)
data3.domain["Activity"].values.append("OTHER")
for val in ["6","7","8","9"]:
    data3.domain["Attr3"].values.append(val)
data3.save("BinClass_No_metas_UnusedValues_Train.tab")



#============ Regression ====================
domainR = orange.Domain([attr for attr in data.domain if attr.name != "Measure" and attr.name != "Activity"],data.domain["Measure"])
domainR.addmetas(data.domain.getmetas())
dataR = dataUtilities.DataTable(domainR,data)
random.seed(1)
for ex in dataR:
    ex["Measure"] = ex["Measure"] + random.random()
# With Meta
fileName = "Reg"+"_W_metas_"
dataR.save(fileName + "Train.tab")
os.system("head -n "+str(int(NEX*0.8))+" "+fileName + "Train.tab > "+fileName + "Test.tab")

# No Meta
domain = orange.Domain(dataR.domain.attributes,dataR.domain.classVar)
data2 = dataUtilities.DataTable(domain,dataR)
fileName = "Reg"+"_No_metas_"
data2.save(fileName + "Train.tab")
os.system("head -n "+str(int(NEX*0.8))+" "+fileName + "Train.tab > "+fileName + "Test.tab")





#================== Full Numeric ======================
DataDesc = "BinClass"
data = dataUtilities.DataTable("baseData.txt",noMeta = False)
data.domain["activity"].values.append("POS")
data.domain["activity"].values.append("NEG")
for idx,ex in enumerate(data):
  if ex["activity"].value == "1":
        ex["activity"] = "POS"
  else:
        ex["activity"] = "NEG"

NEx = len(data)

data = dataUtilities.getDataWithoutUnusedValues(data,True)
# With Meta
fileName = DataDesc+"_W_metas_FullNumeric_"
data.save(fileName + "Train.tab")
os.system("head -n "+str(int(NEx*0.8))+" "+fileName + "Train.tab > "+fileName + "Test.tab")

# No Meta
domain = orange.Domain(data.domain.attributes,data.domain.classVar)
data = dataUtilities.DataTable(domain,data)
fileName = DataDesc+"_No_metas_FullNumeric_"
data.save(fileName + "Train.tab")
os.system("head -n "+str(int(NEx*0.8))+" "+fileName + "Train.tab > "+fileName + "Test.tab")
#---- Regression ---
DataDesc = "Reg"
data = dataUtilities.DataTable("baseData.txt",noMeta = False)
data.domain["[C]([C]=[C])"].name = "Measure"
domain = orange.Domain([attr for attr in data.domain.attributes if attr.name not in ["activity","Measure"]],data.domain["Measure"])
data = orange.ExampleTable(domain,data) 
random.seed(6)
for ex in data:
    ex["Measure"] = ex["Measure"] + random.random()
data.save(DataDesc+"_No_metas_FullNumeric_Train.tab")

#================ Create small test set =================
dataFile = "BinClass_No_metas_SmallTest.tab"
os.system("head -n 33 BinClass_No_metas_Test.tab > "+dataFile)
data = dataUtilities.DataTable(dataFile)
var = orange.StringVariable("Comments")
data.domain.addmeta(-1,var)
idxs=[2,5,12,25]
for idx,ex in enumerate(data):
    if idx in idxs:
        ex["Comments"] = "notok"
    else:
        ex["Comments"] = "ok"
data.save("BinClass_W_metas_SmallTest.tab")
#================ Test data for BAD thinks :) =================
fileH = open(dataFile,"r")
lines = fileH.readlines()
fileH.close()
#-----------------
saveFile = "BinClass_BadVarType.tab"
fileH = open(saveFile,"w")
for idx,line in enumerate(lines):
    line = line.split("\t")
    if idx == 1:
        line[1] = "discrete"
    elif idx == 3:
        line[1] = "unknown_symbol"
    elif idx == 4:
        line[1] = "bad_error"
    elif idx == 5:
        line[1] = "something"
    line = string.join(line,"\t")
    fileH.write(line)
fileH.close()
#-----------------
saveFile = "BinClass_BadVarName.tab"
data = dataUtilities.DataTable(dataFile)
data.domain[3].name = "OtherName"
data.save(saveFile)
#-----------------
saveFile = "BinClass_BadVarOrder.tab"
data = dataUtilities.DataTable(dataFile)
varsN = ["Measure","[Br]([C])","[N]([N])","[O]([C])","[C]([C][F])","Level","DiscAttr1","DiscAttr2","YetOther","Attr3","Activity"]
data= dataUtilities.DataTable(orange.Domain([data.domain[attr] for attr in varsN]),data)
data.save(saveFile)
#-----------------
saveFile = "BinClass_BadVarOrderValues.tab"
fileH = open(saveFile,"w")
for idx,line in enumerate(lines):
    if idx == 1:
        line = line.replace("1 2 3 4 5","2 1 3 5 4")
        line = line.replace("POS NEG","NEG POS")
    fileH.write(line)
fileH.close()
#-----------------
saveFile = "BinClass_BadVarCount.tab"
data = dataUtilities.DataTable(dataFile)
varsN = ["Measure","[Br]([C])","[N]([N])","[O]([C])","[C]([C][F])","Level","DiscAttr1","DiscAttr2","Attr3","Activity"]
data= dataUtilities.DataTable(orange.Domain([data.domain[attr] for attr in varsN]),data)
data.save(saveFile)
#-----------------
saveFile = "BinClass_BadVarType2.tab"
fileH = open(saveFile,"w")
for idx,line in enumerate(lines):
    line = line.split("\t")
    if idx == 1:
        line[1] = "discrete"
        line[3] = "discrete"
        line[5] = "discrete"
    line = string.join(line,"\t")
    fileH.write(line)
fileH.close()
#-----------------
BarVarOrder = "BinClass_BadVarOrder.tab"
fileH = open(BarVarOrder,"r")
linesBVO = fileH.readlines()
fileH.close()

saveFile = "BinClass_BadVarTypeSOME.tab"
fileH = open(saveFile,"w")
for idx,line in enumerate(linesBVO):
    line = line.split("\t")
    if idx == 1:
        line[0] = "discrete"
        line[2] = "discrete"
    line = string.join(line,"\t")
    fileH.write(line)
fileH.close()

#===================== Reg Data for testing Imputation  =================
data = dataUtilities.DataTable("Reg_No_metas_Train.tab")
data.domain["DiscAttr2"].values = ["NO","YES"]
for ex in data:
    if ex["Level"].value > 2:
        ex["DiscAttr2"] = "YES" 
    else:
        ex["DiscAttr2"] = "NO"
train = orange.ExampleTable(data[0:int(len(data)/5)])
test = orange.ExampleTable(data[int(len(data)/5)+1:])
train.save("Reg_No_metas_Imp_Train.tab")
test.save("Reg_No_metas_Imp_Test.tab")

os.system("rm baseDataTmp.txt")
os.system("rm baseData.tab")
