"""
<name>Test Optimized Learners</name>
<description>Tests the accuracy of optimized Learners on a data set.</description>
<icon>icons/TestOptLearners.png</icon>
<contact>Pedro Almeida</contact>
<priority>200</priority>
"""
#
# OWAZTestOptLearners.py
#
from OWWidget import *
import OWGUI
from AZutilities import  getAccWOptParam
import time
import warnings
import AZLearnersParamsConfig
warnings.filterwarnings("ignore", "'id' is not a builtin attribute",
                        orange.AttributeWarning)

##############################################################################
#         Return s a confusion matrix in the form of a vector:
#                For Binary classifiers:  
#                        [[TP, FN],
#                         [FP, TN]]
#                For classifiers with class having N values:
#                                             Predicted class
#                                        |   A       B       C
#                                     ---------------------------
#                           known     A  |  tpA     eAB     eAC
#                           class     B  |  eBA     tpB     eBC
#                                     C  |  eCA     eCB     tpC
#
#                    [[tpA, eAB, ..., eAN],
#                     [eBA, tpB, ..., eBN],
#                      ...,
#                     [eNA, eNB, ..., tpN ]]
#
#                 where A, B, C are the class values in the same order as testData.domain.classVar.values
#
#  from http://www.compumine.com/web/public/newsletter/20071/precision-recall
#       Sensitivity = tp/(tp+fn)
#       SensitivityA = tpA/(tpA+eAB+eAC)
#
#       Specificity = tn/(tn+fp)
#       SpecificityA = tnA/(tnA+eBA+eCA), where tnA = tpB + eBC + eCB + tpC 

def sens(cm,idx):
    TP = float(cm[idx][idx])
    FN = float(sum(cm[idx]) - TP)
    return TP/(TP+FN)
def spec(cm,idx):
    TN = 0
    FP = 0
    for row in [x for (i,x) in enumerate(cm) if i != idx]:
            FP += row[idx] 
            TN += sum(row)-row[idx]

    return TN/(TN+FP)

class Learner:
    def __init__(self, learner, id):
        learner.id = id
        self.learner = learner
        self.name = learner.name
        self.id = id
        self.scores = []
        self.results = None
        self.evaluator = None
        self.time = time.time() # used to order the learners in the table

#return:  {"RMSE":0.2,"R2":0.1,"CA":0.98,"CM":[[TP, FP],[FN,TN]]}
class Score:
    def __init__(self, name, label, res, show=True, cmBased=False):
        self.name = name
        self.label = label
        self.res = res        # if is cmBased, then res represent a function to be applyed to the CM
        self.show = show
        self.cmBased = cmBased

class OWAZTestOptLearners(OWWidget):
    settingsList = ["nInnerFolds", "nOuterFolds", "precision",
                    "selectedCScores", "selectedRScores", "applyOnAnyChange"]
    contextHandlers = {"": DomainContextHandler("", ["targetClass"])}
    callbackDeposit = []


    cStatistics = [Score(*s) for s in [\
        ('Classification accuracy', 'CA', "CA", True),
        ('Sensitivity', 'Sens', 'sens(CM,classIndex)', True, True),
        ('Specificity', 'Spec', 'spec(CM,classIndex)', True, True),
        ]]

    rStatistics = [Score(*s) for s in [\
        ("Root mean squared error", "RMSE", "RMSE",True),
        ("R-squared", "R2", "R2",True)
        ]]


    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "TestOptLearners")

        self.inputs = [("Data", ExampleTable, self.setData, Default), ("Learner", orange.Learner, self.setLearner, Default)]
        #Output the conf. Matrix if exists
        #self.outputs = [("Evaluation Results", orngTest.ExperimentResults)]

        # Settings
        self.learnerType = None
        self.nInnerFolds = 5            # Inner folds
        self.nOuterFolds = 5            # Outer folds
        self.precision = 4
        self.applyOnAnyChange = False
        self.selectedCScores = [i for (i,s) in enumerate(self.cStatistics) if s.show]
        self.SelectedParams = []
        self.paramsLabels = []
        self.selectedRScores = [i for (i,s) in enumerate(self.rStatistics) if s.show]
        self.targetClass = 0
        self.loadSettings()

        self.stat = self.cStatistics

        self.data = None                # input data set
        self.learners = {}              # set of learners (input)

        # GUI
        box = OWGUI.widgetBox(self.controlArea, "Sampling", orientation='vertical', addSpace=False)
        OWGUI.spin(box, self, 'nInnerFolds', 2, 100, step=1, label='Number of inner folds:',
                   callback=lambda p=0: self.conditionalRecompute(p))
        OWGUI.spin(box, self, 'nOuterFolds', 2, 100, step=1, label='Number of outer folds:',
                   callback=lambda p=0: self.conditionalRecompute(p))

        self.OptParam = OWGUI.listBox(self.controlArea, self, 'SelectedParams',
                                     'paramsLabels', box = "Parameters optimization",
                                     selectionMode = QListWidget.MultiSelection,
                                     callback=self.newParamSelection)

        #OWGUI.separator(self.sBtns, height = 3)
        
        
        boxA= OWGUI.widgetBox(box, orientation='vertical', addSpace=False)
        OWGUI.separator(boxA)
        OWGUI.checkBox(boxA, self, 'applyOnAnyChange',
                       label="Apply on any change", callback=self.applyChange)
        self.applyBtn = OWGUI.button(box, self, "&Apply",
                                     callback=lambda f=True: self.recompute(f))
        self.applyBtn.setDisabled(True)


        OWGUI.separator(self.controlArea)

        # statistics
        self.cbox = OWGUI.widgetBox(self.controlArea)
        self.cStatLabels = [s.name for s in self.cStatistics]
        self.cstatLB = OWGUI.listBox(self.cbox, self, 'selectedCScores',
                                     'cStatLabels', box = "Performance scores",
                                     selectionMode = QListWidget.MultiSelection,
                                     callback=self.newscoreselection)
        OWGUI.separator(self.cbox)
        self.targetCombo=OWGUI.comboBox(self.cbox, self, "targetClass", orientation=0,
                                        callback=[self.changedTarget],
                                        box="Target class")

        self.rStatLabels = [s.name for s in self.rStatistics]
        self.rstatLB = OWGUI.listBox(self.controlArea, self, 'selectedRScores',
                                     'rStatLabels', box = "Performance scores",
                                     selectionMode = QListWidget.MultiSelection,
                                     callback=self.newscoreselection)
        self.rstatLB.box.hide()

        ##scPA
        self.filePath=os.path.join(os.getcwd(),"TOLResults.txt")
        BoxF = OWGUI.widgetBox(self.controlArea,"File for saving results")
        # Set location of model file 
        OWGUI.lineEdit(BoxF, self, 'filePath',tooltip=\
                'Enter full file path to save the displayed results to a tab separated file.')

        button = OWGUI.button(BoxF, self, 'Browse...', callback = self.browseFile, disabled=0,tooltip = "Browse for a location...")
        #button.setMaximumWidth(25)
        # Save the Results
        OWGUI.button(self.controlArea, self,"&Save Results", callback=self.saveRes)
        ##ecPA


        # score table
        # table with results
        #self.g = OWGUI.widgetBox(self.mainArea, 'Evaluation Results')
        self.tab = OWGUI.table(self.mainArea, selectionMode = QTableWidget.NoSelection)

        # Confusion Matrix Table
        #self.CM = OWGUI.widgetBox(self.mainArea, 'Confusion Matrix')
        self.CMtab = OWGUI.table(self.mainArea, selectionMode = QTableWidget.NoSelection)


        import sip   # (row, col)
        sip.delete(self.mainArea.layout())
        self.layout = QGridLayout(self.mainArea)
        self.layout.addWidget(OWGUI.widgetLabel(self.mainArea, "Evaluation Results"), 0, 0, Qt.AlignLeft)        
        self.layout.addWidget(self.tab, 1, 1)
        self.layout.addWidget(OWGUI.widgetLabel(self.mainArea, "Confusion Matrix"), 2, 0, Qt.AlignLeft)        
        self.layout.addWidget(OWGUI.widgetLabel(self.mainArea, "Prediction"), 3, 1, Qt.AlignCenter)
        self.layout.addWidget(OWGUI.widgetLabel(self.mainArea, "Correct <br>Class  "), 4, 0, Qt.AlignRight)
        self.layout.addWidget(self.CMtab, 4, 1)
        #self.layout.setColumnStretch(1, 100)
        #self.layout.setRowStretch(2, 100)

        #self.resize(680,470)

    ##scPA
    def browseFile(self):
        # Possible modes:
        var = os.path.realpath(str(self.filePath))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        else:
            startfile=os.getcwd()

        filename = QFileDialog.getSaveFileName(self, "Save Results", startfile)
        if filename:
            self.filePath = str(filename)

    def saveRes(self):
        if self.tab.horizontalHeader().count() and self.tab.rowCount()>0 and self.tab.columnCount()>0:
                if os.path.isdir(os.path.split(self.filePath)[0]) and self.filePath!="" and self.filePath!=None:
                        try:
                            line=""
                            # Open the file
                            f = open(self.filePath, 'w')
                            # Write the Header of the table 
                            for nCol in range(self.tab.columnCount()):
                                line += str(self.tab.horizontalHeaderItem(nCol).text()) + "\t"
                            f.write(line[0:-1]+"\r\n")
                            # Write all lines of the table
                            for nRow in range(self.tab.rowCount()):
                                line = ""
                                for nCol in range(self.tab.columnCount()):
                                    line += str(self.tab.item(nRow, nCol).text()) + "\t"
                                f.write(line[0:-1]+"\r\n")
                            #Save the Confusion Matrix
                            if self.isclassification():
                                f.write("\r\nConfusion Matrix\r\n")
                                # Write the Header of the table 
                                line = "\t"
                                for nCol in range(self.CMtab.columnCount()):
                                    line += str(self.CMtab.horizontalHeaderItem(nCol).text()) + "\t"
                                f.write(line[0:-1]+"\r\n")
                                # Write all lines of the table
                                for nRow in range(self.CMtab.rowCount()):
                                    line = str(self.CMtab.verticalHeaderItem(nRow).text()) + "\t"
                                    for nCol in range(self.CMtab.columnCount()):
                                        line += str(self.CMtab.item(nRow, nCol).text()) + "\t"
                                    f.write(line[0:-1]+"\r\n")
                            # Close the File
                            f.close()
                            print "Results saved to: " + self.filePath
                        except:
                            type, val, traceback = sys.exc_info()
                            sys.excepthook(type, val, traceback)  # print the exception
                            self.error ("ERROR: Failed to save the results.")
                else:
                        print "Failed to save the results: Invalid  path."
                        self.error ("Failed to save the results: Invalid  path.")
        else:
                print "There are no results yet to save!"
                self.error("There are no results yet to save!")



    # scoring and painting of score table
    def isclassification(self):
        if not self.data or not self.data.domain.classVar:
            return True
        return self.data.domain.classVar.varType == orange.VarTypes.Discrete
        
    def paintscores(self):
        """paints the table with evaluation scores"""

        self.tab.setColumnCount(len(self.stat)+1)
        self.tab.setHorizontalHeaderLabels(["Method"] + [s.label for s in self.stat])
        
        prec="%%.%df" % self.precision

        learners = [(l.time, l) for l in self.learners.values()]
        learners.sort()
        learners = [lt[1] for lt in learners]

        self.tab.setRowCount(len(self.learners))
        for (i, l) in enumerate(learners):
            OWGUI.tableItem(self.tab, i,0, l.name)
            
        for (i, l) in enumerate(learners):
            if l.scores:
                for j in range(len(self.stat)):
                    if l.scores[j] is not None:
                        OWGUI.tableItem(self.tab, i, j+1, prec % l.scores[j])
                    else:
                        OWGUI.tableItem(self.tab, i, j+1, "N/A")
            else:
                for j in range(len(self.stat)):
                    OWGUI.tableItem(self.tab, i, j+1, "")
        
        # adjust the width of the score table cloumns
        self.tab.resizeColumnsToContents()
        self.tab.resizeRowsToContents()

        #Call the method for painting the Confusion Matrix
        self.paintConfMat()

        usestat = [self.selectedRScores, self.selectedCScores][self.isclassification()]
        for i in range(len(self.stat)):
            if i not in usestat:
                self.tab.hideColumn(i+1)

    def paintConfMat(self):
        if not self.learners or not self.data or not self.data.domain.classVar or not self.isclassification():
            self.CMtab.setColumnCount(0)
            self.CMtab.setRowCount(0)
            return
        #Use only Confusion M;Atrix CM from the first and unique learner 
        learners = [(l.time, l) for l in self.learners.values()]
        learners.sort()
        learners = [lt[1] for lt in learners]
        CM = learners[0].results["CM"]
        nClasses = len(self.data.domain.classVar.values)

        #Create CM Vertical and Horizontal Headers
        self.CMtab.setColumnCount(nClasses + 1)
        self.CMtab.setHorizontalHeaderLabels([str(x) for x in self.data.domain.classVar.values]+[" "])
        self.CMtab.setRowCount(nClasses + 1)
        self.CMtab.setVerticalHeaderLabels([str(x) for x in self.data.domain.classVar.values]+[" "])

        #Fill CM data and last Col Sums
        for (r, row) in enumerate(CM):
            for (c, val) in enumerate(row):
                OWGUI.tableItem(self.CMtab, r, c, val) 
            OWGUI.tableItem(self.CMtab, r, c+1, sum(row))
        #last row with sums
        for idx in range(nClasses): 
            OWGUI.tableItem(self.CMtab, r+1, idx, sum([x[idx] for x in CM ]))
        # Sum od Sums
        OWGUI.tableItem(self.CMtab, r+1, idx+1, sum([sum(x) for x in CM]))

        boldf = self.CMtab.item(0, nClasses).font()
        boldf.setBold(True)
        for ri in range(nClasses+1):
            self.CMtab.item(ri, nClasses).setFont(boldf)
            self.CMtab.item(nClasses, ri).setFont(boldf)


        self.CMtab.resizeColumnsToContents()
        self.CMtab.resizeRowsToContents()


    def sendReport(self):
        exset = [("Number of Inner Folds", self.nInnerFolds), ("Number of Outer Folds",self.nOuterFolds)]
        
        self.reportData(self.data)
        
        self.reportSection("Results")
        learners = [(l.time, l) for l in self.learners.values()]
        learners.sort()
        learners = [lt[1] for lt in learners]
        usestat = [self.selectedRScores, self.selectedCScores][self.isclassification()]
        
        res = "<table><tr><th></th>"+"".join("<th><b>%s</b></th>" % hr for hr in [s.label for i, s in enumerate(self.stat) if i in usestat])+"</tr>"
        for i, l in enumerate(learners):
            res += "<tr><th><b>%s</b></th>" % l.name
            if l.scores:
                for j in usestat:
                    scr = l.scores[j]
                    res += "<td>" + ("%.4f" % scr if scr is not None else "") + "</td>"
            res += "</tr>"
        res += "</table>"
        self.reportRaw(res)

        if self.learners and self.data and self.data.domain.classVar and self.isclassification():
            self.reportSection("Confusion Matrix")
            classVals = [str(x) for x in self.data.domain.classVar.values]
            nClassVals = len(classVals)
            res = "<table>\n<tr><td></td>" + "".join('<td align="center"><b>&nbsp;&nbsp;%s&nbsp;&nbsp;</b></td>' % cv for cv in classVals) + "</tr>\n"
            for i, cv in enumerate(classVals):
                res += '<tr><th align="right"><b>%s</b></th>' % cv + \
                   "".join('<td align="center">%s</td>' % self.CMtab.item(i, j).text() for j in range(nClassVals)) + \
                   '<th align="right"><b>%s</b></th>' % self.CMtab.item(i, nClassVals).text() + \
                   "</tr>\n"
            res += '<tr><th></th>' + \
               "".join('<td align="center"><b>%s</b></td>' % self.CMtab.item(nClassVals, j).text() for j in range(nClassVals+1)) + \
               "</tr>\n"
            res += "</table>\n<p><b>Note:</b> columns represent predictions, row represent true classes</p>"
            self.reportRaw(res)

            
    def score(self, ids):
        """compute scores for the list of learners"""
        if (not self.data):
            for id in ids:
                self.learners[id].results = None
            return

        pb = OWGUI.ProgressBar(self, iterations=1+len(self.learners))
        # test which learners can accept the given data set
        # e.g., regressions can't deal with classification data
        learners = []
        n = len(self.data.domain.attributes)*2
        ##scPA     Changed from "p0min(n, len(self.data))="  which could result in empty dataset , to  "p0=len(self.data) - min(n, len(self.data))"       ##ecPA
        indices = orange.MakeRandomIndices2(p0=len(self.data) - min(n, len(self.data)), stratified=orange.MakeRandomIndices2.StratifiedIfPossible)
        new = self.data.selectref(indices(self.data))
        #new = self.data.selectref([1]*min(n, len(self.data)) +
        #                          [0]*(len(self.data) - min(n, len(self.data))))
        self.warning(0)
        for l in [self.learners[id] for id in ids]:
            try:
                predictor = l.learner(new)
                if predictor(new[0]).varType == new.domain.classVar.varType:
                    learners.append(l)
                else:
                    l.scores = []
            except Exception, ex:
                self.warning(0, "Learner %s ends with exception: %s" % (l.name, str(ex)))
                l.scores = []
        pb.advance() 
        if not learners:
            pb.finish()
            return

        # computation of results (res, and cm if classification)
        for l in learners:
            #Set the selected parameters to optimize
            paramList = []
            for paramIdx in self.SelectedParams:
                paramList.append(self.paramsLabels[paramIdx])
            #print "getAccWOptParam command: "
            #for x in (self.data,l.learner, paramList, self.nOuterFolds,self.nInnerFolds):
            #    print x
            l.evaluator = getAccWOptParam.AccWOptParamGetter(data = self.data, learner = l.learner, paramList = paramList, nExtFolds = self.nOuterFolds, nInnerFolds = self.nInnerFolds)
            l.results = l.evaluator.getAcc()
            pb.advance()

        #if self.isclassification():
        #    cm = orngStat.computeConfusionMatrices(res, classIndex = self.targetClass)

        #self.error(range(len(self.stat)))
                
        for (i, l) in enumerate(learners):
            scores = []
            for ii, s in enumerate(self.stat):
                try:
                    if s.cmBased:
                        classIndex = self.targetClass
                        CM = l.results["CM"]
                        scores.append(eval(s.res))
                    else:
                        scores.append(l.results[s.res])
                except Exception, ex:
                    self.error(ii, "An error occurred while evaluating " + s.res + " on %s due to %s" % \
                              (l.name, ex.message))
                    scores.append(None)
            self.learners[l.id].scores = [s if s else None for s in scores]
        pb.finish()    
        self.sendResults()

    def recomputeCM(self):
        if not self.learners or not self.stat:
            return
        scores = [(indx, s.res)
                  for (indx, s) in enumerate(self.stat) if s.cmBased]
        for (indx, score) in scores:
            for (i, l) in enumerate([l for l in self.learners.values() if l.scores]):
                classIndex = self.targetClass
                CM = l.results["CM"]
                l.scores[indx] = eval(score) 
        self.paintscores()
        

    # handle input signals

    def setData(self, data):
        """handle input train data set"""
        self.closeContext()
        self.data = data
        self.fillClassCombo()
        if not self.data:
            # data was removed, remove the scores
            for l in self.learners.values():
                l.scores = []
                l.results = None
            #self.send("Evaluation Results", None)
        else:
            # new data has arrived
            self.data = orange.Filter_hasClassValue(self.data)
            if self.isclassification():
                self.rstatLB.box.hide()
                self.cbox.show()
            else:
                self.cbox.hide()
                self.rstatLB.box.show()
            self.stat = [self.rStatistics, self.cStatistics][self.isclassification()]
            
            if self.learners:
                self.score([l.id for l in self.learners.values()])

        self.openContext("", data)
        self.paintscores()


    def fillOptParams(self):
        self.error(0)
        paramsLabels = []
        #Get the Learner Name: Use only Confusion M;Atrix CM from the first and unique learner 
        learners = [(l.time, l) for l in self.learners.values()]
        learners.sort()
        learner = [lt[1] for lt in learners][0].learner
       
        # check if learner is defined in AZLearnersParamsConfig
        self.learnerType = str(learner).split()[0]
        if self.learnerType[0] == "<":
            self.learnerType = self.learnerType[1:]
        if self.learnerType.rfind('.') >=0:
            self.learnerType = self.learnerType[self.learnerType.rfind('.')+1:]

        if not hasattr(AZLearnersParamsConfig, self.learnerType):
            self.error("No configuration for specified learner " + str(self.learnerType) +" in AZLearnersParamsConfig")
            self.learnerType = None
            self.paramsLabels = []
            return
        else:
            optAPI = AZLearnersParamsConfig.API(self.learnerType)
            paramsLabels = optAPI.getParameterNames() 

        self.paramsLabels = paramsLabels
        #Select the params optimized by default
        SelectedParams = []
        for idx,param in enumerate(paramsLabels):
            if optAPI.getParameter(param,"optimize"):
                SelectedParams.append(idx)
        self.SelectedParams = SelectedParams    


    def fillClassCombo(self):
        """upon arrival of new data appropriately set the target class combo"""
        self.targetCombo.clear()
        if not self.data or not self.data.domain.classVar or not self.isclassification():
            return

        domain = self.data.domain
        self.targetCombo.addItems([str(v) for v in domain.classVar.values])
        
        if self.targetClass<len(domain.classVar.values):
            self.targetCombo.setCurrentIndex(self.targetClass)
        else:
            self.targetCombo.setCurrentIndex(0)
            self.targetClass=0

    def setLearner(self, learner, id=None):
        """add/remove a learner"""
        if learner: # a new or updated learner
            if id in self.learners: # updated learner
                time = self.learners[id].time
                self.learners[id] = Learner(learner, id)
                self.learners[id].time = time
            else: # new learner
                self.learners[id] = Learner(learner, id)
            if self.applyBtn.isEnabled():
                self.recompute(True)
            else:
                self.score([id])
            self.fillOptParams()
        else: # remove a learner and corresponding results
            if id in self.learners:
                res = self.learners[id].results
                if res and res.numberOfLearners > 1:
                    indx = [l.id for l in res.learners].index(id)
                    res.remove(indx)
                    del res.learners[indx]
                del self.learners[id]
            self.sendResults()
        self.paintscores()

    # handle output signals

    def sendResults(self):
        """commit evaluation results"""
        # for each learner, we first find a list where a result is stored
        # and remember the corresponding index

        return # Cannot send anything

        valid = [(l.results, [x.id for x in l.results.learners].index(l.id))
                 for l in self.learners.values() if l.scores and l.results]
            
        if not (self.data and len(valid)):
            self.send("Evaluation Results", None)
            return

        # find the result set for a largest number of learners
        # and remove this set from the list of result sets
        rlist = dict([(l.results,1) for l in self.learners.values() if l.scores]).keys()
        rlen = [r.numberOfLearners for r in rlist]
        results = rlist.pop(rlen.index(max(rlen)))
        
        for (i, l) in enumerate(results.learners):
            if not l.id in self.learners:
                results.remove(i)
                del results.learners[i]
        for r in rlist:
            for (i, l) in enumerate(r.learners):
                if (r, i) in valid:
                    results.add(r, i)
                    results.learners.append(r.learners[i])
                    self.learners[r.learners[i].id].results = results
        self.send("Evaluation Results", results)
        self.results = results

    # signal processing

    def newsampling(self):
        """handle change of evaluation method"""
        if not self.applyOnAnyChange:
            self.applyBtn.setDisabled(self.applyOnAnyChange)
        else:
            if self.learners:
                self.recompute()

    def newParamSelection(self):
        """handle change in set of Parameters to optimize"""
        print "Selected: ", self.SelectedParams
        print "Labels: ", self.paramsLabels
        self.recompute()

    def newscoreselection(self):
        """handle change in set of scores to be displayed"""
        usestat = [self.selectedRScores, self.selectedCScores][self.isclassification()]
        for i in range(len(self.stat)):
            if i in usestat:
                self.tab.showColumn(i+1)
                self.tab.resizeColumnToContents(i+1)
            else:
                self.tab.hideColumn(i+1)

    def recompute(self, forced=False):
        """recompute the scores for all learners,
           if not forced, will do nothing but enable the Apply button"""
        if self.applyOnAnyChange or forced:
            self.score([l.id for l in self.learners.values()])
            self.paintscores()
            self.applyBtn.setDisabled(True)
        else:
            self.applyBtn.setEnabled(True)

    def conditionalRecompute(self, option):
        """calls recompute only if specific sampling option enabled"""
        self.recompute(False)

    def applyChange(self):
        if self.applyOnAnyChange:
            self.applyBtn.setDisabled(True)
        
    def changedTarget(self):
        self.recomputeCM()

##############################################################################
# Test the widget, run from DOS prompt

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWAZTestOptLearners()
    ow.show()
    a.exec_()

    data1 = orange.ExampleTable(r'../../doc/datasets/voting')
    data2 = orange.ExampleTable(r'../../golf')
    datar = orange.ExampleTable(r'../../auto-mpg')
    data3 = orange.ExampleTable(r'../../sailing-big')
    data4 = orange.ExampleTable(r'../../sailing-test')

    l1 = orange.MajorityLearner(); l1.name = '1 - Majority'

    l2 = orange.BayesLearner()
    l2.estimatorConstructor = orange.ProbabilityEstimatorConstructor_m(m=10)
    l2.conditionalEstimatorConstructor = \
        orange.ConditionalProbabilityEstimatorConstructor_ByRows(
        estimatorConstructor = orange.ProbabilityEstimatorConstructor_m(m=10))
    l2.name = '2 - NBC (m=10)'

    l3 = orange.BayesLearner(); l3.name = '3 - NBC (default)'

    l4 = orange.MajorityLearner(); l4.name = "4 - Majority"

    import orngRegression as r
    r5 = r.LinearRegressionLearner(name="0 - lin reg")

    testcase = 4

    if testcase == 0: # 1(UPD), 3, 4
        ow.setData(data2)
        ow.setLearner(r5, 5)
        ow.setLearner(l1, 1)
        ow.setLearner(l2, 2)
        ow.setLearner(l3, 3)
        l1.name = l1.name + " UPD"
        ow.setLearner(l1, 1)
        ow.setLearner(None, 2)
        ow.setLearner(l4, 4)
#        ow.setData(data1)
#        ow.setData(datar)
#        ow.setData(data1)
    if testcase == 1: # data, but all learners removed
        ow.setLearner(l1, 1)
        ow.setLearner(l2, 2)
        ow.setLearner(l1, 1)
        ow.setLearner(None, 2)
        ow.setData(data2)
        ow.setLearner(None, 1)
    if testcase == 2: # sends data, then learner, then removes the learner
        ow.setData(data2)
        ow.setLearner(l1, 1)
        ow.setLearner(None, 1)
    if testcase == 3: # regression first
        ow.setData(datar)
        ow.setLearner(r5, 5)
    if testcase == 4: # separate train and test data
        ow.setData(data3)
        ow.setTestData(data4)
        ow.setLearner(l2, 5)
        ow.setTestData(None)

    ow.saveSettings()
