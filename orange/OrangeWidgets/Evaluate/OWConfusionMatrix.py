"""
<name>Confusion Matrix</name>
<description>Shows confusion matrices.</description>
<contact>Janez Demsar</contact>
<icon>icons/ConfusionMatrix.png</icon>
<priority>1001</priority>
"""
from OWWidget import *
import OWGUI
import orngStat, orngTest
import statc, math
from operator import add

class OWConfusionMatrix(OWWidget):
    settings = ["shownQuantity", "autoApply", "appendPredictions", "appendProbabilities"]

    quantities = ["Number of examples", "Observed and expected examples", "Proportions of predicted", "Proportions of true"]
    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "Confusion Matrix", 1)

        # inputs
        self.inputs=[("Evaluation Results", orngTest.ExperimentResults, self.setTestResults, Default)]
        self.outputs=[("Selected Examples", ExampleTable, 8)]

        self.selectedLearner = []
        self.learnerNames = []
        self.selectionDirty = 0
        self.autoApply = True
        self.appendPredictions = True
        self.appendProbabilities = False
        self.shownQuantity = 0

        self.learnerList = OWGUI.listBox(self.controlArea, self, "selectedLearner", "learnerNames", box = "Learners", callback = self.learnerChanged)
        self.learnerList.setMinimumHeight(100)
        OWGUI.separator(self.controlArea)

        OWGUI.comboBox(self.controlArea, self, "shownQuantity", items = self.quantities, box = "Show", callback=self.reprint)

        box = OWGUI.widgetBox(self.controlArea, "Selection", addSpace=True)
        OWGUI.button(box, self, "Correct", callback=self.selectCorrect)
        OWGUI.button(box, self, "Misclassified", callback=self.selectWrong)
        OWGUI.button(box, self, "None", callback=self.selectNone)

        box = OWGUI.widgetBox(self.controlArea, "Output")
        OWGUI.checkBox(box, self, "appendPredictions", "Append class predictions", callback = self.sendIf)
        OWGUI.checkBox(box, self, "appendProbabilities", "Append predicted class probabilities", callback = self.sendIf)
        applyButton = OWGUI.button(box, self, "Commit", callback = self.sendData)
        autoApplyCB = OWGUI.checkBox(box, self, "autoApply", "Commit automatically")
        OWGUI.setStopper(self, applyButton, autoApplyCB, "selectionDirty", self.sendData)

        import sip
        sip.delete(self.mainArea.layout())
        self.layout = QGridLayout(self.mainArea)

        self.layout.addWidget(OWGUI.widgetLabel(self.mainArea, "Prediction"), 0, 1, Qt.AlignCenter)
        self.layout.addWidget(OWGUI.widgetLabel(self.mainArea, "Correct Class  "), 2, 0, Qt.AlignCenter)
        self.table = OWGUI.table(self.mainArea, rows = 0, columns = 0, selectionMode = QTableWidget.MultiSelection, addToLayout = 0)
        self.layout.addWidget(self.table, 2, 1)
        self.layout.setColumnStretch(1, 100)
        self.layout.setRowStretch(2, 100)
        self.connect(self.table, SIGNAL("itemSelectionChanged()"), self.sendIf)

        ##scPA
        # Get location of model file
        self.filePath=os.path.join(os.getcwd(),"ConfusionMat.txt")
        boxFile = OWGUI.widgetBox(self.controlArea, "File for saving results", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "filePath", labelWidth=80,  orientation = "horizontal", tooltip = "Enter full file path to save the displayed results to a tab separated file.")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Browse for a location...")
        button.setMaximumWidth(25)
        # Save the model
        OWGUI.button(self.controlArea, self,"&Save Results", callback=self.saveRes)
        ##ecPA
        self.resize(700,450)


    def browseFile(self):
        # Possible modes:
        var = os.path.realpath(str(self.filePath))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        else:
            startfile=os.getcwd()

        filename = QFileDialog.getSaveFileName(self, "Save Confusion Matrix", startfile)
        if filename:
            self.filePath = str(filename)


    ##scPA
    def saveRes(self):
            if self.table!=None and self.table.rowCount()>0 and self.table.columnCount()>0:
                if os.path.isdir(os.path.split(self.filePath)[0]) and self.filePath!="" and self.filePath!=None:
                        try:
                            # Open the file
                            f = open(self.filePath, 'w')
                            nLearners =  len(self.matrix) # the number of learners connected
                            for idx in range(nLearners):
                                PMatrix = self.getPrintableConfMat(idx)
                                # Write the Learner name of the file
                                f.write(str(self.learnerNames[idx])+"\n")
                                #Save the selected learner's table to file
                                for ri, r in enumerate(PMatrix):
                                    line = ""
                                    for ci, c in enumerate(r):
                                        line += c + "\t"
                                    f.write(line[0:-1]+"\n")
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
    ##ecPA


    def setTestResults(self, res):
        ##scPA
        self.warning(0)
        self.error(0)
        ##ecPA

        self.res = res
        if not res:
            self.matrix = None
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            return

        self.matrix = orngStat.confusionMatrices(res, -2)

        dim = len(res.classValues)

        self.table.setRowCount(dim+1)
        self.table.setColumnCount(dim+1)

        self.table.setHorizontalHeaderLabels(res.classValues+[""])
        self.table.setVerticalHeaderLabels(res.classValues+[""])

        for ri in range(dim+1):
            for ci in range(dim+1):
                it = QTableWidgetItem()
                it.setFlags(Qt.ItemIsEnabled | (ri<dim and ci<dim and Qt.ItemIsSelectable or Qt.NoItemFlags))
                it.setTextAlignment(Qt.AlignRight)
                self.table.setItem(ri, ci, it)

        boldf = self.table.item(0, dim).font()
        boldf.setBold(True)
        for ri in range(dim+1):
            self.table.item(ri, dim).setFont(boldf)
            self.table.item(dim, ri).setFont(boldf)
            
        self.learnerNames = res.classifierNames[:]
        if not self.selectedLearner and self.res.numberOfLearners:
            self.selectedLearner = [0]
        self.learnerChanged()
        self.table.clearSelection()


    def learnerChanged(self):
        if not (self.res and self.res.numberOfLearners):
            return
        
        if self.selectedLearner and self.selectedLearner[0] > self.res.numberOfLearners:
            self.selectedLearner = [0]
        if not self.selectedLearner:
            return
        
        cm = self.matrix[self.selectedLearner[0]]
        
        self.isInteger = " %i "
        for r in reduce(add, cm):
            if int(r) != r:
                self.isInteger = " %5.3f "
                break

        self.reprint()
        self.sendIf()


    def getPrintableConfMat(self,idx_Learner):
        """Compute the full confusion Matrix including the sums for the specifyed learner
           and returns a Matrix with the formated strings for preatty printing
           [[ a       ,b       ,SumRow1]
            [ c       ,d       ,SumRow2]
            [ SumCol1 ,SumCol2 ,Total  ]]        
        """
        cm = self.matrix[idx_Learner]
        dim = len(cm)
        
        PMatrix = [["N/A"]*(dim+1) for x in range(dim+1)]

        rowSums = [sum(r) for r in cm]
        colSums = [sum([r[i] for r in cm]) for i in range(dim)]
        total = sum(rowSums)
        if self.shownQuantity == 1:
            if total > 1e-5:
                rowPriors = [r/total for r in rowSums]
                colPriors = [r/total for r in colSums]
            else:
                rowPriors = [0 for r in rowSums]
                colPriors = [0 for r in colSums]
        for ri, r in enumerate(cm):
            for ci, c in enumerate(r):
                if self.shownQuantity == 0:
                    PMatrix[ri][ci] = self.isInteger % c
                elif self.shownQuantity == 1:
                    PMatrix[ri][ci] = (self.isInteger + "/ %5.3f ") % (c, total*rowPriors[ri]*colPriors[ci])
                elif self.shownQuantity == 2:
                    PMatrix[ri][ci] = colSums[ci] > 1e-5 and (" %2.1f %%  " % (100 * c / colSums[ci])) or " "+"N/A"+" "
                elif self.shownQuantity == 3:
                    PMatrix[ri][ci] = rowSums[ri] > 1e-5 and (" %2.1f %%  " % (100 * c / rowSums[ri])) or " "+"N/A"+" "

        for ci in range(dim):
            PMatrix[dim][ci] = self.isInteger % colSums[ci]
            PMatrix[ci][dim] = self.isInteger % rowSums[ci]
        PMatrix[dim][dim] = self.isInteger % total

        return PMatrix

    def reprint(self):
        if not self.matrix or not self.selectedLearner: 
            return
        
        PMatrix = self.getPrintableConfMat(self.selectedLearner[0])
        for ri, r in enumerate(PMatrix):
            for ci, c in enumerate(r):
                #if not elf.table.item:
                #    continue
                self.table.item(ri, ci).setText(c)

        self.table.resizeColumnsToContents()

    def sendReport(self):
        self.reportSettings("Contents",
                            [("Learner", self.learnerNames[self.selectedLearner[0]]),
                             ("Data", self.quantities[self.shownQuantity])])
        
        self.reportSection("Matrix")
        classVals = self.res.classValues
        nClassVals = len(classVals)
        res = "<table>\n<tr><td></td>" + "".join('<td align="center"><b>&nbsp;&nbsp;%s&nbsp;&nbsp;</b></td>' % cv for cv in classVals) + "</tr>\n"
        for i, cv in enumerate(classVals):
            res += '<tr><th align="right"><b>%s</b></th>' % cv + \
                   "".join('<td align="center">%s</td>' % self.table.item(i, j).text() for j in range(nClassVals)) + \
                   '<th align="right"><b>%s</b></th>' % self.table.item(i, nClassVals).text() + \
                   "</tr>\n"
        res += '<tr><th></th>' + \
               "".join('<td align="center"><b>%s</b></td>' % self.table.item(nClassVals, j).text() for j in range(nClassVals+1)) + \
               "</tr>\n"
        res += "</table>\n<p><b>Note:</b> columns represent predictions, row represent true classes</p>"
        self.reportRaw(res)
            

    def selectCorrect(self):
        if not self.res:
            return

        sa = self.autoApply
        self.autoApply = False
        self.table.clearSelection()
        for i in range(len(self.matrix[0])):
            self.table.setRangeSelected(QTableWidgetSelectionRange(i, i, i, i), True)
        self.autoApply = sa
        self.sendIf()


    def selectWrong(self):
        if not self.res:
            return

        sa = self.autoApply
        self.autoApply = False
        self.table.clearSelection()
        dim = len(self.matrix[0])
        self.table.setRangeSelected(QTableWidgetSelectionRange(0, 0, dim-1, dim-1), True)
        for i in range(len(self.matrix[0])):
            self.table.setRangeSelected(QTableWidgetSelectionRange(i, i, i, i), False)
        self.autoApply = sa
        self.sendIf()


    def selectNone(self):
        self.table.clearSelection()


    def sendIf(self):
        if self.autoApply:
            self.sendData()
        else:
            self.selectionDirty = True


    def sendData(self):
        self.selectionDirty = False

        selected = [(x.row(), x.column()) for x in self.table.selectedIndexes()]
        res = self.res
        if not res or not selected or not self.selectedLearner:
            self.send("Selected Examples", None)
            return

        learnerI = self.selectedLearner[0]
        selectionIndices = [i for i, rese in enumerate(res.results) if (rese.actualClass, rese.classes[learnerI]) in selected]
        data = res.examples.getitemsref(selectionIndices)
        
        if self.appendPredictions or self.appendProbabilities:
            domain = orange.Domain(data.domain.attributes, data.domain.classVar)
            domain.addmetas(data.domain.getmetas())
            data = orange.ExampleTable(domain, data)
        
            if self.appendPredictions:
                cname = self.learnerNames[learnerI]
                predVar = type(domain.classVar)("%s(%s)" % (domain.classVar.name, cname.encode("utf-8") if isinstance(cname, unicode) else cname))
                if hasattr(domain.classVar, "values"):
                    predVar.values = domain.classVar.values
                predictionsId = orange.newmetaid()
                domain.addmeta(predictionsId, predVar)
                for i, ex in zip(selectionIndices, data):
                    ex[predictionsId] = res.results[i].classes[learnerI]
                    
            if self.appendProbabilities:
                probVars = [orange.FloatVariable("p(%s)" % v) for v in domain.classVar.values]
                probIds = [orange.newmetaid() for pv in probVars]
                domain.addmetas(dict(zip(probIds, probVars)))
                for i, ex in zip(selectionIndices, data):
                    for id, p in zip(probIds, res.results[i].probabilities[learnerI]):
                        ex[id] = p
    
        self.send("Selected Examples", data)


if __name__ == "__main__":
    a = QApplication(sys.argv)
    owdm = OWConfusionMatrix()
    owdm.show()
    a.exec_()
