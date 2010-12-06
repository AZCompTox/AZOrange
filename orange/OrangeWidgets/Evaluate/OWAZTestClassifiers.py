"""
<name>Test Classifiers</name>
<description>Tests the accuracy of Classifiers(trained models) on a separate (temporal/external) test data set.</description>
<icon>icons/TestClassifiers.png</icon>
<contact>Pedro Almeida (Based on OWTestLearners)</contact>
<priority>200</priority>
"""
#
# OWAZTestClassifiers.py
#
from OWWidget import *
import orngTest, orngStat, OWGUI
import time
import warnings
from AZutilities import evalUtilities

warnings.filterwarnings("ignore", "'id' is not a builtin attribute",
                        orange.AttributeWarning)

##############################################################################

class Classifier:
    def __init__(self, classifier, id):
        classifier.id = id
        self.classifier = classifier
        self.name = classifier.name
        self.id = id
        self.scores = []
        self.results = None
        self.isCompatible = False   #Flag indicating if the classifier is compatible with the input dataset
        self.time = time.time() # used to order the classifiers in the table

class Score:
    def __init__(self, name, label, f, show=True, cmBased=False, evalModule="orngStat"):
        self.name = name
        self.label = label
        self.f = f
        self.show = show
        self.cmBased = cmBased
        self.evalModule = evalModule

class OWAZTestClassifiers(OWWidget):
    settingsList = ["precision",
                    "selectedCScores", "selectedRScores"]
    contextHandlers = {"": DomainContextHandler("", ["targetClass"])}
    callbackDeposit = []


    cStatistics = [Score(*s) for s in [\
        ('Classification accuracy', 'CA', 'CA(res)', True),
        ('Sensitivity', 'Sens', 'sens(cm)', True, True),
        ('Specificity', 'Spec', 'spec(cm)', True, True),
        ('Area under ROC curve', 'AUC', 'AUC(res)', True),
        ('Information score', 'IS', 'IS(res)', False),
        ('F-measure', 'F1', 'F1(cm)', False, True),
        ('Precision', 'Prec', 'precision(cm)', False, True),
        ('Recall', 'Recall', 'recall(cm)', False, True),
        ('Brier score', 'Brier', 'BrierScore(res)', True),
        ('Matthews correlation coefficient', 'MCC', 'MCC(cm)', False, True)]]

    rStatistics = [Score(*s) for s in [\
        ("Mean squared error", "MSE", "MSE(res)", False),
        ("Root mean squared error", "RMSE", "RMSE(res)"),
        ("Mean absolute error", "MAE", "MAE(res)", False),
        ("Relative squared error", "RSE", "RSE(res)", False),
        ("Root relative squared error", "RRSE", "RRSE(res)"),
        ("Relative absolute error", "RAE", "RAE(res)", False),
        ("Q-squared", "Q2", "Q2(res,trainMeans)",True,False,"evalUtilities"),
        ("R-squared", "R2", "R2(res)")]]


    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "TestClassifiers")

        self.inputs = [("Test Data", ExampleTable, self.setTestData), ("Classifier", orange.Classifier, self.setClassifier, Multiple)]
        self.outputs = [("Evaluation Results", orngTest.ExperimentResults)]

        # Settings
        self.selectedCScores = [i for (i,s) in enumerate(self.cStatistics) if s.show]
        self.selectedRScores = [i for (i,s) in enumerate(self.rStatistics) if s.show]
        self.targetClass = 0
        self.loadSettings()

        self.stat = self.cStatistics

        self.testdata = None            # separate test data set
        self.classifiers = {}              # set of classifiers (input)
        self.results = None             # from orngTest
        self.precision = 4

        # GUI
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
        self.filePath=os.path.join(os.getcwd(),"TCResults.txt")
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
        self.g = OWGUI.widgetBox(self.mainArea, 'Evaluation Results')
        self.tab = OWGUI.table(self.g, selectionMode = QTableWidget.NoSelection)

        #self.lab = QLabel(self.g)

        self.resize(680,470)

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
        if self.results!=None and self.tab.horizontalHeader().count() and self.tab.rowCount()>0 and self.tab.columnCount()>0:
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

    # scoring and painting of score table
    def isclassification(self):
        if not self.testdata or not self.testdata.domain.classVar:
            return True
        return self.testdata.domain.classVar.varType == orange.VarTypes.Discrete


    def paintscores(self):
        """paints the table with evaluation scores"""

        self.tab.setColumnCount(len(self.stat)+1)
        self.tab.setHorizontalHeaderLabels(["Method"] + [s.label for s in self.stat])
        
        prec="%%.%df" % self.precision

        classifiers = [(l.time, l) for l in self.classifiers.values()]
        classifiers.sort()
        classifiers = [lt[1] for lt in classifiers]

        self.tab.setRowCount(len(self.classifiers))
        for (i, l) in enumerate(classifiers):
            OWGUI.tableItem(self.tab, i,0, l.name)
            
        for (i, l) in enumerate(classifiers):
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
        usestat = [self.selectedRScores, self.selectedCScores][self.isclassification()]
        for i in range(len(self.stat)):
            if i not in usestat:
                self.tab.hideColumn(i+1)


    def sendReport(self):
        self.reportData(self.testdata)
        
        self.reportSection("Results")
        classifiers = [(l.time, l) for l in self.classifiers.values()]
        classifiers.sort()
        classifiers = [lt[1] for lt in classifiers]
        usestat = [self.selectedRScores, self.selectedCScores][self.isclassification()]
        
        res = "<table><tr><th></th>"+"".join("<th><b>%s</b></th>" % hr for hr in [s.label for i, s in enumerate(self.stat) if i in usestat])+"</tr>"
        for i, l in enumerate(classifiers):
            res += "<tr><th><b>%s</b></th>" % l.name
            if l.scores:
                for j in usestat:
                    scr = l.scores[j]
                    res += "<td>" + ("%.4f" % scr if scr is not None else "") + "</td>"
            res += "</tr>"
        res += "</table>"
        self.reportRaw(res)
            
    def score(self, ids):
        """compute scores for the list of classifiers"""
        if (not self.testdata):
            for id in ids:
                self.classifiers[id].results = None
            return
        # test which classifiers can accept the given data set
        # e.g., regressions can't deal with classification data
        classifiers = []
        self.warning(0)
        for l in [self.classifiers[id] for id in ids]:
            try:
                l.isCompatible = False
                if l.classifier(self.testdata[0]).varType == self.testdata.domain.classVar.varType:
                    classifiers.append(l.classifier)
                    l.isCompatible = True
                else:
                    l.scores = []
            except Exception, ex:
                self.warning(0, "Classifier %s ends with exception: %s" % (l.name, str(ex)))
                l.scores = []

        if not classifiers:
            return


        # computation of results (res, and cm if classification)
        
        pb = OWGUI.ProgressBar(self, iterations=len(classifiers))
        res = orngTest.testOnData(classifiers, self.testdata,None,1,callback=pb.advance)
        pb.finish()
        trainMeans = []
        for c in classifiers:
            if c.classVar and hasattr(c,"basicStat") and c.classVar.name in c.basicStat and c.basicStat[c.classVar.name] != None:
                trainMeans.append(c.basicStat[c.classVar.name]["avg"])
            else:
                trainMeans.append(None)

        if self.isclassification():
            cm = orngStat.computeConfusionMatrices(res, classIndex = self.targetClass)

        res.classifiers = classifiers
        for l in [self.classifiers[id] for id in ids]:
            if l.classifier in classifiers:
                l.results = res

        self.error(range(len(self.stat)))
        scores = []
        for i, s in enumerate(self.stat):
            try:
                scores.append(eval(s.evalModule+"." + s.f))
            except Exception, ex:
                self.error(i, "An error occurred while evaluating " + s.evalModule + "." + s.f + "on %s due to %s" % \
                           (" ".join([l.name for l in classifiers]), ex.message))
                scores.append([None] * len(self.classifiers))
                
        for (i, l) in enumerate(classifiers):
            self.classifiers[l.id].scores = [s[i] if s else None for s in scores]
            
        self.sendResults()

    def recomputeCM(self):
        if not self.results:
            return
        cm = orngStat.computeConfusionMatrices(self.results, classIndex = self.targetClass)
        scores = [(indx, eval(s.evalModule+"." + s.f))
                  for (indx, s) in enumerate(self.stat) if s.cmBased]
        for (indx, score) in scores:
            for (i, l) in enumerate([l for l in self.classifiers.values() if l.scores]):
                l.scores[indx] = score[i]
        self.paintscores()
        

    # handle input signals

    def setTestData(self, data):
        """handle test data set"""
        self.closeContext()
        if not data:
            self.testdata = None
            # data was removed, remove the scores
            for l in self.classifiers.values():
                l.scores = []
                l.results = None
            self.send("Evaluation Results", None)
        else:
            self.testdata = orange.Filter_hasClassValue(data)
            if len(self.testdata) != 0:
                self.fillClassCombo()
                # new data has arrived
                if self.isclassification():
                    self.rstatLB.box.hide()
                    self.cbox.show()
                else:
                    self.cbox.hide()
                    self.rstatLB.box.show()
                self.stat = [self.rStatistics, self.cStatistics][self.isclassification()]

                if self.classifiers:
                    self.score([l.id for l in self.classifiers.values()])
                self.recompute()
            else:
                self.testdata = None
                for l in self.classifiers.values():
                    l.scores = []
                    l.results = None
                self.send("Evaluation Results", None)
                self.warning("The dataset does not have class responses for evauation.")

        self.openContext("", data)
        self.paintscores()





    def fillClassCombo(self):
        """upon arrival of new data appropriately set the target class combo"""
        self.targetCombo.clear()
       
        if not self.testdata or not self.isclassification():
            return
    
        domain = self.testdata.domain
        self.targetCombo.addItems([str(v) for v in domain.classVar.values])
        
        if self.targetClass<len(domain.classVar.values):
            self.targetCombo.setCurrentIndex(self.targetClass)
        else:
            self.targetCombo.setCurrentIndex(0)
            self.targetClass=0

    def setClassifier(self,classifier, id=None):
        """add/remove a Classifier"""
        if classifier: # a new or updated classifier
            if id in self.classifiers: # updated classifier
                time = self.classifiers[id].time
                self.classifiers[id] = Classifier(classifier, id)
                self.classifiers[id].time = time
            else: # new classifier
                self.classifiers[id] = Classifier(classifier, id)
            self.recompute()
        else: # remove a classifier and corresponding results
            if id in self.classifiers:
                res = self.classifiers[id].results
                if res and res.numberOfLearners > 1:
                    indx = [l.id for l in res.classifiers].index(id)
                    res.remove(indx)
                    del res.classifiers[indx]
                del self.classifiers[id]
            self.sendResults()
        self.paintscores()

    # handle output signals

    def sendResults(self):
        """commit evaluation results"""
        # for each classifier, we first find a list where a result is stored
        # and remember the corresponding index

        valid = [(l.results, [x.id for x in l.results.classifiers].index(l.id))
                 for l in self.classifiers.values() if l.scores and l.results]
            
        if not (self.testdata and len(valid)):
            self.send("Evaluation Results", None)
            return

        # find the result set for a largest number of classifiers
        # and remove this set from the list of result sets
        rlist = dict([(l.results,1) for l in self.classifiers.values() if l.scores]).keys()
        rlen = [r.numberOfLearners for r in rlist]
        results = rlist.pop(rlen.index(max(rlen)))
        
        for (i, l) in enumerate(results.classifiers):
            if not l.id in self.classifiers:
                results.remove(i)
                del results.classifiers[i]
        for r in rlist:
            for (i, l) in enumerate(r.classifiers):
                if (r, i) in valid:
                    results.add(r, i)
                    results.classifiers.append(r.classifiers[i])
                    self.classifiers[r.classifiers[i].id].results = results
        self.send("Evaluation Results", results)
        self.results = results

    # signal processing

    def newscoreselection(self):
        """handle change in set of scores to be displayed"""
        usestat = [self.selectedRScores, self.selectedCScores][self.isclassification()]
        for i in range(len(self.stat)):
            if i in usestat:
                self.tab.showColumn(i+1)
                self.tab.resizeColumnToContents(i+1)
            else:
                self.tab.hideColumn(i+1)

    def recompute(self):
        self.score([l.id for l in self.classifiers.values()])
        self.paintscores()
        #self.applyBtn.setDisabled(True)

        
    def changedTarget(self):
        self.recomputeCM()

##############################################################################
# Test the widget, run from DOS prompt

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWAZTestClassifiers()
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
