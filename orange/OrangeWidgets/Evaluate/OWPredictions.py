"""
<name>Predictions</name>
<description>Calculate predictions of Classifiers for a separate data set. Connect to 'Data Table' to see the predictions.</description>
<icon>icons/Predictions.png</icon>
<contact>Blaz Zupan (blaz.zupan(@at@)fri.uni-lj.si), Modifyed ny Pedro Almeida</contact>
<priority>300</priority>
"""
from AZutilities import dataUtilities
from OWWidget import *
import OWGUI
import statc
import time

class OWPredictions(OWWidget):
    settingsList = ["ShowAttributeMethod"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "Predictions")

        self.callbackDeposit = []
        self.inputs = [("Examples", ExampleTable, self.setData), ("Predictors", orange.Classifier, self.setPredictor, Multiple)]
        self.outputs = [("Predictions", ExampleTable)]
        self.predictors = {}

        # saveble settings
        self.ShowAttributeMethod = 0
        self.classes = []
        self.selectedClasses = []
        self.loadSettings()
        self.datalabel = "N/A"
        self.predictorlabel = "N/A"
        self.tasklabel = "N/A"
        self.outvar = None # current output variable (set by the first predictor/data set send in)
        self.classifications = []
        self.rindx = None        
        self.data = None
        self.verbose = 0
        self.nVarImportance = 0

        # GUI - Options

        # Options - classification
        ibox = OWGUI.widgetBox(self.controlArea, "Info")
        OWGUI.label(ibox, self, "Data: %(datalabel)s")
        OWGUI.label(ibox, self, "Predictors: %(predictorlabel)s")
        OWGUI.label(ibox, self, "Task: %(tasklabel)s")
        OWGUI.separator(ibox)
        OWGUI.label(ibox, self, "Predictions can be viewed with the 'Data Table'\nand saved with the 'Save' widget!")

        OWGUI.separator(self.controlArea)
        
        self.copt = OWGUI.widgetBox(self.controlArea,"Probabilities (classification)")
        self.copt.setDisabled(1)

        self.lbcls = OWGUI.listBox(self.copt, self, "selectedClasses", "classes",
                                   selectionMode=QListWidget.MultiSelection)
        self.lbcls.setFixedHeight(50)
        OWGUI.separator(self.controlArea)
        self.VarImportanceBox = OWGUI.doubleSpin(self.controlArea, self, "nVarImportance", 0,9999999,1, label="Variable importance",  orientation="horizontal", tooltip="The number of variables to report for each prediction.")

        OWGUI.checkBox(self.controlArea, self, 'verbose', 'Verbose', tooltip='Show detailed info while predicting. This will slow down the predictions process!') 
        OWGUI.separator(self.controlArea)
        self.outbox = OWGUI.widgetBox(self.controlArea,"Output")
        self.apply = OWGUI.button(self.outbox, self, "Apply", callback=self.sendpredictions)

        OWGUI.rubber(self.controlArea)
        self.adjustSize()


    ##scPA
    def onDeleteWidget(self):
        pass
    ##ecPA


    ##############################################################################
    # Contents painting
 
    def checkenable(self):
        # following should be more complicated and depends on what data are we showing
        cond = (self.outvar != None) and (self.data != None)
        #self.att.setEnabled(cond)
        self.copt.setEnabled(cond)
        e = (self.data and (self.data.domain.classVar <> None) + len(self.predictors)) >= 2
        # need at least two classes to compare predictions
    
    def clear(self):
        self.send("Predictions", None)
        self.checkenable()
        if len(self.predictors) == 0:
            self.outvar = None
            self.classes = []
            self.selectedClasses = []
            self.predictorlabel = "N/A"
    
    ##############################################################################
    # Input signals

    def setData(self, data):
        if not data:
            self.data = data
            self.datalabel = "N/A"
            self.clear()
            self.sendpredictions()
        else:
            vartypes = {1:"discrete", 2:"continuous"}
            self.data = data
            self.rindx = range(len(self.data))
            self.sendpredictions()
            self.datalabel = "%d instances" % len(data)
        self.checkenable()

    def setPredictor(self, predictor, id):
        """handles incoming classifier (prediction, as could be a regressor as well)"""

        def getoutvar(predictors):
            """return outcome variable, if consistent among predictors, else None"""
            if not len(predictors):
                return None
            ov = predictors[0].classVar
            for predictor in predictors[1:]:
                if str(ov) != str(predictor.classVar):
                    self.warning(0, "Mismatch in class variable (e.g., predictors %s and %s)" % (predictors[0].name, predictor.name))
                    return None
            return ov
        # remove the classifier with id, if empty
        if not predictor:
            if self.predictors.has_key(id):
                del self.predictors[id]
                if len(self.predictors) == 0:
                    self.clear()
                else:
                    self.predictorlabel = "%d" % len(self.predictors)
                self.sendpredictions()
            return

        # set the classifier
        self.predictors[id] = predictor
        self.predictorlabel = "%d" % len(self.predictors)

        # set the outcome variable
        ov = getoutvar(self.predictors.values())
        if len(self.predictors) and not ov:
            self.tasklabel = "N/A (type mismatch)"
            self.classes = []
            self.selectedClasses = []
            self.clear()
            self.outvar = None
            print "NO outvar"
            return
        self.warning(0) # clear all warnings

        if ov != self.outvar:
            self.outvar = ov
            # regression or classification?
            if self.outvar.varType == orange.VarTypes.Continuous:
                self.copt.hide();
                self.tasklabel = "Regression"
            else:
                self.copt.show()
                self.classes = [str(v) for v in self.outvar.values]
                self.selectedClasses = []
                self.tasklabel = "Classification"

        if self.data:
            self.sendpredictions()
        self.checkenable()

    ##############################################################################
    # Ouput signals
    ##scPA
    def sendpredictions(self):
        def getValue(c,ex,getProb = None):
            """ Get the predicted value of ex using classifier c and gets the probability of symbol of order getProb"""
            if getProb != None:
                theValue = c(ex,orange.GetProbabilities)
                if hasattr(c,"isRealProb") and not c.isRealProb():
                    self.warning(0,"The probabilities are not available in this particular case.")
                    return orange.Value('?')
            else:
                theValue = c(ex)
                #print "theValue: ",theValue
            if theValue:
                if getProb != None:
                    return orange.Value(theValue[getProb])
                else:
                    return orange.Value(theValue)
            else:
                self.warning(0,"Some example(s) were not able to be predicted. Check the domain compatibility of train and test datasets!")
                return orange.Value("?")

        self.warning(0)
        if not self.data or not self.outvar:
            self.send("Predictions", None)
            return
        messages = []
        # predictions, data set with class predictions
        classification = self.outvar.varType == orange.VarTypes.Discrete

        # create a new domain for the new data handling the predictions
        domain = orange.Domain(self.data.domain.attributes + [self.data.domain.classVar])

        # Add to the predictions the original meta attributes present in Data
        domain.addmetas(self.data.domain.getmetas())

        # Create the new Data Table containing the Data and the Predictions
        predictions = dataUtilities.DataTable(domain, self.data)

        # The number of examples to be predicted
        nEx = len(self.data)
        # the number of Learners
        nL = len(self.predictors)
        # The number of calculated iteractions
        nIter = 1
        # the number of iterations Done
        iter = 0

        self.progressBarSet(0)
        self.progressBarInit()
        if self.verbose:
            for c in self.predictors.values():
                c.verbose = int(self.verbose)
        if classification:
            if len(self.selectedClasses):
                nIter = (nEx * nL * len(self.selectedClasses)) + (nEx * nL)
                for c in self.predictors.values():
                    for i in self.selectedClasses:
                        m = orange.FloatVariable(name="%s(%s)" % (c.name, str(self.outvar.values[i])))
                        domain.addmeta(orange.newmetaid(), m)
                        for ex in predictions:
                            ex[m.name] = getValue(c,orange.Example(self.data.domain,ex),i)
                            self.progressBarSet((iter*100)/nIter)
                            iter += 1
            else:
                iter = 0
                nIter = nEx * nL    
       
            for c in self.predictors.values():
                if hasattr(c,'examplesFixedLog'):
                    c.examplesFixedLog={}
                m = orange.EnumVariable(name="%s" % c.name, values = self.outvar.values)
                domain.addmeta(orange.newmetaid(), m)
                for ex in predictions:
                    ex[m.name] = getValue(c,orange.Example(self.data.domain,ex))                    
                    self.progressBarSet((iter*100)/nIter)
                    iter += 1

        else:
            # regression
            nIter = nEx * nL
            for c in self.predictors.values():
                if hasattr(c,'examplesFixedLog'):
                    c.examplesFixedLog={}
                m = orange.FloatVariable(name="%s" % c.name)
                domain.addmeta(orange.newmetaid(), m)
                for ex in predictions:
                    ex[m.name] = getValue(c,orange.Example(self.data.domain,ex))
                    self.progressBarSet((iter*100)/nIter)
                    iter += 1

        if self.verbose:
            for c in self.predictors.values():
                c.verbose = 0

        #Compute and return individual Var Importance
        iter = 0
        nIter = nEx * nL
        if self.nVarImportance > 0:
            for c in self.predictors.values():
                if hasattr(c,'getTopImportantVars'):
                    m = orange.StringVariable(name="%s" % c.name+"(Top Important Vars)")
                    domain.addmeta(orange.newmetaid(), m)
                    for ex in predictions:
                        topVars = c.getTopImportantVars(ex,self.nVarImportance)
                        if topVars:
                            if len(topVars) == 1:
                                topVars = str(topVars[0])
                            else:
                                topVars = str(topVars)
                            #if topVars not in m.values: 
                            #    m.values.append(topVars)
                            ex[m.name] = topVars
                        else:
                            ex[m.name] = "?"
                        self.progressBarSet((iter*100)/nIter)
                        iter += 1


        self.progressBarFinished()
        for c in self.predictors.values():
                if hasattr(c,'examplesFixedLog') and ('Missing Attributes' in c.examplesFixedLog) and c.examplesFixedLog['Missing Attributes']:
                         missingAttrs = ""
                         for attr in c.examplesFixedLog['Missing Attributes']:
                                missingAttrs += "  " + attr + "\\n"
                         messages.append("QMessageBox.warning( None, \"Missing Attributes\" ,"+\
                                "\"The following attributes were missing in the examples to be predicted:\\n" + \
                                missingAttrs + "\", QMessageBox.Ok)")
                if hasattr(c,'examplesFixedLog') and ('Fixed Types of variables' in c.examplesFixedLog):
                        if 'Vars needing type fix' in c.examplesFixedLog:
                            msg = "Some variable types were fixed while predicting with "+c.name+"!\\nTypes Fixed: \\n"
                            for var in c.examplesFixedLog['Vars needing type fix']:
                                msg+="  "+var+": "+str(c.examplesFixedLog['Vars needing type fix'][var])+'\\n'
                        else:
                            msg = "Some variable types were fixed while predicting with "+c.name+"!"
                        messages.append("QMessageBox.warning( None, \"Fixed Types of variables in "+\
                                str(c.examplesFixedLog['Fixed Types of variables'])+
                                " examples\",\""+ msg+"\", QMessageBox.Ok)")
                if hasattr(c,'examplesFixedLog') and ('Fixed Order of variables' in c.examplesFixedLog):
                        messages.append("QMessageBox.warning( None, \"Fixed Order of variables in "+\
                                str(c.examplesFixedLog['Fixed Order of variables'])+
                                " examples\",\"The order of variables in test data was not the same has in the original\\n"+\
                                "training set used on "+c.name+", so they were fixed.\", QMessageBox.Ok)")
                if hasattr(c,'examplesFixedLog') and ('Fixed Number of variables' in c.examplesFixedLog):
                        messages.append("QMessageBox.warning( None, \"Fixed Number of variables in "+\
                                str(c.examplesFixedLog['Fixed Number of variables'])+
                                " examples\",\"The Number of variables in test data were not the same has in the original\\n"+\
                                "training set used on "+c.name+", so only the variables\\npresent in the training set were used .\", QMessageBox.Ok)")
        predictions.name = self.data.name
        self.send("Predictions", predictions)

        for msg in messages:
            exec(msg) in globals()
        ##ecPA


##############################################################################
# Test the widget, run from DOS prompt

if __name__=="__main__":    
    a = QApplication(sys.argv)
    ow = OWPredictions()
    a.setMainWidget(ow)
    ow.show()

    import orngTree

    dataset = dataUtilities.DataTable('../../doc/datasets/iris.tab')
#    dataset = dataUtilities.DataTable('../../doc/datasets/auto-mpg.tab')
    ind = orange.MakeRandomIndices2(p0=0.5)(dataset)
    data = dataset.select(ind, 0)
    test = dataset.select(ind, 1)
    testnoclass = dataUtilities.DataTable(orange.Domain(test.domain.attributes, False), test)        
    tree = orngTree.TreeLearner(data)
    tree.name = "tree"
    maj = orange.MajorityLearner(data)
    maj.name = "maj"
    knn = orange.kNNLearner(data, k = 10)
    knn.name = "knn"

    if 0: # data set only
        ow.setData(test)
    if 0: # two predictors, test data with class
        ow.setPredictor(maj, 1)
        ow.setPredictor(tree, 2)
        ow.setData(test)
    if 0: # two predictors, test data with no class
        ow.setPredictor(maj, 1)
        ow.setPredictor(tree, 2)
        ow.setData(testnoclass)
    if 1: # three predictors
        ow.setPredictor(tree, 1)
        ow.setPredictor(maj, 2)
        ow.setData(data)
        ow.setPredictor(knn, 3)
    if 0: # just classifier, no data
        ow.setData(data)
        ow.setPredictor(maj, 1)
        ow.setPredictor(knn, 2)
    if 0: # change data set
        ow.setPredictor(maj, 1)
        ow.setPredictor(tree, 2)
        ow.setData(testnoclass)
        data = dataUtilities.DataTable('../../doc/datasets/titanic.tab')
        tree = orngTree.TreeLearner(data)
        tree.name = "tree"
        ow.setPredictor(tree, 2)
        ow.setData(data)

    a.exec_loop()
    ow.saveSettings()
