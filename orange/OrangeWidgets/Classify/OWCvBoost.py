"""
<name>Cv-Boost</name>
<description>OpenCV Boost learner/classifier. Boosted Tree implementation in OpenCV. Please observe that this learner can only be used for binary classification.</description>
<icon>icons/Boost.png</icon>
<contact>Pedro Almeida</contact>
<priority>5</priority>
"""
from OWWidget import *
import OWGUI
from exceptions import Exception
from trainingMethods import AZorngCvBoost
import AZOrangeConfig as AZOC

class OWCvBoost(OWWidget):
    def __init__(self, parent=None, signalManager = None, name='kNN'):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea = 0, resizingEnabled = 0)

        self.callbackDeposit = []

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Learner", orange.Learner),("Classifier", orange.Classifier)]

        self.boostType = [("Discrete Boosting","DISCRETE"),
                       ("Real Boosting","REAL"),
                       ("Logit Boosting","LOGIT"),
                       ("Gentle Boosting","GENTLE")
                            ]
        self.splitCrit = [("Default","DEFAULT"),
                       ("Gini","GINI"),
                       ("MissClass","MISCLASS"),
                       ("Squared Error","SQERR")
                            ]

        self.currentSplitCrits = []
        # Settings
        self.name = 'CvBoost'
        self.classifier = None
        self.learner = None
        self.modelFile = os.path.join(os.getcwd(),"Boost.model")
        #CVBOOSTTYPE = { "DISCRETE":0, "REAL":1, "LOGIT":2, "GENTLE":3 }
        #CVBOOSTSPLITCRIT{ "DEFAULT":0, "GINI":1, "MISCLASS":3, "SQERR":4 }
        #CVBOOSTDEFAULTDICT = {"boost_type":"DISCRETE","weak_count":100,"split_criteria":"DEFAULT","weight_trim_rate":0.95, "max_depth":1, "use_surrogates":True, "priors":None}
        for par in ("boost_type","weak_count","split_criteria","weight_trim_rate", "max_depth", "use_surrogates","priors"):
            setattr(self, par, AZOC.CVBOOSTDEFAULTDICT[par])

        self.data = None                    # input data set
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Learner/Classifier Name', \
                 tooltip='Name to be used by other widgets to identify the learner/classifier.')

        OWGUI.separator(self.controlArea)

        pars = OWGUI.widgetBox(self.controlArea, "Parameters")
        self.bType = [x[1] for x in self.boostType].index(self.boost_type)
        self.splitC = 0
        OWGUI.comboBox(pars, self, "bType", items = [x[0] for x in self.boostType], label = "Boost Type", orientation="horizontal", valueType = str,callback= self.updateSplitCritCombo)
        self.splitCombo = OWGUI.comboBox(pars, self, "splitC", items = [x[0] for x in self.splitCrit], label = "Split Criteria", orientation="horizontal", valueType = str, tooltip = "Node splitting criteria.")
        self.updateSplitCritCombo()
        if self.split_criteria not in [x[1] for x in self.currentSplitCrits]:
            print "ERROR: Invalid CvBoost default parameters in the AZOrangeConfig.py file"
        self.splitC = [x[1] for x in self.currentSplitCrits].index(self.split_criteria)

        OWGUI.spin(pars, self, "weak_count", 0, 1000, 1, None,     "Weak Count         ", orientation="horizontal", tooltip = "The number of trees used.")
        wtrBox = OWGUI.doubleSpin(pars, self, "weight_trim_rate", 0.0, 1.0, 0.01, label="Weight Trim Rate",  orientation="horizontal", tooltip = "Removal of well classified examples. Used only for computational efficiency. The greater the number, the fewer the number of examples removed.")
        wtrBox.control.setDecimals(2)
        OWGUI.spin(pars, self, "max_depth", 1, 1000, 1, None,      "Max Depth          ", orientation="horizontal", tooltip = "Maximal branching depth")
        OWGUI.checkBox(pars, self, "use_surrogates", "Use Surrogates", tooltip = "Use surrogate nodes for missing values")
        

        OWGUI.separator(self.controlArea)

        self.setLearner()
        OWGUI.button(self.controlArea, self, "&Apply", callback = self.setLearner, disabled=0)
        
        # Get desired location of model file
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving Model", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "modelFile", labelWidth=80,  orientation = "horizontal", \
        tooltip = "Once a model is created (connect this widget with a data widget), \nit can be saved by giving a \
file name here and clicking the save button.\nPlease observe that model names should not contain an extention.")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        button.setMaximumWidth(25)

        # Save the model
        OWGUI.button(self.controlArea, self,"&Save model to file", callback=self.saveModel)


        self.adjustSize()
        #self.resize(100,250)

    def browseFile(self):
        # Possible modes:
        var = os.path.realpath(str(self.modelFile))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        else:
            startfile=os.getcwd()

        filename = QFileDialog.getSaveFileName(self, "Save Model", startfile)
        if filename:
            self.modelFile = str(filename)


    def saveModel(self):
        """Write a Boost classifier instance to disk """
        self.warning(0)
        if self.classifier and self.modelFile and self.data:
            if not self.classifier.write(str(self.modelFile)):
                self.warning("Cannot save model. Please check the output window.")
                print "ERROR: model was NOT saved!"
            else:
                print "Saved CvBoost type model to ",self.modelFile
        else:
            self.warning("Cannot save model. Please check the output window.")
            print "ERROR: Something is missing:"
            print "   classifier:",self.classifier
            print "   modelFile:",self.modelFile
            print "   data:",self.data


    def setWindowTitle(self, caption):
        """ This is an override of the function setWindowTitle in OWBaseWidget.py 
            It is used to update the name of the Learner with a suffix corresponding to the number given by DOC 
            orngDoc.py to this Learner widget when several widgets like this are present in the Canvas.
            This way, when signals are sent, and the destination widget uses the name property, there will 
            be different names identifying each one.
        """
        # First call the base class function in order to the captionTitle is updated with a numbered suffix
        OWWidget.setWindowTitle(self, caption)
        # Get the suffix that is something like: "(3)" and add it to the name specified in the text box
        suffix = str(self.windowTitle().split(" ")[-1]).strip()
        if len(suffix)>=3 and suffix[0]=="(" and suffix[-1]==")":
            self.name = self.name + " " + suffix
            self.setLearner()


    def showEvent(self, ev):
        self.refreshParams()

    def sendReport(self):
        self.reportSettings("Learning parameters",
                            [("Boost Type", self.boost_type),
                             ("Split Criteria", self.split_criteria),
                             ("Used Surrogates", OWGUI.YesNo[self.use_surrogates]),
                             ("Weak Count", self.weak_count),
                             ("Weight Trim Rate", self.weight_trim_rate),
                             ("Max Depth", self.max_depth),
                             ])
        self.reportData(self.data)
        
            
    def setData(self,data):
        self.data = self.isDataWithClass(data, orange.VarTypes.Discrete) and data or None
        self.setLearner()

    def refreshParams(self):
        if self.learner:
            for par in ("boost_type","weak_count","split_criteria","weight_trim_rate", "max_depth", "use_surrogates","priors"):
                setattr(self, par, getattr(self.learner,par))

        #Convert some parameters to be used in GUI
        self.bType = [x[1] for x in self.boostType].index(self.boost_type)
        self.updateSplitCritCombo()
        if self.split_criteria not in [x[1] for x in self.currentSplitCrits]:
            self.error("Invalid configuration from Learner?!?!")
            raise(Exception("Invalid configuration from Learner?!?!"))
        self.splitC= [x[1] for x in self.currentSplitCrits].index(self.split_criteria)
        print "Parameters received from the Learner:"
        self.printLearnerPars()

    def printLearnerPars(self):
        for par in ("boost_type","weak_count","split_criteria","weight_trim_rate", "max_depth", "use_surrogates","priors"):
            print par,":",getattr(self, par)
        print "="*40

    def updateSplitCritCombo(self):
        #Refresh from boost type Combo
        self.boost_type = self.boostType[self.bType][1]
        #Fill the split combo acordingly
        self.splitCombo.clear()
        if self.boost_type == "DISCRETE":
            self.currentSplitCrits = [x for x in self.splitCrit if x[1] in ["MISCLASS", "GINI"]]
            self.splitCombo.insertItems(0,[x[0] for x in self.currentSplitCrits])
            self.splitC = [x[1] for x in self.currentSplitCrits].index("MISCLASS")
        elif self.boost_type == "REAL":
            self.currentSplitCrits = [x for x in self.splitCrit if x[1] in  ["GINI", "MISCLASS"]]
            self.splitCombo.insertItems(0,[x[0] for x in self.currentSplitCrits])
            self.splitC = [x[1] for x in self.currentSplitCrits].index("GINI")
        elif self.boost_type in ["LOGIT","GENTLE"]:
            self.currentSplitCrits = [x for x in self.splitCrit if x[1] in  ["SQERR"]]
            self.splitCombo.insertItems(0,[x[0] for x in self.currentSplitCrits])
            self.splitC = 0
        else:
            self.error("Invalid boostType?!?!")
            raise(Exception("Invalid boostType?!?!"))


    def setLearner(self):
        self.error(0)
        self.warning(0)
        self.boost_type = self.boostType[self.bType][1]
        self.split_criteria = self.currentSplitCrits[self.splitC][1] 
        print "Parameters sent to Learner:"
        self.printLearnerPars()
        self.learner = AZorngCvBoost.CvBoostLearner(boost_type = self.boost_type, weak_count = self.weak_count, split_criteria = self.split_criteria, weight_trim_rate = self.weight_trim_rate, max_depth = self.max_depth, use_surrogates = self.use_surrogates)
        self.learner.name = self.name
        self.send("Learner", self.learner)
        self.learn()


    def learn(self):
        if self.data and self.learner:
            try:
                self.classifier = self.learner(self.data)
                self.classifier.name = self.name
            except Exception, (errValue):
                self.classifier = None
                self.error(str(errValue))
            self.send("Classifier", self.classifier)
        else:
            self.send("Classifier", None)

##############################################################################
# Test the widget, run from DOS prompt
# > python OWDataTable.py)
# Make sure that a sample data set (adult_sample.tab) is in the directory

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWCvBoost()

##    dataset = orange.ExampleTable('adult_sample')
##    ow.setData(dataset)

    ow.show()
    a.exec_()
    ow.saveSettings()
