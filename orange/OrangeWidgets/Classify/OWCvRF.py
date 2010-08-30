from AZutilities import dataUtilities
"""
<name>CV-RF</name>
<description>OpenCv Random Forests learner/classifier.</description>
<icon>icons/RandomForest.png</icon>
<contact>Jonna Stalring and Pedro Almeida</contact>
<priority>2</priority>
"""

from OWWidget import *
import orngTree, OWGUI
import orngEnsemble
from exceptions import Exception
import os

import AZOrangeConfig as AZOC
from trainingMethods import AZorngRF


class OWCvRF(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Random Forests'):
        OWWidget.__init__(self, parent, signalManager, name)
        # Define the input and output channels
        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Learner", orange.Learner),("Classifier", orange.Classifier)]

        self.name = name

        # Set default values of the random forest parameters
        self.maxDepth = AZOC.RFDEFAULTDICT["maxDepth"]
        self.minSample = AZOC.RFDEFAULTDICT["minSample"]
        self.useSurrogates = AZOC.RFDEFAULTDICT["useSurrogates"]
        self.getVarVariance = AZOC.RFDEFAULTDICT["getVarVariance"]
        self.nActVars = AZOC.RFDEFAULTDICT["nActVars"]
        self.nTrees = AZOC.RFDEFAULTDICT["nTrees"]
        self.forestAcc = AZOC.RFDEFAULTDICT["forestAcc"]
        self.termCrit = AZOC.RFDEFAULTDICT["termCrit"]
        self.priorsGUI = ""
        self.priors = AZOC.RFDEFAULTDICT["priors"]
        self.stratify = AZOC.RFDEFAULTDICT["stratify"] == "true"
        self.useBuiltInMissValHandling = int(AZOC.RFDEFAULTDICT["useBuiltInMissValHandling"])
        self.modelFile = os.path.join(os.getcwd(),"RF.model")

        self.data = None

	##scPA
        self.classifier = None
        self.learner = None
        ##ecPA

        self.defineGUI()
        self.applySettings()

    ##scPA
    def setWindowTitle(self, caption):
        """ This is an override of the function setWindowTitle in OWBaseWidget.py 
            It is used to update the name of the Learner with a suffix corresponding to the number given by DOC 
            orngDoc.py to this Learner widget when several widgets like this are present in the Canvas.
            This way, when signals are sent, and the destination widget uses the name property, there will 
            be different names identifying each one.
        """
        #print "WinTitleBefore:",self.windowTitle()
        # First call the base class function in order to the captionTitle is updated with a numbered suffix
        OWWidget.setWindowTitle(self, caption)
        #print "WinTitleAfter:",self.windowTitle()
        # Get the suffix that is something like: "(3)" and add it to the name specified in the text box
        suffix = str(self.windowTitle().split(" ")[-1]).strip()
        if len(suffix)>=3 and suffix[0]=="(" and suffix[-1]==")":
            self.name = self.name + " " + suffix
            self.applySettings()
    
    ##ecPA

    def defineGUI(self):        

        OWGUI.lineEdit(self.controlArea, self, 'name', box='Learner/Classifier Name', \
                       tooltip='Name to be used by other widgets to identify your learner/classifier.<br>This should be a unique name!')

        OWGUI.separator(self.controlArea, 0, 30)

        self.bBox = self.controlArea# QGroupBox('Model parameters for the Random Forest algorithm',self.controlArea)
        #self.bBox.setTitle('Model parameters for the Random Forest algorithm')

        OWGUI.lineEdit(self.bBox, self, 'nActVars', box='The number of active attributes', \
                       tooltip='The number of randomly sampled attributes at each node.\nThe default value 0, sets the number of samples \nto the square root of the total number of samples.')
        OWGUI.lineEdit(self.bBox, self, 'maxDepth', box='Maximum branching depth.', \
                       tooltip='The maximal number of nodes along one branch')
        OWGUI.lineEdit(self.bBox, self, 'minSample', box='Minimum number of samples', \
                       tooltip='The number of samples at which to stop branching.')

        self.trainingBox=b=QButtonGroup(self.bBox)#"Termination criteria", self.bBox)
        self.trainingBox=b=OWGUI.widgetBox(self.bBox, "Termination criteria")
        self.trainingRadio = OWGUI.radioButtonsInBox(b, self, "termCrit", btnLabels= \
                             ["Number of trees in the forest", "Out of bag error to achive before termination          "])
        OWGUI.lineEdit(self.bBox, self, 'nTrees', box='The number of trees', \
                       tooltip='This is the termination criteria used if "Number of trees in the forest" is checked')
        OWGUI.lineEdit(self.bBox, self, 'forestAcc', box='OOB error', \
                       tooltip='This is the termination criteria used if "Out of bag error" is checked')
        self.priorBox = OWGUI.lineEdit(self.bBox, self, 'priorsGUI', box='Class weights ', \
                       tooltip='Ex: POS:0.8 , NEG:0.2')
        #OWGUI.separator(self.bBox, 0, 5)
        #OWGUI.checkBox(self.bBox, self, 'stratify','Stratify',  tooltip='Stratify?')
        OWGUI.checkBox(self.bBox, self, 'useBuiltInMissValHandling','Use surrogate nodes for missing value handling',  tooltip='Use the RF built in missing value handling (as originally implemented by Breiman and Cutler) instead of the AZorange imputer. The training will be more time consuming and the model larger.')

        # Apply the settings and send the learner to the output channel
        OWGUI.button(self.controlArea, self,"&Apply settings", callback=self.applySettings)

        # Get desired location of model file
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving Model", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "modelFile", labelWidth=80,  orientation = "horizontal", \
        tooltip = "Once a model is created (connect this widget with a data widget), \nit can be saved by giving a \
file name here and clicking the save button.")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        button.setMaximumWidth(25)

        # Save the model
        OWGUI.button(self.controlArea, self,"&Save model to file", callback=self.saveModel)

        self.resize(200,400)


    ##scPA
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


    def refreshParams(self):
        if not self.learner:
            return
        for attr in ('maxDepth','minSample','nActVars','nTrees','forestAcc','termCrit'):
            setattr(self, attr, str(getattr(self.learner, attr)))
        self.termCrit = int(self.termCrit)
        #The priors
        self.priorsGUI = self.priors and str(self.priors) or ""
        #The stratify parameter have to be converted
        self.stratify = self.learner.stratify=="true" 
        self.adjustSize()

    def showEvent(self, ev):
        self.refreshParams()
    ##ecPA


    def saveModel(self):
        """Write a RF classifier instance to disk """
        if self.classifier and self.modelFile and self.data:
            if not self.classifier.write(str(self.modelFile)):
                print "ERROR: model was NOT saved!"
            else:
                print "Saved CvSVM type model to ",self.modelFile
        else:
            print "ERROR: Something is missing:"
            print "   classifier:",self.classifier
            print "   modelFile:",self.modelFile
            print "   data:",self.data




    def applySettings(self):
        """Create the learner with selected settings and send to the output channel. """
        self.error(0)
        self.warning(0)
        if self.priorsGUI:
            self.priors = str(self.priorsGUI)
        else:
            self.priors = None
        self.learner = AZorngRF.RFLearner(maxDepth = str(self.maxDepth), minSample = str(self.minSample), useSurrogates = str(self.useSurrogates), \
                                          getVarVariance = str(self.getVarVariance), nActVars = str(self.nActVars), nTrees = str(self.nTrees), \
                                          forestAcc = str(self.forestAcc), termCrit = str(self.termCrit), stratify = self.stratify and "true" or "false",\
                                          priors = self.priors, useBuiltInMissValHandling = bool(self.useBuiltInMissValHandling))
        if self.learner:
            self.learner.name = self.name
        else:
            self.learner = None
            self.error(0,"ERROR: It was not possible to create a learner. Check any previous errors.")
            self.send("Learner", self.learner)
            return
	#print time.asctime(), " -> Send learner"
        self.send("Learner", self.learner)
        self.createClassifier()
        self.refreshParams()

    def setData(self,data):
        self.error(0)
        self.warning(0)

  
        tooltip = "Ex: POS:2 , NEG:1"
        if data:
            if data.domain.classVar:
                self.data=data
                # Update the tooltip for class weights
                if  data.domain.classVar.varType == orange.VarTypes.Discrete:
                    tooltip += "<br>&nbsp;Class values available:"
                    for cls_v in data.domain.classVar.values:
                        tooltip += "<br>&nbsp;&nbsp;"+str(cls_v)
                else:
                    tooltip += "<br>&nbsp;&nbsp;Continuous class does not support weights!"
            else:
                self.data=None
                self.error(0, "The dataset does not contain a class variable")
        else:
            self.data=None

        # Change the tooltip
        self.priorBox.setToolTip(tooltip)
        self.createClassifier()

    def createClassifier(self):
        """ Output a classifier. Set to None if there is no input data.  """
        if self.data:
            if self.priorsGUI:
                self.priors = str(self.priorsGUI)
            else:
                self.priors = None
            self.learner.priors = self.priors
	    #print time.asctime(), " -> self.classifier = self.learner(self.data)"
            self.classifier = self.learner(self.data)
            if self.classifier:
                self.classifier.name = self.name
            else:
                self.classifier = None
                self.error(0,"ERROR: It was not possible to create a classifyer. Check any previous errors.")
     	    #print time.asctime(), " -> self.data = self.data"
        else:
            self.classifier = None
	#print time.asctime(), " -> self.send(...)"
        self.send("Classifier", self.classifier)


    def pbchange(self, val):
	self.progressBarSet(val*100)



##############################################################################
# Test the widget, run from DOS prompt
# > python OWDataTable.py)
# Make sure that a sample data set (adult_sample.tab) is in the directory

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWCvRF()
    a.setMainWidget(ow)

    d = dataUtilities.DataTable('adult_sample')
    ow.setData(d)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
