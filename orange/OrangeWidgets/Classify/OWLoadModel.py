"""
<name>Load model</name>
<description>Reads CvANN, CvSVM, CvRF and PLS models from disc.</description>
<icon>icons/LoadModel.png</icon>
<contact>Jonna Stalring & Pedro Almeida</contact> 
<priority>1</priority>
"""

#
# OWFile.py
# The File Widget
# A widget for opening orange data files
#

from OWWidget import *
import OWGUI
import os,string
import orange
#import orngSVM_Jakulin
from trainingMethods import AZBaseClasses

class OWLoadModel(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Loaded model'):
        OWWidget.__init__(self, parent, signalManager, name)

        self.resize(400,200)

        # Define the input and output channels
        self.inputs = []
        self.outputs = [("Classifier", orange.Classifier)]
        self.name = name
        self.modelFile = os.getcwd()
        self.modelType = None
        
        #GUI

        OWGUI.lineEdit(self.controlArea, self, 'name', box='Classifier Name', \
                       tooltip='Name to be used by other widgets to identify your classifier.<br>This should be a unique name!')


        # Define the model type. SVM, SVM_Jakulin or ANN.
        #self.kernelBox=b=OWGUI.widgetBox(self.controlArea, "Select the model type you wish to load")
        #self.kernelradio = OWGUI.radioButtonsInBox(b, self, "modelType", btnLabels=["SVM", \
        #                   "ANN                                                  ", \
        #                   "Random Forest","PLS"])

        # Get location of model file
        boxFile = OWGUI.widgetBox(self.controlArea, "Model file location", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "modelFile", labelWidth=80,  orientation = "horizontal", tooltip = "Enter full path to the model file location")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browse, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        button.setMaximumWidth(25)

        # Load the model
        OWGUI.button(self.controlArea, self,"&Load model from file", callback=self.loadModel)

        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info = OWGUI.widgetLabel(box, 'No model loaded.')


        self.adjustSize()
        

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

    ##ecPA

    def browse(self):
        #if self.modelType in [0,2,3]:
        #    self.browseDir()
        #else:
        #    self.browseFile() 
        self.browseDir()
              

    def browseFile(self):
        var = os.path.realpath(str(self.modelFile))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        else:
            startfile=os.getcwd()

        filename = QFileDialog.getOpenFileName(self, "Load Model", startfile) 
        if filename:
            self.modelFile = str(filename)


    def browseDir(self):
        var = os.path.realpath(str(self.modelFile))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        else:
            startfile=os.getcwd()
        filename = QFileDialog.getExistingDirectory(self, "Load Model", startfile,QFileDialog.ShowDirsOnly)
        if filename:
            self.modelFile = str(filename)


    def loadModel(self):
        self.modelFile = str(self.modelFile)

        self.classifier = AZBaseClasses.modelRead(self.modelFile)
        self.modelType = AZBaseClasses.modelRead(self.modelFile, retrunClassifier = False)
        if not self.classifier:
            QMessageBox.information(None, "Invalid Model", "There is no valid model in the specifyed location,\nor the model is not compatible with the accepted ones.\nPlease choose the correct model path.\nThis widget can only load " + string.join(AZBaseClasses.modelRead()[:-1],", ") + " and " + AZBaseClasses.modelRead()[-1] + " models.", QMessageBox.Ok + QMessageBox.Default)  
            self.info.setText("No model Loaded") 
            self.send("Classifier", None)
            return
 
        # Model 
        self.classifier.name = self.modelType + " - "  + self.name
        print self.modelType + " model loaded from " + self.modelFile
        self.info.setText("Loaded model: "+self.modelType)
        self.send("Classifier", self.classifier)


