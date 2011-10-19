from AZutilities import dataUtilities
from AZutilities import miscUtilities
"""
<name>Combi-QSAR</name>
<description>Widget to automatically select the most accurate AZOrange machine learning algorithm (with optimized model hyper-parameters) and to calculated an unbiased assessemnt of the expected generalization accuracy of the resulting QSAR model. For a detailed description, please see the Poster; http://svn.seml.astrazeneca.net/trac/CC-AZ-Orange/browser/trunk/doc/AZOrangeOpenTox.pdf </description>
<icon>icons/CombiQSAR.png</icon>
<contact>Pedro Rafael Almeida</contact>
<priority>8</priority>
"""
import string

from OWWidget import *
import OWGUI
import orange
import AZOrangeConfig as AZOC
from AZutilities import competitiveWorkflow
import pprint


class OWCombiQSAR(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Combi-QSAR model'):
        OWWidget.__init__(self, parent, signalManager, name)

        # Define the input and output channels
        self.inputs = [("Classified Examples", ExampleTable, self.setData)]
        self.outputs = [("Classifier", orange.Classifier)]

        self.queueTypes = ["NoSGE","batch.q","quick.q"] 
        self.outputModes = ["Model and statistics (unbiased wrt model selection)","Statistics for all available algorithms. Please note, no model selection."]

        self.name = name
	self.dataset = None
	self.classifier = None
        self.modelFile = ""

        self.statPath = ""
        self.queueType = 0
        self.outputSel = 0
        self.defineGUI()

    def setWindowTitle(self, caption):
        """ This is an override of the function setWindowTitle in OWBaseWidget.py 
            It is used to update the name of the classifier with a suffix corresponding to the number given by DOC 
            orngDoc.py to this widget when several widgets like this are present in the Canvas.
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
            #self.Apply()

        

    def saveModel(self):
        self.warning(0)
        if self.classifier and self.modelFile and self.dataset:
            if not hasattr(self.classifier , "write"):
                self.warning("This widget does not support the save functionality for the connected lerner.")
                return
            if not self.classifier.write(self.modelFile):
                self.warning("ERROR: model was NOT saved!")
            else:
                print "Saved model to ",self.modelFile
        else:
            if not self.dataset:
                self.warning("ERROR: Missing data")
            elif not self.modelFile:
                self.warning("ERROR: Missing a save path")
            else:
                self.warning("ERROR: Unexpected error!")

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


    def destroy(self, dw = 1, dsw = 1):
	self.linksOut.clear()
	if self.dataset:
	    del self.dataset
	if self.classifier:
	    del self.classifier



    def defineGUI(self):
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Name', \
                       tooltip='Name to be used by other widgets to identify your classifier.<br>This should be a unique name!')

        
        # Queue radio buttons
        self.queueBox = OWGUI.radioButtonsInBox(self.controlArea, self, "outputSel", box="Output mode",
                                             btnLabels=self.outputModes,
                                             callback=self.changeOutputMode)

        # Queue radio buttons
        self.queueBox = OWGUI.radioButtonsInBox(self.controlArea, self, "queueType", box="Execution Mode",
                                             btnLabels=self.queueTypes,
                                             callback=None)

        # Set location of statistics file
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving the statistics results", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "statPath", labelWidth=80,  orientation = "horizontal", tooltip = "Please use full path to results file to be created.")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Choose the dir where to save.")
        button.setMaximumWidth(25)

        

        # Apply the settings and send the model to the output channel
        OWGUI.button(self.controlArea, self,"&Apply settings", callback=self.Apply)
 
        # Set location of model file
        self.boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving Model", addSpace = True, orientation=0)
        self.L1 = OWGUI.lineEdit(self.boxFile, self, "modelFile", labelWidth=80,  orientation = "horizontal", tooltip = "Once a model is created (connect this widget with a data widget), \nit can be saved by giving a file name here and clicking the save button.")
        self.L1.setMinimumWidth(200)
        self.button = OWGUI.button(self.boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        self.button.setMaximumWidth(25)

        # Save the model
        self.saveButton = OWGUI.button(self.controlArea, self,"&Save model", callback=self.saveModel)

        
        # Statistics show
        statBox = OWGUI.widgetBox(self.mainArea, "Statistics", addSpace = True, orientation=0)
        self.statInfo = OWGUI.widgetLabel(statBox, '')

        self.changeOutputMode()
	self.adjustSize()

    def changeOutputMode(self):
        self.warning(0)
        self.error(0)
        if self.outputSel != 0:
            state = False
        else:
            state = True
        self.boxFile.setEnabled(state)
        self.L1.box.setEnabled(state)
        self.button.setEnabled(state)
        self.saveButton.setEnabled(state)


    def Apply(self):
        self.warning(0)
        self.error(0)

        if not self.dataset:
            self.warning("Missing input data")
            return

        progressSteps = 2
        progress = QProgressDialog("Running Combi-QSAR.\nThis may take a while. Please wait....", str(), 0, progressSteps , None, Qt.Dialog )
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.forceShow()
        progress.setValue(1)

        if self.outputSel == 0:
            self.getModel()
        else:   
            self.getStatistics()
        progress.close()


    def getModel(self):
        self.warning(0)
        """Get the model and send it to the output"""
        #if not os.path.isdir(str(self.statPath)):
        #    statPath = None
        #else:
        #    statPath = os.path.join(str(self.statPath),"modelStat.pkl")
        statPath = None
        

        model = competitiveWorkflow.getModel(self.dataset, savePath = statPath, queueType = self.queueTypes[self.queueType], verbose = 0, getAllModels = False)
        
        self.classifier = model[model.keys()[0]] 
        self.classifier.name = str(self.name)
        self.send("Classifier", self.classifier)


    def getStatistics(self):
        self.warning(0)
        """Get the statistics and save to desired place if specifyes"""
        if not os.path.isdir(str(self.statPath)):
            statPath = None
        else:
            statPath =  os.path.join(str(self.statPath),"statistics.pkl")
        runPath = miscUtilities.createScratchDir(desc = "CombiQSAR")
        statistics = competitiveWorkflow.getStatistics(self.dataset, runPath, statPath, queueType = self.queueTypes[self.queueType], getAllModels = False)

        print pprint.pformat(statistics)
        self.statInfo.setText(pprint.pformat(statistics))
        self.send("Classifier", None)


    def setData(self, dataset): 
        self.warning(0)
        if dataset:
            self.dataset = dataset
            self.Apply()
        else:
            self.dataset = None


if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWCombiQSAR() 
    appl.setMainWidget(ow) 
    ow.show() 
    dataset = dataUtilities.DataTable('iris.tab') 
    ow.data(dataset) 
    appl.exec_loop()


