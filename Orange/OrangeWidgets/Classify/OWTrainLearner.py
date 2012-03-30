from AZutilities import dataUtilities
"""
<name>Train Learner</name>
<description>Train any learner obtaining a classifier</description>
<icon>icons/Train.png</icon>
<contact>Pedro Rafael Almeida</contact>
<priority>8</priority>
"""
import string

from OWWidget import *
import OWGUI
import orange
import AZOrangeConfig as AZOC

class OWTrainLearner(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Trained learner'):
        OWWidget.__init__(self, parent, signalManager, name)

        # Define the input and output channels
        self.inputs = [("Classified Examples", ExampleTable, self.setData), ("Learner", orange.Learner, self.setLearner)]
        self.outputs = [("Classifier", orange.Classifier)]

        self.name = name
	self.dataset = None
	self.classifier = None
	self.learner = None
        self.modelFile = os.path.join(os.getcwd(),"TrainedLearner.model")

        self.defineGUI()

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
            self.createClassifier()

        

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
            elif not self.classifier or not self.learner:
                self.warning("ERROR: Missing Learner")
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
	if self.learner:
	    del self.learner
	if self.classifier:
	    del self.classifier



    def defineGUI(self):
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Name', \
                       tooltip='Name to be used by other widgets to identify your classifier.<br>This should be a unique name!')

        # Apply the settings and send the learner to the output channel
        OWGUI.button(self.controlArea, self,"&Apply settings", callback=self.createClassifier)

        # Set location of model file
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving Model", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "modelFile", labelWidth=80,  orientation = "horizontal", tooltip = "Once a model is created (connect this widget with a data widget), \nit can be saved by giving a file name here and clicking the save button.")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        button.setMaximumWidth(25)

        # Save the model
        OWGUI.button(self.controlArea, self,"&Save model", callback=self.saveModel)

	self.adjustSize()



    def createClassifier(self):
        self.warning(0)
        """Create the classifier and send it to the output"""
        
        if not self.learner or not self.dataset:
            self.classifier = None
            self.send("Classifier", self.classifier)
            return

	# Impute shall be made on learner python layer 
        self.classifier = self.learner(self.dataset)
        self.classifier.name = str(self.name)
        self.send("Classifier", self.classifier)


    def setLearner(self, learner):
        self.warning(0) 
        if learner:
            self.learner = learner
        else:
            self.learner = None
        self.createClassifier()

    def setData(self, dataset): 
        self.warning(0)
        if dataset:
            self.dataset = dataset
        else:
            self.dataset = None
        self.createClassifier()


if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWTrainLearner() 
    appl.setMainWidget(ow) 
    ow.show() 
    dataset = dataUtilities.DataTable('iris.tab') 
    ow.data(dataset) 
    appl.exec_loop()


