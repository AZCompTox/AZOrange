"""
<name>PLS</name>
<description>Partial Least Squares learner/classifier</description>
<icon>icons/PLS.png</icon>
<contact>Pedro Rafael Almeida</contact>
<priority>5</priority>
"""
import string
from AZutilities import dataUtilities
from OWWidget import *
import OWGUI
import orange
from trainingMethods import AZorngPLS
import AZOrangeConfig as AZOC
import os

class OWPLS(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='PLS'):
        OWWidget.__init__(self, parent, signalManager, name)

        # Define the input and output channels
        self.inputs = [("Examples", ExampleTable, self.createClassifier)]
        self.outputs = [("Learner", orange.Learner),("Classifier", orange.Classifier)]

        self.name = 'Partial Least Squares'
        self.dataset = None
        self.classifier = None
        self.learner = None
        self.modelPath = os.path.join(os.getcwd(),"PLS.model")

        # Define the parameters readable through the GUI and set their default values
        self.optMethodDict = AZOC.OPTMETHODDICT
        self.k = AZOC.PLSDEFAULTDICT["k"] 
        self.precision = AZOC.PLSDEFAULTDICT["precision"]
        self.methodGUI = AZOC.PLSDEFAULTDICT["method"]
        self.transformOptMethod()
       
        self.defineGUI()
        self.applySettings()

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


    def onDeleteWidget(self):
        self.linksOut.clear()
        if self.dataset:
            del self.dataset
        if self.learner:
            del self.learner
        if self.classifier:
            del self.classifier
    def defineGUI(self):

        # Set the precision throught the GUI 
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Learner/Classifier Name', \
                       tooltip='Name to be used by other widgets to identify your learner/classifier.<br>This should be a unique name!')

        OWGUI.lineEdit(self.controlArea, self, 'precision', box='Precision', \
        tooltip='Desired precision.')

        # Set the number of components of PLS throught the GUI 
        OWGUI.lineEdit(self.controlArea, self, 'k', box='The number of components (factors)', \
        tooltip='The number of components to include.')

        # Set the training PLS methodGUI throught the GUI 
        self.trainingBox=b=OWGUI.widgetBox(self.controlArea, "PLS method")
        self.trainingRadio = OWGUI.radioButtonsInBox(b, self, "methodGUI", btnLabels= \
        ["KERNEL - The NIPALS-kernel method.", \
        "SIMPLS - The SIMPLS method.", \
        "PLS1 - The original PLS method."])

        # Apply the settings and send the learner to the output channel
        OWGUI.button(self.controlArea, self,"&Apply settings", callback=self.applySettings)

        #Refresh Learner Parameters
        #OWGUI.button(self.controlArea, self,"&Refresh parameters", callback=self.refreshParams)
        
        # Set location of model file
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving Model", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "modelPath", labelWidth=80,  orientation = "horizontal", tooltip = "Once a model is created (connect this widget with a data widget), \nit can be saved by giving a \
file name here and clicking the save button.")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        button.setMaximumWidth(25)

        # Save the model
        OWGUI.button(self.controlArea, self,"&Save PLS model", callback=self.saveModel)

        self.adjustSize()

    def browseFile(self):
        # Possible modes:
        var = os.path.realpath(str(self.modelPath))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        else:
            startfile=os.getcwd()

        filename = QFileDialog.getSaveFileName(self, "Save Model", startfile)
        if filename:
            self.modelPath = str(filename)


    def showEvent(self, ev):
        self.refreshParams()
 
    def refreshParams(self):
        if not self.learner:
            return
        self.precision = self.learner.precision 
        self.methodGUI = self.optMethodDict[str(self.learner.method)]
        self.k = self.learner.k

        self.adjustSize()

    def saveModel(self):
        if self.classifier!=None:
                if os.path.isdir(os.path.split(self.modelPath)[0]) and self.dataset:
                        # Model of PLS type
                        try:
                            #self.classifier.write(self.modelPath, self.dataset)
                            self.classifier.write(str(self.modelPath)) 
                            print "PLS model saved to "+str(self.modelPath)
                        except:
                            print "ERROR: Failed to Save PLS model."
                else:
                        print "Failed to Save the model. Invalid  path.";
        else:
                print "No trained model for saving!"

    def transformOptMethod(self):
        """ Transform the methodGUI parameter. OWGUI.radioButtonsInBox requires integers, 
            while AZorngPLS takes strings and in the default settings the method is defined by a string.""" 

        # Default value is the string
        if self.methodGUI in self.optMethodDict.keys():
            for key, value in self.optMethodDict.iteritems():
                if self.methodGUI == key: 
                        self.methodGUI = value
                        self.method = key

        # Value read from GUI is an int
        if self.methodGUI in self.optMethodDict.values():
            for key, value in self.optMethodDict.iteritems():
                if self.methodGUI == value: self.method = key


    def transformSettings(self):
        """Modifiy some of the parameters read by the GUI to values required by AZorngPLS. """

        # Transform the method parameter. OWGUI.radioButtonsInBox requires integers, while AZorngPLS takes strings. 
        self.transformOptMethod()
     
        # Assure that k (components) is an integer
        if type(self.k) is str: self.k = atoi(self.k)

        # Assure that precision is a float
        if type(self.precision) is str: self.precision = atof(self.precision)


    def applySettings(self):
        """Create the learner with selected settings and send to the output channel. """
        # Transform settings to those appropriate for AZorngPLS
        self.transformSettings()

        # Output a learner regardless of whether input data is provided
        self.learner = AZorngPLS.PLSLearner(k = str(self.k), method = str(self.method), \
                        precision = str(self.precision))
        self.learner.name = str(self.name)
        self.send("Learner", self.learner)
        if self.dataset:
            # Impute is made on AZorngPLS.py 
            self.classifier = self.learner(self.dataset)
            self.classifier.name = str(self.name)
        else:
            self.classifier = None
        self.send("Classifier", self.classifier)



    def createClassifier(self, dataset):
        """ Output a classifier. Set to None if there is no input data.  """
        if dataset:
            # Impute is made on AZorngPLS.py
            self.classifier = self.learner(dataset)
            self.classifier.name = str(self.name)
            #print time.asctime(), "- self.dataset = dataset"
            self.dataset = dataset
        else:
            self.classifier = None
            self.dataset = None
        #print time.asctime(), "- self.send(...)"
        self.send("Classifier", self.classifier) 

if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWPLS() 
    appl.setMainWidget(ow) 
    ow.show() 
    dataset = dataUtilities.DataTable('iris.tab') 
    ow.data(dataset) 
    appl.exec_loop()


