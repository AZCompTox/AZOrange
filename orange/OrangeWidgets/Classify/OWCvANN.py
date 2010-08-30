"""
<name>Cv-ANN</name>
<description>ANN- OpenCV Artificial Neural Network learner/classifier</description>
<icon>icons/ANN.png</icon>
<contact>Jonna Stalring & Pedro Almeida</contact>
<priority>4</priority>
"""
import string,os

from OWWidget import *
import OWGUI
import orange
from trainingMethods import AZorngCvANN
import AZOrangeConfig as AZOC
import orngImpute
from AZutilities import dataUtilities
from AZutilities import miscUtilities

class OWCvANN(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Artificial Neural Network'):
        OWWidget.__init__(self, parent, signalManager, name)

        # Define the input and output channels
        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Learner", orange.Learner),("Classifier", orange.Classifier)]

        self.name = name
        self.data = None

        # Define the parameters readable through the GUI and set their default values
        #optAlg
        self.optAlgDict = AZOC.CVANNOPTALGDICT
        self.optAlgGUI = None
        self.optAlg = AZOC.CVANNDEFAULTDICT["optAlg"]

        self.stopCritDict = AZOC.CVANNSTOPCRITDICT
        self.stopCritGUI = None
        self.stopCrit = AZOC.CVANNDEFAULTDICT["stopCrit"]
        self.maxIter = AZOC.CVANNDEFAULTDICT["maxIter"]
        self.eps = AZOC.CVANNDEFAULTDICT["eps"]
        self.stopValueGUI = None


        self.scaleData = AZOC.CVANNDEFAULTDICT["scaleData"]
        self.scaleDataGUI = None
        self.scaleClass = AZOC.CVANNDEFAULTDICT["scaleClass"]
        self.scaleClassGUI = None
        
        self.nHiddenGUI = None
        self.nHidden = AZOC.CVANNDEFAULTDICT["nHidden"]
        self.modelFile = os.path.join(os.getcwd(),"ANN.model")

        self.priors = AZOC.CVANNDEFAULTDICT["priors"] or ""

        self.setGUIVars() 
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

        # Set the number of hidden neurons throught the GUI
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Learner/Classifier Name', \
                       tooltip='Name to be used by other widgets to identify your learner/classifier.<br>This should be a unique name!')

        OWGUI.lineEdit(self.controlArea, self, 'nHiddenGUI', box='The number of hidden neurons', \
        tooltip='The number of neurons in each hidden layer should be given as a comma separated string. \n\
For example, 2,3,4 will create 3 hidden layers (5 layers in total, including in and out),\n\
with 2, 3 and 4 neurons respectively. If no value is given, the default, which is one layer of 5 neurons, \n\
will be used. ')

        # Set the training function throught the GUI 
        # CVANNOPTALGDICT = {"FANN_TRAIN_BACKPROP":0, "FANN_TRAIN_RPROP":1}
        self.trainingBox=b=OWGUI.widgetBox(self.controlArea,"Training function")
        self.trainingRadio = OWGUI.radioButtonsInBox(b, self, "optAlgGUI", btnLabels= \
        ["TRAIN_BACKPROP - Backpropagation", \
        "TRAIN_RPROP - iRPROP algorithm (Recommended)"],callback = self.changedOptAlg)

        #CVANNSTOPCRITDICT = {"ITER":1, "EPS":2}
        sc = OWGUI.widgetBox(self.controlArea,"Stop criteria")
        self.trainingRadio = OWGUI.radioButtonsInBox(sc, self, "stopCritGUI", btnLabels= \
        ["Max number of Iterations", \
         "Epsilon"], callback = self.OnChangedStopCrit) 

        OWGUI.lineEdit(OWGUI.indentedBox(sc), self, 'stopValueGUI', 'Value:', orientation="horizontal", tooltip="Value to use as stop Criteria")

        self.priorBox = OWGUI.lineEdit(b, self, 'priorsGUI', box='Class weights', \
                       tooltip='Ex: POS:2 , NEG:1')
        b=OWGUI.widgetBox(self.controlArea, "Scaling", addSpace = True)
        self.scaleDataBox=OWGUI.checkBox(b, self, "scaleDataGUI", "Scale Data", tooltip="Scales the training data, and the input examples",callback = self.setLearnerVars)
        self.scaleClassBox=OWGUI.checkBox(b, self, "scaleClassGUI", "Scale Class variable", tooltip="Scales the class variable of training data, and the output predictions",callback = self.setLearnerVars)

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
        self.setLearnerVars()
        self.applySettings()
        self.adjustSize()
##scPA
      
    def changedOptAlg(self): 
        if self.optAlgGUI == 1:
            self.priorBox.setDisabled(0)
            self.priorsGUI = ""
        else:
            self.priorBox.setDisabled(1)
            self.priorsGUI = ""        

    def setGUIVars(self):
        #CVANNOPTALGDICT = {"FANN_TRAIN_BACKPROP":0, "FANN_TRAIN_RPROP":1}
        self.optAlgGUI =  int(self.optAlg)
        self.priorsGUI = self.priors and str(self.priors) or ""
        
        #CVANNSTOPCRITDICT = {"ITER":1, "EPS":2}
        self.stopCritGUI = int(self.stopCrit) - 1     
        if self.stopCrit == 1: #ITER
            self.stopValueGUI = str(self.maxIter)
        else: #EPS
            self.stopValueGUI = str(self.eps)

        #Scale
        self.scaleDataGUI = int(self.scaleData)
        self.scaleClassGUI = int(self.scaleClass)
#        if not self.scaleData:
#            self.scaleClassBox.setDisabled(1)
#        else:
#            self.scaleClassBox.setDisabled(0)

        if type(self.nHidden) != list:
            self.nHidden = [self.nHidden]
        self.nHiddenGUI = str(self.nHidden)[1:-1]


    def setLearnerVars(self):
        if self.priorsGUI:
            self.priors = str(self.priorsGUI)
        else:   
            self.priors = None
        self.optAlg = int(self.optAlgGUI)
        self.stopCrit = int(self.stopCritGUI) + 1
        if miscUtilities.isNumber(self.stopValueGUI):
            if self.stopCrit == 1: #ITER
                self.maxIter = int(self.stopValueGUI)
            else: #EPS
                self.eps = float(self.stopValueGUI)
        else:
            self.warning(0, "Value for stop criteria must be a number! Using default value.")
            self.maxIter = AZOC.CVANNDEFAULTDICT["maxIter"]
            self.eps = AZOC.CVANNDEFAULTDICT["eps"]
        self.scaleData = bool(self.scaleDataGUI)
        self.scaleClass = bool(self.scaleClassGUI)
        # Transform the nHidden to a list
        try:
            self.nHidden = string.split(str(self.nHiddenGUI),",")
            self.nHidden = [int(elem.strip()) for elem in self.nHidden]
        except:
            self.warning(0, "Bad values for hidden neuronsi! Using default value.")
            self.nHidden = AZOC.CVANNDEFAULTDICT["nHidden"]
        self.setGUIVars()           
 
    def OnChangedStopCrit(self):
        if self.stopCritGUI == 0: #ITER
            self.stopValueGUI = AZOC.CVANNDEFAULTDICT["maxIter"]
        else: #EPS
            self.stopValueGUI = AZOC.CVANNDEFAULTDICT["eps"]
        self.setLearnerVars()


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
        for attr in ("scaleClass","scaleData","nHidden", "stopCrit", "maxIter", "optAlg", "eps", "priors"):
            setattr(self, attr, getattr(self.learner, attr))
        self.setGUIVars()
        self.adjustSize()

    def showEvent(self, ev):
        self.refreshParams()


    def saveModel(self):
        """Write an ann classifier instance to disk """
        if self.classifier and self.modelFile and self.data:
            if not self.classifier.write(self.modelFile):
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
        # Transform settings to those appropriate for AZorngCvANN
        self.setLearnerVars()

        # Output a learner regardless of whether input data is provided
        self.learner = AZorngCvANN.CvANNLearner(\
                        nHidden = self.nHidden,\
                        stopCrit = self.stopCrit, \
                        maxIter = self.maxIter, \
                        eps = self.eps,
                        optAlg = self.optAlg, \
                        priors = self.priors, \
                        scaleClass = self.scaleClass,\
                        scaleData = self.scaleData)


        self.learner.name = self.name
        self.send("Learner", self.learner)
        self.createClassifier() 
    
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
            self.classifier = self.learner(self.data)
            self.classifier.name = self.name
        else:
            self.classifier = None
            self.data = None
        self.send("Classifier", self.classifier) 


if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWCvANN() 
    appl.setMainWidget(ow) 
    ow.show() 
    self.data = dataUtilities.DataTable('iris.tab') 
    ow.data(self.data) 
    appl.exec_loop()


