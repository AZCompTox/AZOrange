"""
<name>Cv-Bayes</name>
<description>OpenCV Naive/Normal Bayes learner/classifier.</description>
<icon>icons/NaiveBayes.png</icon>
<contact>Pedro Almeida</contact>
<priority>5</priority>
"""
from OWWidget import *
import OWGUI
from exceptions import Exception
from trainingMethods import AZorngCvBayes
import AZOrangeConfig as AZOC

class OWCvBayes(OWWidget):
    parameters = ["scale"]
    def __init__(self, parent=None, signalManager = None, name='kNN'):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea = 0, resizingEnabled = 0)

        self.callbackDeposit = []

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Learner", orange.Learner),("Classifier", orange.Classifier)]


        # Settings
        self.name = 'CvBayes'
        self.classifier = None
        self.learner = None
        self.modelFile = os.path.join(os.getcwd(),"Bayes.model")
        for par in self.parameters:
            setattr(self, par, AZOC.CVBAYESDEFAULTDICT[par])


        self.data = None                    # input data set
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Learner/Classifier Name', \
                 tooltip='Name to be used by other widgets to identify the learner/classifier.')

        OWGUI.separator(self.controlArea)

        pars = OWGUI.widgetBox(self.controlArea, "Parameters")
        OWGUI.checkBox(pars, self, "scale", "Scale Data")
        

        OWGUI.separator(self.controlArea)

        self.setLearner()
        OWGUI.button(self.controlArea, self, "&Apply", callback = self.setLearner, disabled=0)
        
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
        """Write a Bayes classifier instance to disk """
        self.warning(0)
        if self.classifier and self.modelFile and self.data:
            if not self.classifier.write(str(self.modelFile)):
                print "ERROR: model was NOT saved!"
                self.warning("Cannot save model. Please check the output window.")
            else:
                print "Saved CvBayes type model to ",self.modelFile
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
                            [("Scale Data", OWGUI.YesNo[self.scale]),
                             ])
        self.reportData(self.data)
        
            
    def setData(self,data):
        self.data = self.isDataWithClass(data, orange.VarTypes.Discrete) and data or None
        self.setLearner()

    def refreshParams(self):
        if self.learner:
            for par in self.parameters:
                setattr(self, par, getattr(self.learner,par))
        #self.printLearnerPars()

    def printLearnerPars(self):
        for par in self.parameters:
            print par,":",getattr(self, par)


    def setLearner(self):
        self.error(0)
        self.warning(0)
        #self.printLearnerPars()
        self.learner = AZorngCvBayes.CvBayesLearner(scale = self.scale)
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
    ow=OWCvBayes()

##    dataset = orange.ExampleTable('adult_sample')
##    ow.setData(dataset)

    ow.show()
    a.exec_()
    ow.saveSettings()
