"""
<name>Consensus</name>
<description>Generalized Consensus Models. Combines learners/classifiers into majority or average consensus models. </description>
<icon>icons/consensus.png</icon>
<contact>Pedro Almeida</contact>
<priority>6</priority>
"""
import string,os

from OWWidget import *
import OWGUI
import orange
from trainingMethods import AZorngConsensus
import AZOrangeConfig as AZOC
import orngImpute
from AZutilities import dataUtilities
from AZutilities import miscUtilities



class Classifier:
    def __init__(self, classifier, id):
        classifier.id = id
        self.classifier = classifier
        self.name = classifier.name
        self.id = id
        self.time = time.time() # Time stamp

class Learner:
    def __init__(self, learner, id):
        learner.id = id
        self.learner = learner
        self.name = learner.name
        self.id = id
        self.time = time.time() # Time stamp


class OWConsensus(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Consensus'):
        OWWidget.__init__(self, parent, signalManager, name)

        # Define the input and output channels
        self.inputs = [("Examples", ExampleTable, self.setData), ("Classifiers", orange.Classifier, self.setClassifiers, Multiple), ("Learners", orange.Learner, self.setLearners, Multiple)]
        self.outputs = [("Learner", orange.Learner),("Classifier", orange.Classifier)]

        self.name = name
        self.expr = ""
        self.data = None
        self.classifiers = {}
        self.learners = {}
        self.classifier = None
        self.learner = None
        # Define the parameters readable through the GUI and set their default values
        self.modelFile = os.path.join(os.getcwd(),"Consensus.model")

        self.defineGUI()
        self.applySettings()

    def setWindowTitle(self, caption):
        """ This is an override of the function setWindowTitle in OWBaseWidget.py 
            It is used to update the name of the Learner with a suffix corresponding to the number given by DOC 
            orngDoc.py to this Learner widget when several widgets like this are present in the Canvas.
            This way, when signals are sent, and the destination widget uses the name property, there will 
            be different names identifying each one.
        """
        # First call the base class function in order to the captionTitle is updated with a numbered suffix
        OWWidget.setWindowTitle(self, caption)
        suffix = str(self.windowTitle().split(" ")[-1]).strip()
        if len(suffix)>=3 and suffix[0]=="(" and suffix[-1]==")":
            self.name = self.name + " " + suffix
            self.applySettings() 

    def defineGUI(self):

        # Set the number of hidden neurons throught the GUI
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Learner/Classifier Name', \
                       tooltip='Name to be used by other widgets to identify your learner/classifier.<br>This should be a unique name!')

        info = OWGUI.widgetBox(self.controlArea, "Constituting Consensus Learners", addSpace = True, orientation=0)
        self.learnerInfo = OWGUI.widgetLabel(info, '')

        OWGUI.lineEdit(self.controlArea, self, 'expr', box='Expression', \
                       tooltip='Expression')

        # Apply the settings and send the learner to the output channel
        OWGUI.button(self.controlArea, self,"&Apply settings", callback=self.applySettings)

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

        # Status Label
        #box = OWGUI.widgetBox(self.controlArea, "Info ")
        self.info = OWGUI.widgetLabel(self.controlArea, '')
 

        self.applySettings()
        self.adjustSize()
      
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
        """Write an ann classifier instance to disk """
        if self.classifier and self.modelFile:
            if not self.classifier.write(self.modelFile):
                print "ERROR: model was NOT saved!"
            else:
                print "Saved Consensus model to ",self.modelFile
        else:
            print "ERROR: Something is missing:"
            print "   classifiers:",self.classifiers
            print "   learners:",self.learners
            print "   modelFile:",self.modelFile


    def setLearners(self, learner, id=None):
        """add/remove a Learner"""
        if len(self.classifiers):
            self.warning("To use Learners, please disconmnect first all the classifiers.")
            return
        if learner: # a new or updated learner
            if id in self.learners: # updated learner
                time = self.learners[id].time
                self.learners[id] = Learner(learner, id)
                self.learners[id].time = time
            else: # new learner
                self.learners[id] = Learner(learner, id)
        else: # remove a learner 
            if id in self.learners:
                del self.learners[id]
        self.applySettings()
 


    def setClassifiers(self, classifier, id=None):
        """add/remove a Classifier"""
        if len(self.learners):
            self.warning("To use Classifiers, please disconmnect first all the Learners.")
            return
        if classifier: # a new or updated classifier
            if id in self.classifiers: # updated classifier
                time = self.classifiers[id].time
                self.classifiers[id] = Classifier(classifier, id)
                self.classifiers[id].time = time
            else: # new classifier
                self.classifiers[id] = Classifier(classifier, id)
        else: # remove a classifier 
            if id in self.classifiers:
                del self.classifiers[id]
        self.applySettings()
 
      
    def sendObjects(self):
        self.send("Learner", self.learner)
        self.send("Classifier", self.classifier)

    def applySettings(self):
        """Create the learner with selected settings and send to the output channel. """
        self.error(0)
        self.warning(0)
        self.createLearner()
        self.createClassifier() 
        self.sendObjects()
    
    def setData(self,data):
        self.error(0)
        self.warning(0)
        if data:
            self.data=data
        else:
            self.data=None
        self.applySettings()

    def createLearner(self):
        # Output a learner regardless of whether input data is provided
        myLearnersText = ""
        if len(self.learners) >= 2:
            # Get names of learners used in the consensus
            for idx,l in enumerate(self.learners.values()):
                myLearnersText = myLearnersText + l.name+"\n"
            self.learnerInfo.setText(myLearnersText)
            # Assure that the expression uses the defined learner names
            expression = str(self.expr.replace('\x00','')).strip().split(",")
            if len(expression) <= 1:
                expression = str(self.expr.replace('\x00','')).strip()
            else:
                expression = [e.strip() for e in expression]

            if expression:
                learners = {}
                for l in self.learners.values():
                    learners[str(l.name)] = l.learner
                self.learner = AZorngConsensus.ConsensusLearner(learners = learners, expression = expression)

            else:
                learners = [l.learner for l in self.learners.values()]
                self.learner = AZorngConsensus.ConsensusLearner(learners = [l.learner for l in self.learners.values()])
            if not self.learner:
                self.error("Could not build a Consensus learner. Please check the outpou window for more details")
                self.learner = None
                return
            self.learner.name = str(self.name)
        else:
            self.learner = None

    def createClassifier(self): 
        """ Output a classifier. Set to None if there is no input data.  """
        if self.data and self.learner: 
            self.classifier = self.learner(self.data)
            self.classifier.name = str(self.name)
            self.info.setText(self.classifier.status)
        elif len(self.classifiers) >= 2:
            if self.data:
                self.warning("When using classifiers as inputs, data must be disconnected!")
                self.classifier = None
            else:
                expression = str(self.expr.replace('\x00','')).strip().split(",")
                if len(expression) <= 1:
                    expression = str(self.expr.replace('\x00','')).strip()
                else:
                    expression = [e.strip() for e in expression]
                if expression:
                    classifiers = {}
                    for idx,l in enumerate(self.classifiers.values()):
                        classifiers["Classifier_"+str(idx)] = l
                    self.classifier = AZorngConsensus.ConsensusClassifier(classifiers = classifiers, expression = expression)

                else:
                    self.classifier = AZorngConsensus.ConsensusClassifier(classifiers = [c.classifier for c in self.classifiers.values()], expression = None)
                self.classifier.name = str(self.name)
                self.info.setText(self.classifier.status)
            if not self.classifier:
                self.error("Could not build a Consensus Classifier. Please check the outpou window for more details")
                self.classifier = None
                return
        else:
            self.classifier = None

    def showEvent(self, ev):
        if self.classifier and hasattr(self.classifier,"status"):
            self.info.setText(self.classifier.status)
        

if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWConsensus() 
    appl.setMainWidget(ow) 
    ow.show() 
    self.data = dataUtilities.DataTable('iris.tab') 
    ow.data(self.data) 
    appl.exec_loop()


