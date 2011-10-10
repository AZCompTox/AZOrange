"""
<name>Cv-SVM</name>
<description>SVM - OpenCV Support Vector Machine learner/classifier.</description>
<icon>icons/BasicSVM.png</icon>
<contact>Pedro Almeida</contact>
<priority>3</priority>
"""
from AZutilities import dataUtilities
import orange, OWGUI, sys
from OWWidget import *
from trainingMethods import AZorngCvSVM
from exceptions import SystemExit
import AZOrangeConfig as AZOC
import os

class OWCvSVM(OWWidget):
    #settingsList=["C","nu","p","probability","shrinking","gamma","degree", "coef0", "kernel_type", "name", "useNu"]
    def __init__(self, parent=None, signalManager=None, name="SVM"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs=[("Example Table", ExampleTable, self.setData)]
        self.outputs=[("Learner", orange.Learner),("Classifier", orange.Classifier)]   #  ,("Support Vectors", ExampleTable)]

        self.data = None
        self.name="Support Vector Machine"
        self.nMin = -1
        self.nMax = 1
        self.nClassMin = -1
        self.nClassMax = 1
        self.svm_typeGUI = None
        self.stopCritGUI = None
        self.priorsGUI = ""
        #Read default parameters from AZOrangeConfig.py file
        for par in ("kernel_type", "svm_type","gamma","C","p","epsC","epsR","stopCrit",\
                        "maxIter","nu","scaleData","scaleClass","priors","coef0","degree"):
            setattr(self, par, AZOC.CVSVMDEFAULTDICT[par])
        if self.svm_type in (103,104):
            self.eps = self.epsR
        else:
            self.eps = self.epsC
        self.setSVMTypeGUI(self.svm_type)
        self.setStopCritGUI(self.stopCrit)
        self.setPriorsGUI()#self.priors)

##scJS Added parameter needed for saving an SVM model
        self.modelFile = os.path.join(os.getcwd(),"SVM.model")
##ecJS
        OWGUI.lineEdit(self.controlArea, self, 'name', box='Learner/Classifier Name', \
                       tooltip='Name to be used by other widgets to identify your learner/classifier.<br>This should be a unique name!')
        ## scLC
        # Do not use the one class svm (type 2).
        self.svmTypeBox=b=OWGUI.widgetBox(self.controlArea,"SVM Type")
        self.svmTyperadio = OWGUI.radioButtonsInBox(b, self, "svm_typeGUI", btnLabels=["C-SVC", "nu-SVC","epsilon-SVR", "nu-SVR"], callback=self.setSVMTypeLearner)
        ## ecLC

        self.kernelBox=b=OWGUI.widgetBox(self.controlArea, "Kernel")
        self.kernelradio = OWGUI.radioButtonsInBox(b, self, "kernel_type", btnLabels=["Linear,   x.y", "Polynomial,   (g*x.y+c)^d",
                    "RBF,   exp(-g*(x-y).(x-y))", "Sigmoid,   tanh(g*x.y+c)"], callback=self.changeKernel)

        self.gcd = OWGUI.widgetBox(b, orientation="horizontal")
        self.leg = OWGUI.doubleSpin(self.gcd, self, "gamma",0.0,10.0,0.0001, label="gamma: ", orientation="horizontal", callback=self.changeKernel)
        self.led = OWGUI.doubleSpin(self.gcd, self, "coef0", 0.0,10.0,0.0001, label="  coef0: ", orientation="horizontal", callback=self.changeKernel)
        self.lec = OWGUI.doubleSpin(self.gcd, self, "degree", 0.0,10.0,0.5, label="  degree: ", orientation="horizontal", callback=self.changeKernel)

        OWGUI.separator(self.controlArea)


        #self.layout=QVBoxLayout(self.mainArea)
        self.optionsBox=b=OWGUI.widgetBox(self.mainArea, "Options", addSpace = True)
        self.cBox = OWGUI.doubleSpin(b,self, "C", 0.0, 512.0, 0.5, label="Model complexity (C)",  orientation="horizontal", tooltip="Penalty cost")
        self.cBox.control.setDecimals(3)
        self.pBox = OWGUI.doubleSpin(b,self, "p", 0.0, 10.0, 0.1, label="Tolerance (p)",  orientation="horizontal", tooltip="zero-loss width zone")
        self.pBox.control.setDecimals(3)
        self.nuBox = OWGUI.doubleSpin(b, self, "nu", 0.0,1.0,0.1, label="Complexity bound (nu)",  orientation="horizontal", tooltip="Upper bound on the ratio of support vectors")
        self.nuBox.control.setDecimals(3)        

        self.priorBox = OWGUI.lineEdit(b, self, 'priorsGUI', box='Class weights', \
                       tooltip='Ex: POS:2 , NEG:1')

        self.stopCritBox=OWGUI.widgetBox(self.optionsBox,"Stop Criteria")
        self.stopCritradio = OWGUI.radioButtonsInBox(self.stopCritBox, self, "stopCritGUI", btnLabels=["Maximum number of iterations (maxIter)", "Numeric precision (eps)"], callback=self.setStopCritLearner)
        OWGUI.doubleSpin(self.stopCritBox,self, "eps", 0.0, 0.5, 0.001, label="Numeric precision (eps)",  orientation="horizontal")
        OWGUI.doubleSpin(self.stopCritBox,self, "maxIter", 0, 100000, 1, label="Maximum number of iterations (maxIter)",  orientation="horizontal")


        self.scaleDataBox=OWGUI.checkBox(b, self, "scaleData", "Scale Data", tooltip="Scales the training data, and the input examples",callback = self.changeScaling)
        self.scaleClassBox=OWGUI.checkBox(OWGUI.indentedBox(b), self, "scaleClass", "Scale Class variable", tooltip="Scales the class variable of training data, and the output predictions",callback = self.changeScaling)
        #self.scalenMin = OWGUI.lineEdit(OWGUI.indentedBox(b), self, 'nMin', 'Var min:', orientation="horizontal", tooltip='Variables scaling minimum')
        #self.scalenMax = OWGUI.lineEdit(OWGUI.indentedBox(b), self, 'nMax', 'Var max:', orientation="horizontal", tooltip='Variables scaling maximum')
        #OWGUI.separator(OWGUI.indentedBox(b))
        self.scalenClassMin = OWGUI.lineEdit(OWGUI.indentedBox(b), self, 'nClassMin', 'Class min:', orientation="horizontal", tooltip='Class scaling minimum')
        self.scalenClassMax = OWGUI.lineEdit(OWGUI.indentedBox(b), self, 'nClassMax', 'Class max:', orientation="horizontal", tooltip='Class scaling maximum')
        OWGUI.separator(self.controlArea)

        #self.paramButton=OWGUI.button(self.controlArea, self, "Automatic parameter search", callback=self.parameterSearch,
        #                             tooltip="Automaticaly searches for parameters that optimize classifier acuracy")
        self.optionsBox.adjustSize()
        #self.layout.add(self.optionsBox)

        #OWGUI.separator(self.controlArea)

        OWGUI.button(self.controlArea, self,"&Apply", callback=self.applySettings)
##scJS Adding saving of svm models 

        # Get desired location of model file
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving Model", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "modelFile", labelWidth=80,  orientation = "horizontal", \
        tooltip = "Once a model is created (connect this widget with a data widget), \nit can be saved by giving a \
file name here and clicking the save button.")
        #L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        #button.setMaximumWidth(25)

        # Save model
        OWGUI.button(self.controlArea, self,"&Save SVM model", callback=self.saveModel)

        #self.adjustSize()
        #self.loadSettings()
        self.setStopCritLearner()
        self.setSVMTypeLearner()
        self.changeKernel() 
        self.applySettings()
        self.resize(500,400)

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


    def showEvent(self, ev):
        self.refreshParams()

    def refreshParams(self):
        if not self.learner:
            return
        self.setSVMTypeGUI(self.learner.svm_type)
        self.setStopCritGUI(self.learner.stopCrit)
        self.setPriorsGUI()#self.learner.priors)
        self.kernel_type = self.learner.kernel_type
        self.degree = self.learner.degree
        self.maxIter = self.learner.maxIter
        self.gamma = self.learner.gamma
        self.coef0 = self.learner.coef0
        self.C = self.learner.C
        self.p = self.learner.p
        self.eps = self.learner.eps
        self.nu = self.learner.nu
        self.scaleData = self.learner.scaleData
        self.scaleClass = self.learner.scaleClass
        self.changeScaling()
        self.adjustSize()

    def setPriorsGUI(self):#,learnerPriors): 
        #Set the GUI priors accordingly to the LEarner
        self.priorsGUI = self.priors and str(self.priors) or ""
        
    def setStopCritGUI(self,stopCrit):
        #Set the GUI stopCrit  accordingly to the Learner
        #    GUI:       Max=0   EPS=1
        #    Learner:   Max=1   EPS=2
        self.stopCritGUI = stopCrit-1

    def setStopCritLearner(self):
        #Set the self.stopCrit  accordingly to the GUI defined by user
        self.stopCrit = self.stopCritGUI+1

    
    def changeScaling(self):
        if not self.scaleData:
            self.scaleClassBox.setDisabled(1)
            #self.scalenMin.setDisabled(1)
            #self.scalenMax.setDisabled(1)
        else:
            self.scaleClassBox.setDisabled(0)
            #self.scalenMin.setDisabled(0)
            #self.scalenMax.setDisabled(0)
        if not self.scaleClass or not self.scaleData:
            self.scalenClassMin.setDisabled(1)
            self.scalenClassMax.setDisabled(1)
        else:
            self.scalenClassMin.setDisabled(0)
            self.scalenClassMax.setDisabled(0)

    def setSVMTypeGUI(self,svm_type):
        #Set the GUI svm_typeGUI  accordingly to the Learner
        if svm_type == 100:
            self.svm_typeGUI = 0
        elif svm_type == 101:
            self.svm_typeGUI = 1
        elif svm_type == 103:
            self.svm_typeGUI = 2
        elif svm_type == 104:
            self.svm_typeGUI = 3
        else:
            self.svm_typeGUI = 0
            svm_type = 100

    def setSVMTypeLearner(self):
        #Set the self.svm_type  accordingly to the GUI defined by user
        if self.svm_typeGUI == 0:
            svm_type = 100 #C_SVC
        elif self.svm_typeGUI == 1:
            svm_type = 101 #NU_SVC
        elif self.svm_typeGUI == 2:
            svm_type = 103 #EPS_SVR
        elif self.svm_typeGUI == 3:
            svm_type = 104 #NU_SVR
        else:
            svm_type = 100
            self.svm_typeGUI = 0
        #enable/disable parameters specific
        if svm_type==100:
            for a,b in zip([self.cBox, self.pBox, self.nuBox, self.priorBox.box], [0,1,1,0]):
                a.setDisabled(b)
        elif svm_type==101:
            for a,b in zip([self.cBox, self.pBox, self.nuBox, self.priorBox.box], [1,1,0,0]):
                a.setDisabled(b)
        elif svm_type==103:
            for a,b in zip([self.cBox, self.pBox, self.nuBox, self.priorBox.box], [0,0,1,0]):
                a.setDisabled(b)
        elif svm_type==104:
            for a,b in zip([self.cBox, self.pBox, self.nuBox, self.priorBox.box], [0,1,0,0]):
                a.setDisabled(b)


    def changeKernel(self):
        if self.kernel_type==0:
            for a,b in zip([self.leg, self.led, self.lec], [1,1,1]):
                a.setDisabled(b)
        elif self.kernel_type==1:
            for a,b in zip([self.leg, self.led, self.lec], [0,0,0]):
                a.setDisabled(b)
        elif self.kernel_type==2:
            for a,b in zip([self.leg, self.led, self.lec], [0,1,1]):
                a.setDisabled(b)
        elif self.kernel_type==3:
            for a,b in zip([self.leg, self.led, self.lec], [0,0,1]):
                a.setDisabled(b)


    def setData(self, data=None):
        self.error(0)
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
        self.applySettings()

    def applySettings(self):
        self.error(0)
        self.warning(0)
        if self.priorsGUI:
            self.priors = str(self.priorsGUI)
        else:
            self.priors = None

        #Update the variables handling the learner specific values
        self.setSVMTypeLearner()
        self.setStopCritLearner()
        self.changeKernel()
        #Create the Learner
        self.learner=AZorngCvSVM.CvSVMLearner()
        #Set the several learner parameters
        self.learner.priors = self.priors
        for attr in ("nMin","nMax","nClassMin","nClassMax","scaleClass","scaleData"):
            setattr(self.learner, attr, getattr(self, attr))

        for attr in ("name", "kernel_type", "degree", "svm_type", "stopCrit", "maxIter"):
            setattr(self.learner, attr, getattr(self, attr))

        for attr in ("gamma", "coef0", "C", "p", "eps", "nu"):
            setattr(self.learner, attr, float(getattr(self, attr)))

        self.classifier=None
        self.supportVectors=None
        #Create the classifier
        if self.data:
            self.classifier=self.learner(self.data)
        if self.learner:
            self.learner.name=str(self.name)
            self.send("Learner", self.learner)
        if self.classifier:
            self.classifier.name=str(self.name)
            #self.supportVectors=self.classifier.classifier.supportVectors
            self.send("Classifier", self.classifier)
        #if self.supportVectors:
        #    self.send("Support Vectors", self.supportVectors)
        self.refreshParams()


    def saveModel(self):
        """Write a SVM classifier instance to disk."""

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


from exceptions import Exception        
class UnhandledException(Exception):
    pass

import sys
if __name__=="__main__":
    app=QApplication(sys.argv)
    w=OWSVM()
    app.setMainWidget(w)
    w.show()
    d=dataUtilities.DataTable("../../doc/datasets/iris.tab")
    w.setData(d)
    app.exec_loop()
    w.saveSettings()
