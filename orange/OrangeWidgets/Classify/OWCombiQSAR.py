from AZutilities import dataUtilities
from AZutilities import miscUtilities
"""
<name>Combi-QSAR</name>
<description>Widget to automatically select the most accurate AZOrange machine learning algorithm (with optimized model hyper-parameters) and to calculated an unbiased assessemnt of the expected generalization accuracy of the resulting QSAR model. For a detailed description, please see the Poster; http://svn.seml.astrazeneca.net/trac/CC-AZ-Orange/browser/trunk/doc/AZOrangeOpenTox.pdf </description>
<icon>icons/CombiQSAR.png</icon>
<contact>Pedro Rafael Almeida</contact>
<priority>8</priority>
"""
import string,time

from OWWidget import *
import OWGUI
import orange
import AZOrangeConfig as AZOC
from AZutilities import competitiveWorkflow

#DEBUG
#import pickle
#from trainingMethods import AZBaseClasses

class OWCombiQSAR(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Combi-QSAR model'):
        OWWidget.__init__(self, parent, signalManager, name)

        # Define the input and output channels
        self.inputs = [("Classified Examples", ExampleTable, self.setData)]
        self.outputs = [("Classifier", orange.Classifier), ("Examples", ExampleTable)]

        self.queueTypes = ["NoSGE","batch.q","quick.q"] 
        self.outputModes = ["Statistics for all available algorithms. Please note, no model selection.", "Model and statistics (unbiased wrt model selection)"]
        self.name = name
	self.dataset = None
	self.classifier = None
        self.statInfo = ""
        self.statistics = None
        self.queueType = 0
        self.outputSel = 0
        self.isClassDiscrete = None

        #Paths
        self.modelFile = ""
        self.statPath = ""
        self.lastPath = os.getcwd()
       
        #progressBar
        self.last_pDone = 0
        self.startTime = 0
        self.dT_ppBuffer = []

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

        

    def saveStat(self, file2save = None):
        self.warning(0)
        self.error(0)
        if not self.statistics:
            self.warning("No statisitcs to save yet!")
            return
        if file2save:
            filename = file2save
        else:
            filename = self.browseFile(mode = "file")
        print filename 
        if filename:
            self.statistics.save(filename)

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

    def setStatSaveDir(self):
        self.statPath = self.browseFile(mode = "dir")

    def setModelSavePath(self):
        self.modelFile = self.browseFile(mode = "file")

    def browseFile(self, var="", mode = "file"):
        # Possible modes:
        #dir   -  Browse for save  dir
        #file  -  Browse for save file
        if var:
            var = os.path.realpath(str(var))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        elif os.path.isdir(self.lastPath):
            startfile = self.lastPath
        elif os.path.isdir(os.path.split(self.lastPath)[0]):
            startfile = os.path.split(self.lastPath)[0]
        else:
            startfile=os.getcwd()

        fileDialog = QFileDialog(self)

        if mode == "dir":
            filename = fileDialog.getExistingDirectory(self,"Select save location",startfile,QFileDialog.ShowDirsOnly)
        else: #file
            filename = str(fileDialog.getSaveFileName(self,"Select save Path",startfile))
        if str(filename):
            self.lastPath = str(filename)
        return str(filename)


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
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving the results", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "statPath", labelWidth=80,  orientation = "horizontal", tooltip = "Please use full path to results file to be created.")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.setStatSaveDir, disabled=0,tooltip = "Choose the dir where to save.")
        button.setMaximumWidth(25)

        

        # Apply the settings and send the model to the output channel
        OWGUI.button(self.controlArea, self,"&Apply settings", callback=self.Apply)
 
        # Set location of model file
        self.boxFile = OWGUI.widgetBox(self.controlArea, "Path for saving Model", addSpace = True, orientation=0)
        self.L1 = OWGUI.lineEdit(self.boxFile, self, "modelFile", labelWidth=80,  orientation = "horizontal", tooltip = "Once a model is created (connect this widget with a data widget), \nit can be saved by giving a file name here and clicking the save button.")
        self.L1.setMinimumWidth(200)
        self.button = OWGUI.button(self.boxFile, self, '...', callback = self.setModelSavePath, disabled=0,tooltip = "Choose the dir where to save. After chosen, add a name for the model file!")
        self.button.setMaximumWidth(25)

        # Save the model
        self.saveButton = OWGUI.button(self.controlArea, self,"&Save model", callback=self.saveModel)

        
        # Statistics show
        statBox = OWGUI.widgetBox(self.mainArea, "Statistics", addSpace = True, orientation=0)
        self.statInfoBox = OWGUI.widgetLabel(statBox, '')
        
        # Save the model
        OWGUI.button(self.mainArea, self,"&Save statistics", callback=self.saveStat)




        self.changeOutputMode()
	self.adjustSize()



    def changeOutputMode(self):
        self.warning(0)
        self.error(0)
        if self.outputSel == 0:
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

        progressSteps = 100
        self.progress = QProgressDialog("Running Combi-QSAR\nThis may take a while. Please wait....", "Cancel", 0, progressSteps , None, Qt.Dialog )
        self.progress.setWindowModality(Qt.WindowModal)
        self.progress.setMinimumDuration(0)
        self.progress.forceShow()
        self.progress.setValue(0)
 
        self.startTime = time.time()
        self.dT_ppBuffer = []
        self.last_pDone = 0

        if self.outputSel != 0:
            self.getModel()
        else:   
            self.getStatistics()

        self.progress.close()
        if self.statistics and os.path.isdir(str(self.statPath)): 
            self.saveStat(os.path.join(str(self.statPath),"statistics.txt"))


    def advance(self, pDone):
        now = time.time()
        LowPassFilterBuffer = 1
        if self.progress.wasCanceled():
            return False
        if pDone > self.last_pDone:
            dT_pp = (now - self.startTime)/((pDone - self.last_pDone)  * 1.0)
            self.dT_ppBuffer.append(dT_pp)
            self.startTime = now
            if len(self.dT_ppBuffer) >= LowPassFilterBuffer:
                estTime = (sum(self.dT_ppBuffer[-LowPassFilterBuffer:])/(1.0 * LowPassFilterBuffer)) * (100.0 - pDone)
                if estTime < 120:# < 2 min, count in sec
                    strEstTime = str(int(round(estTime)))+" sec."
                elif estTime < 7200: # 2 Hours, count in min
                    strEstTime = str(int(round(estTime/60)))+" min."
                elif estTime < 172800: # 2 Days, count in hours
                    strEstTime = str(round(estTime/3600,1))+" hours"
                else: #count in days
                    strEstTime = str(round(estTime/86400,1))+" days"

                self.progress.setLabelText("Running Combi-QSAR\nEstimated time left: "+strEstTime)
        self.progress.setValue(pDone)
        self.last_pDone = pDone
        return True

    def createStatData(self,statistics): 
        specialVars = [orange.StringVariable("Method"), orange.FloatVariable("Fold")]
        classificationVars = [   
                        orange.FloatVariable("CA"),
                        orange.FloatVariable("MCC"),
                        orange.FloatVariable("truePOS"),
                        orange.FloatVariable("trueNEG"),
                        orange.FloatVariable("falsePOS"),
                        orange.FloatVariable("falseNEG")  ]

        regressionVars = [   
                        orange.FloatVariable("Q2"),
                        orange.FloatVariable("RMSE")  ]
        if self.isClassDiscrete:
            allVars = specialVars + classificationVars 
        else:
            allVars = specialVars + regressionVars
        commVars = [  orange.FloatVariable("nTest"),
                       orange.FloatVariable("nTrain")  ]
        allVars += commVars

        self.statistics = orange.ExampleTable(orange.Domain(allVars,orange.FloatVariable("Stability")))
        for ml in statistics:
            # Total row 
            ex = orange.Example(self.statistics.domain)
            if ml == "selectedML": 
                ex["Method"] = "Total"
            else:
                ex["Method"] = ml+" Total"
            ex["Stability"] = statistics[ml]["StabilityValue"]
            # [[TP, FN],[FP, TN]]
            if self.isClassDiscrete:
                ex["CA"] = statistics[ml]["CA"] 
                ex["MCC"] = statistics[ml]["MCC"] 
                ex["truePOS"] = statistics[ml]["CM"][0][0] 
                ex["trueNEG"] = statistics[ml]["CM"][1][1]
                ex["falsePOS"] = statistics[ml]["CM"][1][0]
                ex["falseNEG"] = statistics[ml]["CM"][0][1] 
            else:
                ex["Q2"] = statistics[ml]["Q2"]
                ex["RMSE"] = statistics[ml]["RMSE"]
            self.statistics.append(ex)
            # Fold rows
            for fold,nTest in enumerate(statistics[ml]["foldStat"]["nTestCmpds"]):
                ex = orange.Example(self.statistics.domain)
                if ml == "selectedML":
                    ex["Method"] = statistics[ml]["foldStat"]["foldSelectedML"][fold]
                else:    
                    ex["Method"] = ml
                ex["Fold"] = fold
                ex["nTrain"] = statistics[ml]["foldStat"]["nTrainCmpds"][fold]
                ex["nTest"] = statistics[ml]["foldStat"]["nTestCmpds"][fold]
                # [[TP, FN],[FP, TN]]
                if self.isClassDiscrete:
                    ex["CA"] = statistics[ml]["foldStat"]["CA"][fold]
                    ex["MCC"] = statistics[ml]["foldStat"]["MCC"][fold]
                    ex["truePOS"] = statistics[ml]["foldStat"]["CM"][fold][0][0]
                    ex["trueNEG"] = statistics[ml]["foldStat"]["CM"][fold][1][1]
                    ex["falsePOS"] = statistics[ml]["foldStat"]["CM"][fold][1][0]
                    ex["falseNEG"] = statistics[ml]["foldStat"]["CM"][fold][0][1]
                else:
                    ex["Q2"] = statistics[ml]["foldStat"]["Q2"][fold]
                    ex["RMSE"] = statistics[ml]["foldStat"]["RMSE"][fold]
                self.statistics.append(ex)
        return self.statistics



    def getModel(self):
        self.warning(0)
        self.error(0)
        self.statInfo = ""
        self.statistics = None
        self.statInfoBox.setText("")
        """Get the model and send it to the output"""
        if not os.path.isdir(str(self.statPath)):
            statPath = None
            modelPath = None
        else:
            statPath = os.path.join(str(self.statPath),"statistics.pkl")
            modelPath = os.path.join(str(self.statPath),"Model")
        
        res = competitiveWorkflow.competitiveWorkflow(self.dataset, modelSavePath = modelPath, statisticsSavePath = statPath, runningDir = AZOC.NFS_SCRATCHDIR, queueType = self.queueTypes[self.queueType], callBack = self.advance)
        if not res:
            self.error("Errors occurred. Please check the output window.")
            self.send("Classifier", None)
            self.send("Examples", None)
            return
        self.classifier = res["model"][res["model"].keys()[0]]
        statistics = res["statistics"] 

        # DEBUG
        #fileh = open("/home/palmeida/dev/AZOrange/orange/OrangeWidgets/Classify/model.pkl")
        #statistics = pickle.load(fileh)
        #fileh.close()
        #self.classifier = AZBaseClasses.modelRead("/home/palmeida/dev/AZOrange/orange/OrangeWidgets/Classify/Model")

        if not statistics:
            self.statInfo = "Some error occured. Could not get statistics.\nPlease check the output window"
        else:
            self.statistics = self.createStatData(statistics)
            if statPath and os.path.isfile(statPath):
                self.statInfo = "Statistics were saved to " + statPath+"\n"+\
                                  "You can save the statistics in other place by using \n"
            else:
                self.statInfo +="You can save the statistics by using \n"
        self.statInfo +=    " the button 'Save statistics'\n\n"+\
                              "You can also use or view the statistics by connecting \n"+\
                              " the appropriate widget to this widget output"
        self.statInfoBox.setText(self.statInfo)

        self.classifier.name = str(self.name)
        self.send("Classifier", self.classifier)
        self.send("Examples", self.statistics)


    def getStatistics(self):
        self.warning(0)
        self.error(0)
        self.statInfo = ""
        self.statistics = None
        self.statInfoBox.setText("")
        """Get the statistics and save to desired place if specifyes"""
        if not os.path.isdir(str(self.statPath)):
            statPath = None
        else:
            statPath =  os.path.join(str(self.statPath),"statistics.pkl")
        runPath = miscUtilities.createScratchDir(desc = "CombiQSAR", baseDir = AZOC.NFS_SCRATCHDIR)


        # DEBUG
        #fileh = open("/home/palmeida/dev/AZOrange/orange/OrangeWidgets/Classify/stat.pkl")
        #statistics = pickle.load(fileh)
        #fileh.close()


        statistics = competitiveWorkflow.getMLStatistics(self.dataset, savePath = statPath, queueType = self.queueTypes[self.queueType], verbose = 0, logFile = None, callBack = self.advance)


        if not statistics:
            self.statInfo = "Some error occured. Please check the output window"
        else:
            self.statistics = self.createStatData(statistics)
            if statPath and os.path.isfile(statPath):
                self.statInfo = "Statistics were saved to " + statPath+"\n"+\
                                  "You can save the statistics in other place by using \n"
            else:
                self.statInfo +="You can save the statistics by using \n"
        self.statInfo +=    " the button 'Save statistics'\n\n"+\
                              "You can also use or view the statistics by connecting \n"+\
                              " the appropriate widget to this widget output" 
        self.statInfoBox.setText(self.statInfo)
        self.send("Classifier", None)
        self.send("Examples", self.statistics)


    def setData(self, dataset): 
        self.warning(0)
        self.isClassDiscrete = None
        self.dataset = None
        if dataset:
            self.dataset = dataset
            self.isClassDiscrete = dataset.domain.classVar.varType == orange.VarTypes.Discrete            
            self.Apply()


if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWCombiQSAR() 
    appl.setMainWidget(ow) 
    ow.show() 
    dataset = dataUtilities.DataTable('iris.tab') 
    ow.data(dataset) 
    appl.exec_loop()


