from AZutilities import dataUtilities
"""
<name>Similarity Descriptors</name>
<description>Calculates similarity descriptors</description>
<icon>icons/SimDesc.png</icon>
<contact>Pedro Rafael Almeida</contact>
<priority>14</priority>
"""
import string,time

from OWWidget import *
import OWGUI
import orange
import AZOrangeConfig as AZOC
from AZutilities import SimBoostedQSAR


class OWSimBoostedQSAR(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='Trained learner'):
        OWWidget.__init__(self, parent, signalManager, name)

        # Define the input and output channels
        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples", ExampleTable)]

        self.methods = SimBoostedQSAR.methods
        self.methodsList = SimBoostedQSAR.methods.keys()
        self.selectedMethods = []

        self.last_pDone = 0
        self.startTime = 0
        self.dT_ppBuffer = []

        self.outData = None
	self.dataset = None
        self.actives = None
        self.activesFile = os.getcwd()
        self.defineGUI()

    def calcDesc(self):
        self.outData = None
        self.error()
        self.warning()
        if not self.actives:
            self.warning("No Actives defined!")
        elif not self.dataset:
            self.warning("No input Dataset!")
        elif not self.selectedMethods:
            self.warning("No Methods Selected!")
        else:
            self.progressSteps = 100
            self.progress = QProgressDialog("Calculating similarity descriptors", "Cancel", 0, self.progressSteps , None, Qt.Dialog )
            self.progress.setWindowModality(Qt.WindowModal)
            self.progress.setMinimumDuration(0)
            self.progress.forceShow()
            self.progress.setValue(0)
        
            #set the selected methods
            methods = [SimBoostedQSAR.methods[self.methodsList[mIdx]] for mIdx in self.selectedMethods]

            self.startTime = time.time()
            self.dT_ppBuffer = []
            self.last_pDone = 0
            self.outData = SimBoostedQSAR.getSimDescriptors(self.actives, self.dataset, methods, callBack = self.advance)
                
            self.progress.close()
            if not self.outData or not  len(self.outData):
                self.error("Could not calculate descriptors. Please check the output window for more details.")
                self.outData = None 
        self.send("Examples", self.outData)

    def advance(self, pDone):
        now = time.time()
        LowPassFilterBuffer = 10
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


                self.progress.setLabelText("Calculating similarity descriptors\nEstimated time left: "+strEstTime)
        self.progress.setValue(pDone) 
        self.last_pDone = pDone
        return True

    def loadActives(self):
        self.error()
        self.warning()
        if os.path.isfile(self.activesFile):
            actives = dataUtilities.loadSMI(self.activesFile)
            self.actives = []
            if not actives or len(actives) < 1 or "SMILES" not in actives.domain:
                self.error("Bad formatted file for Actives")
            else:
                for ex in actives:
                    self.actives.append(str(ex["SMILES"].value))
        else:
            self.error("File '"+str(self.activesFile)+"' was not found!")     
            self.actives = None

    def browseFile(self):
        self.error()
        self.warning()
        # Possible modes:
        var = os.path.realpath(str(self.activesFile))
        if os.path.isdir(var):
            startfile = var
        elif os.path.isdir(os.path.split(var)[0]):
            startfile = os.path.split(var)[0]
        else:
            startfile=os.getcwd()
        filename = str(QFileDialog.getOpenFileName(self,"Open Data File",startfile, 'Smiles files (*.smi)'))
        if filename:
            self.activesFile = str(filename)
            self.loadActives() 

    def destroy(self, dw = 1, dsw = 1):
	self.linksOut.clear()
	if self.dataset:
	    del self.dataset


    def defineGUI(self):
        # methods Selection
        self.methodsListBox = OWGUI.listBox(self.controlArea, self, "selectedMethods", "methodsList", box = "Select Chemical Similarity Methods", selectionMode = QListWidget.ExtendedSelection) 
        # Set location of model file
        boxFile = OWGUI.widgetBox(self.controlArea, "Path for loading Actives", addSpace = True, orientation=0)
        L1 = OWGUI.lineEdit(boxFile, self, "activesFile", labelWidth=80,  orientation = "horizontal", tooltip = "Actives must be specified in a .smi file")
        L1.setMinimumWidth(200)
        button = OWGUI.button(boxFile, self, '...', callback = self.browseFile, disabled=0,tooltip = "Browse for a file to load the Actives")
        button.setMaximumWidth(25)

        # Load Actives
        OWGUI.button(boxFile, self,"Load Actives", callback=self.loadActives)

        # Start descriptors calculation
        OWGUI.button(self.controlArea, self,"Calculate descriptors", callback=self.calcDesc)

        self.adjustSize()


    def setData(self, dataset): 
        self.error()
        self.warning()
        if dataset:
            self.dataset = dataset
            self.calcDesc()
        else:
            self.dataset = None
            self.send("Examples", None) 


if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWSimBoostedQSAR() 
    appl.setMainWidget(ow) 
    ow.show() 
    dataset = dataUtilities.DataTable('iris.tab') 
    ow.data(dataset) 
    appl.exec_loop()


