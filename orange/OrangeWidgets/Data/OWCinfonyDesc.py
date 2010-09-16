"""
<name>Cinfony Descriptors</name>
<description>Adds any Cinfony descriptors to a data set. Requires an attribute with smiles ([SMILES, Molecule SMILES, Compound Structure]) </description>
<icon>icons/cinfony.png</icon> Cinfony Descriptors
<contact>Pedro Almeida</contact> 
<priority>13</priority>
"""

#
# OWCinfonyDesc.py
# The File Widget
# A widget to get descriptors from Cinfony
#
from AZutilities import dataUtilities
from OWWidget import *
import OWGUI,OWGUIEx
from AZutilities.AZOWGUI import widgetSelector
import os
import sys
import string
from AZutilities import getCinfonyDesc

class OWCinfonyDesc(OWWidget):
    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "CinfonyDesc",wantMainArea = 0)

        self.data = None
        #self.getCinfonyDescObj = getCinfonyDesc.CinfonyDescGetter()
       
        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples", ExampleTable)]

        # GUI definition
        self.Descriptors = widgetSelector(self.controlArea, "Descriptors",applyText = "Retrieve", callbackOnApply=self.__retrieve)
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info = OWGUI.widgetLabel(box, '')

        self.__initListBox()



    def setData(self, data):
        self.data = data
        self.error(0)
        self.warning(0)
        if not data:
            self.send("Examples", None)


    def __retrieve(self,descList):
        self.error(0)
        self.warning(0)
        resultsData = None
        if not self.data:
            self.warning("Missing input data!")
            self.send("Examples",None)
            return
        else:
            progressSteps = 5 
            progress1 = QProgressDialog("Calculating selected descriptors with Cinfony", "Cancel", 0, progressSteps , None, Qt.Dialog )
            progress1.setWindowModality(Qt.WindowModal)
            progress1.setMinimumDuration(0)
            progress1.forceShow()
            progress1.setValue(1)
            time.sleep(0.1)
            progress1.setValue(2)
            try: 
                resultsData = getCinfonyDesc.getCinfonyDescResults(self.data , descList)
                progress1.setValue(3)
            except:
                progress1.close()
                self.error("Error calculating the descriptors with Cinfony.\n")
                self.send("Examples",None)
                return
            time.sleep(0.1)
            progress1.setValue(4)
            time.sleep(0.1)
            progress1.setValue(5)
            progress1.close()

        if resultsData:
            self.send("Examples",resultsData)
        else:
            self.warning("Nothing retrieved from Cinfony.")
            self.send("Examples",None)



    def __initListBox(self):
        ## Get available descriptors from Cinfony and show progress bar. 
        progressSteps = 3 
        progress1 = QProgressDialog("Retrieving available descriptors from Cinfony.", "Cancel", 0, progressSteps , None, Qt.Dialog )
        progress1.setWindowModality(Qt.WindowModal)
        progress1.setMinimumDuration(0)
        progress1.forceShow()
        progress1.setValue(1)
        time.sleep(0.1)
        progress1.setValue(2)
        try: 
            availableDescriptors = getCinfonyDesc.getAvailableDescs() 
        except:
            progress1.close()
            #self.updateInfo()
            self.error("Error retrieving descriptors from Cinfony.\n")
            return
        progress1.setValue(3)
        progress1.close()

        # Create a sorted list of descriptors available in Cinfony
        self.Descriptors.setInputItems(availableDescriptors)
        #toolkitsEnabled
        if getCinfonyDesc.toolkitsEnabled:
            self.info.setText("Toolkits available for retrieving descriptors: "+ str(getCinfonyDesc.toolkitsEnabled)[1:-1].replace("'",""))
        else:
            self.error("No available toolkits for retrieving descriptors!")
            self.info.setText("Cinfony returned no available scripts!")


if __name__ == "__main__":
    a=QApplication(sys.argv)
    owf=OWIBIS()
    a.setMainWidget(owf)
    owf.show()
    a.exec_loop()
    owf.saveSettings()
