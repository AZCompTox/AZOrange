"""
<name>Merge Data</name>
<description>Merge datasets horizontally based on values of selected attributes.</description>
<icon>icons/MergeData.png</icon>
<priority>1110</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si), changed Ny Pedro Almeida</contact>
"""

import orange
from OWWidget import *
import OWGUI
from AZutilities import dataUtilities


class OWMergeData(OWWidget):

##    settingsList = ["memberName"]

    def __init__(self, parent = None, signalManager = None, name = "Merge data"):
        OWWidget.__init__(self, parent, signalManager, name)  #initialize base class

        # set channels
        self.inputs = [("Examples A", ExampleTable, self.onDataAInput), ("Examples B", ExampleTable, self.onDataBInput)]
        self.outputs = [("Merged Examples A+B", ExampleTable), ("Merged Examples B+A", ExampleTable)]

        # data
        self.dataA = None
        self.dataB = None
        self.varListA = []
        self.varListB = []
        self.varA = None
        self.varB = None
        self.lbAttrAItems = []
        self.lbAttrBItems = []


        # load settings
        self.loadSettings()
        
        # GUI
        w = QWidget(self)
        self.controlArea.layout().addWidget(w)
        grid = QGridLayout()
        grid.setMargin(0)
        w.setLayout(grid)

        # attribute A
        boxAttrA = OWGUI.widgetBox(self, 'Attribute A', orientation = "vertical", addToLayout=0)
        grid.addWidget(boxAttrA, 0,0)
        self.lbAttrA = OWGUI.listBox(boxAttrA, self, "lbAttrAItems", callback = self.lbAttrAChange)

        # attribute  B
        boxAttrB = OWGUI.widgetBox(self, 'Attribute B', orientation = "vertical", addToLayout=0)
        grid.addWidget(boxAttrB, 0,1)
        self.lbAttrB = OWGUI.listBox(boxAttrB, self, "lbAttrBItems", callback = self.lbAttrBChange)

        # info A
        boxDataA = OWGUI.widgetBox(self, 'Data A', orientation = "vertical", addToLayout=0)
        grid.addWidget(boxDataA, 1,0)
        self.lblDataAExamples = OWGUI.widgetLabel(boxDataA, "num examples")
        self.lblDataAAttributes = OWGUI.widgetLabel(boxDataA, "num attributes")

        # info B
        boxDataB = OWGUI.widgetBox(self, 'Data B', orientation = "vertical", addToLayout=0)
        grid.addWidget(boxDataB, 1,1)
        self.lblDataBExamples = OWGUI.widgetLabel(boxDataB, "num examples")
        self.lblDataBAttributes = OWGUI.widgetLabel(boxDataB, "num attributes")

        # general Info
        boxDataB = OWGUI.widgetBox(self, 'Info', orientation = "vertical", addToLayout=0)
        grid.addWidget(boxDataB, 2,0,1,2)
        self.generalInfo = OWGUI.widgetLabel(boxDataB, "The meta attributes will be present at the output table as regular attributes.")

        # icons
        self.icons = self.createAttributeIconDict()

        # resize
        self.resize(400,500)


    ############################################################################################################################################################
    ## Data input and output management
    ############################################################################################################################################################

    def onDataAInput(self, data):
        # set self.dataA, generate new domain if it is the same as of self.dataB.domain 
        self.dataA = data
        # update self.varListA and self.varA
        if self.dataA:
            self.varListA = [var.name for var in self.dataA.domain] + [var.name for var in self.dataA.domain.getmetas().values()] #self.dataA.domain.variables.native() + self.dataA.domain.getmetas().values()
        else:
            self.varListA = []
        if not self.varA in self.varListA:
            self.varA = None
        # update info
        self.updateInfoA()
        # update attribute A listbox
        self.lbAttrA.clear()
        for var in self.varListA:
            if var in [x.name for x in self.dataA.domain.getmetas().values()]:
                icon = OWGUI.createAttributePixmap("M", Qt.blue)
            else:
                icon = self.icons[self.dataA.domain[var].varType]
            self.lbAttrA.addItem(QListWidgetItem(icon, var))
        self.sendData()


    def onDataBInput(self, data):
        # set self.dataB, generate new domain if it is the same as of self.dataA.domain 
        self.dataB = data
        # update self.varListB and self.varB
        if self.dataB:
            self.varListB = [var.name for var in self.dataB.domain]  + [var.name for var in self.dataB.domain.getmetas().values()] #self.dataB.domain.variables.native() + self.dataB.domain.getmetas().values()
        else:
            self.varListB = []
        if not self.varB in self.varListB:
            self.varB = None
        # update info
        self.updateInfoB()
        # update attribute B listbox
        self.lbAttrB.clear()
        for var in self.varListB:
            print var
            print type(var)
            if var in [x.name for x in self.dataB.domain.getmetas().values()]:
                icon = OWGUI.createAttributePixmap("M", Qt.blue)
                print "META"
            else:
                icon = self.icons[self.dataB.domain[var].varType]
            print "icon:",icon,type(icon)
            print "var:",var,type(var)
            self.lbAttrB.addItem(QListWidgetItem(icon, var))
        self.sendData()


    def updateInfoA(self):
        """Updates data A info box.
        """
        if self.dataA:
            self.lblDataAExamples.setText("%s example%s" % self._sp(self.dataA))
            self.lblDataAAttributes.setText("%s attribute%s" % self._sp(self.varListA))
        else:
            self.lblDataAExamples.setText("No data on input A.")
            self.lblDataAAttributes.setText("")
        

    def updateInfoB(self):
        """Updates data B info box.
        """
        if self.dataB:
            self.lblDataBExamples.setText("%s example%s" % self._sp(self.dataB))
            self.lblDataBAttributes.setText("%s attribute%s" % self._sp(self.varListB))
        else:
            self.lblDataBExamples.setText("No data on input B.")
            self.lblDataBAttributes.setText("")


    def sendData(self):
        """Sends out data.
        """
        if self.varA and self.varB:
            etAB = dataUtilities.horizontalMerge(self.dataA,self.dataB,self.varA,self.varB)
            self.send("Merged Examples A+B", etAB)
            
            etBA = dataUtilities.horizontalMerge(self.dataB,self.dataA,self.varB,self.varA)
            self.send("Merged Examples B+A", etBA)
        else:
            self.send("Merged Examples A+B", None)
            self.send("Merged Examples B+A", None)


    ############################################################################################################################################################
    ## Event handlers
    ############################################################################################################################################################

    def lbAttrAChange(self):
        if self.dataA:
            if self.lbAttrA.selectedItems() != []:
                ind = self.lbAttrA.row(self.lbAttrA.selectedItems()[0])
                self.varA = self.varListA[ind]
            else:
                self.varA = None
        else:
            self.varA = None
        self.sendData()


    def lbAttrBChange(self):
        if self.dataB:
            if self.lbAttrB.selectedItems() != []:
                ind = self.lbAttrB.row(self.lbAttrB.selectedItems()[0])
                self.varB = self.varListB[ind]
            else:
                self.varB = None
        else:
            self.varB = None
        self.sendData()



    ############################################################################################################################################################
    ## Utility functions 
    ############################################################################################################################################################

    def _sp(self, l, capitalize=True):
        """Input: list; returns tupple (str(len(l)), "s"/"")
        """
        n = len(l)
        if n == 0:
            if capitalize:                    
                return "No", "s"
            else:
                return "no", "s"
        elif n == 1:
            return str(n), ''
        else:
            return str(n), 's'



if __name__=="__main__":
    import sys
    import OWDataTable, orngSignalManager
    signalManager = orngSignalManager.SignalManager(0)
    #data = dataUtilities.DataTable('dicty_800_genes_from_table07.tab')
##    data = dataUtilities.DataTable(r'..\..\doc\datasets\adult_sample.tab')
##    dataA = dataUtilities.DataTable(r'c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\sterolgene v.0 mouse probeRatios.tab')
##    dataA = dataUtilities.DataTable(r'c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\Copy of sterolgene v.0 mouse probeRatios.tab')
##    dataB = dataUtilities.DataTable(r'c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\sterolgene v.0 mouse probeRatios.tab')
    dataA = dataUtilities.DataTable(r'/home/pedro/DataSets/hMerge/D1.tab')
    dataB = dataUtilities.DataTable(r'/home/pedro/DataSets/hMerge/D2.tab')
    a=QApplication(sys.argv)
    ow=OWMergeData()
    a.setMainWidget(ow)
    ow.show()
    ow.onDataAInput(dataA)
    ow.onDataBInput(dataB)
    # data table
    dt = OWDataTable.OWDataTable(signalManager = signalManager)
    signalManager.addWidget(ow)
    signalManager.addWidget(dt)
    signalManager.setFreeze(1)
    signalManager.addLink(ow, dt, 'Merged Examples A+B', 'Examples', 1)
    signalManager.addLink(ow, dt, 'Merged Examples B+A', 'Examples', 1)
    signalManager.setFreeze(0)
    dt.show()
    a.exec_loop()

