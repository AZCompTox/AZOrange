"""
<name>Concatenate</name>
<description>Vertical merge of examples. </description>
<icon>icons/Concatenate.png</icon>
<priority>1111</priority>
<contact>Janez Demsar (janez.demsar(@at@)fri.uni-lj.si)</contact>
"""

from OWWidget import *
import OWGUI
from AZutilities import dataUtilities

class OWConcatenate(OWWidget):
    settingsList = ["mergeAttributes"]
    minReqVer = 2
    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "FeatureConstructor")
        #To enable the use of primary table, uncomment the next line  and comment the following one
        #self.inputs = [("Primary Table", dataUtilities.DataTable, self.setData), ("Additional Tables", dataUtilities.DataTable, self.setMoreData, Multiple)]
        self.inputs = [("Tables", orange.ExampleTable, self.setMoreData, Multiple)]

        self.outputs = [("Examples", ExampleTable)]

        self.mergeAttributes = 0
        self.primary = None
        self.additional = {}
        
        bg = self.bgMerge = OWGUI.radioButtonsInBox(self.controlArea, self, "mergeAttributes", [], "Domains merging", callback = self.apply)
        OWGUI.widgetLabel(bg, "") # "When there is no primary table, the domain should be")
        OWGUI.appendRadioButton(bg, self, "mergeAttributes", "Union of attributes appearing in all tables")
        OWGUI.appendRadioButton(bg, self, "mergeAttributes", "Intersection of attributes in all tables")
        #OWGUI.widgetLabel(bg, "The resulting table will have class only if there is no conflict betwen input classes.")
        
        infoBox = OWGUI.widgetBox(self,"Concatenate status")
        self.infoStatus = OWGUI.widgetLabel(infoBox, '')

        self.viewBT = OWGUI.button(self.controlArea, self, "&View converted incompatible attributes", callback = self.viewIncompatible, disabled=1)

        self.adjustSize()

    def viewIncompatible(self):
        QMessageBox.warning( None, "Concatenate widget","Some attributes were incompatible. They have the same name but they are represented by different types.\n"+\
                        "\nIncompatible attributes:\n"+\
                        self.attrsMsg)
 

    def setData(self, data):
        self.primary = data
        self.bgMerge.setEnabled(not data)
        self.apply()
        

    def setMoreData(self, data, id):
        if not data:
            if id in self.additional:
                del self.additional[id]
        else:
            self.additional[id] = data
        self.apply()
        
 
    def apply(self):
        incompatibleAttr = []
        if not hasattr(dataUtilities,"version") or (dataUtilities.version < self.minReqVer):
           QMessageBox.critical( None, "Concatenate widget","This widget requires module dataUtilities version " + str(self.minReqVer) + "\nAborting Concatenation now!")
           self.error(0,"This widget requires module dataUtilities version " + str(self.minReqVer) )
           self.infoStatus.setText('This widget requires module dataUtilities version ' + str(self.minReqVer))
           self.send("Examples", None)
           return
        self.infoStatus.setText('')
        self.warning(0)
        self.error(0)
        if self.primary:
            if not self.additional:
                newTable = self.primary
            else:
                allDatasets = [self.primary] + [table for table in self.additional.values()]
                newTable,status = dataUtilities.concatenate(allDatasets, useFirstAsLeader = True)
                for attr in status:
                    if attr not in incompatibleAttr: incompatibleAttr.append(attr)
        else:
            if not self.additional:
                newTable = None     
            else:
                newTable, status = dataUtilities.concatenate([table for table in self.additional.values()], useFirstAsLeader = False, mergeDomains = (self.mergeAttributes==0))
                for attr in status:
                    if attr not in incompatibleAttr: incompatibleAttr.append(attr)
        if not newTable:
            status = "No data"
        else:
            # Update Status
            status = "New concatenated domain attributes: " + str(len(newTable.domain.attributes)) + "\n"
            status = status + "New concatenated Data examples: " + str(len(newTable)) + "\n"
            status = status + "New concatenated class: " + str(newTable.domain.classVar) + "\n"
            if  len(incompatibleAttr) > 0:
                self.warning(0,"Some imcompatible attributes were converted")
                self.attrsMsg=""
                for attr in incompatibleAttr:
                    self.attrsMsg += "  " + attr + "      (used " + str(newTable.domain[attr]).replace("'"+attr+"'","") + ")" + "\n"
                status = status + "Incompatible Attributes converted: " + str(len(incompatibleAttr)) + "    (Check original datasets)\n"
                self.viewBT.setDisabled(False)
            else:
                self.viewBT.setDisabled(True)
                self.attrsMsg=""
        self.infoStatus.setText(status)

        self.send("Examples", newTable)
