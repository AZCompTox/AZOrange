"""
<name>Rename Class</name>
<description>Renames the class variable of a dataset</description>
<icon>icons/RenameClass.png</icon>
<priority>1111</priority>
<contact>Pedro Rafael Almeida (engpedrorafael(@at@)gmail.com)</contact>
"""

from OWWidget import *
import OWGUI
from AZutilities import dataUtilities
import types

class OWRenameClass(OWWidget):
    minReqVer = 4
    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "FeatureConstructor",wantMainArea = 0)
        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples", ExampleTable)]

        self.newName = "CLASS"
        self.data = None
        
        box = OWGUI.widgetBox(self.controlArea, "Rename class")
        newNameBox = OWGUI.lineEdit(box, self, "newName", "New Class name: ", labelWidth=100, orientation="horizontal", tooltip="New name for the Class attribute")
 
        self.BT = OWGUI.button(self.controlArea, self, "&Apply", callback = self.apply, disabled=0)
        infoBox = OWGUI.widgetBox(self,"Info")
        self.info = OWGUI.widgetLabel(infoBox, 'Renames the class attribute of the data.')


        self.adjustSize()

    def setData(self, data):
        self.data = data
        self.error(0)
        self.warning(0)
        if not data:
            self.send("Examples", None)
        else:
            self.apply()
        

    def apply(self):
        self.error(0)
        self.warning(0)
        if not self.data:
            return

        if not self.newName:
            self.warning(0,'No new name specified for the class. The class name was not changed')
            self.send("Examples" , self.data)
            return

        if not self.data.domain.classVar:
            self.warning(0,'The data connected is classless')
            self.send("Examples",self.data)
            return
        self.info.setText("Actual class name: " + self.data.domain.classVar.name)

        newData = dataUtilities.changeClassName(self.data , self.newName)
        if not newData or not newData.domain.classVar:
            self.error(0,'It was not possible to rename the class variable. Invalid data or name.')
            self.send("Examples",None)
            return

        if newData.domain.classVar.name != self.newName:
            self.warning(0,"The name specified already existed. The class was renamed to " + newData.domain.classVar.name)
            

        self.send("Examples", newData)
