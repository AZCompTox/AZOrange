"""
<name>Add Class</name>
<description>Adds a new class to a dataset</description>
<icon>icons/AddClass.png</icon>
<priority>1111</priority>
<contact>Pedro Rafael Almeida (engpedrorafael(@at@)gmail.com)</contact>
"""

from OWWidget import *
import OWGUI
from AZutilities import dataUtilities
import types

class OWAddClass(OWWidget):
    minReqVer = 4
    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "FeatureConstructor",wantMainArea = 0)
        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples", ExampleTable)]

        self.newName = "NEW_CLASS"
        self.value = "?"
        self.values = "[]"
        self.data = None
        
        box = OWGUI.widgetBox(self.controlArea,"Add discrete class")
        newNameBox = OWGUI.lineEdit(box, self, "newName", "Class name: ", labelWidth=70, orientation="horizontal", tooltip="Name of the new Class attribute to be added")
        valueBox = OWGUI.lineEdit(box, self, "value", "Value: ", labelWidth=70, orientation="horizontal", tooltip="Class value to be assigned to each example")
        #valuesBox = OWGUI.lineEdit(box, self, "values", "Class values: ", labelWidth=70, orientation="horizontal", tooltip='Optional, list of possible values of new discrete class attribute ex: ["POS", "NEG"]')
 
        self.BT = OWGUI.button(self.controlArea, self, "&Apply", callback = self.apply, disabled=0)
        infoBox = OWGUI.widgetBox(self,"Info")
        OWGUI.widgetLabel(infoBox, 'Adds a new discrete response variable\nto the data and gives all examples the same value.')


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
        try:
            values = eval(self.values)
            if type(values) != types.ListType:
                raise "Wrong values"
            else:
                for v in values:
                    if type(v) not in  types.StringTypes:
                        raise "Wrong values"
        except:
            values = []
            self.warning(0,'Invalid values list. Usage: ["A","B","C"]')
        #self.newName = self.newName.upper()
        newTable = dataUtilities.addDiscreteClass(self.data, self.newName,self.value,values)
        if not newTable:
            self.error(0,'It was not possible to add the attribute the way you specified. Please check again the options!')
        elif newTable.domain.classVar.name != self.newName:
            self.warning(0,'There was already an attribute named '+ self.newName+' the new class was created with name '+newTable.domain.classVar.name)
        self.send("Examples", newTable)
