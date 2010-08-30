"""
<name>Select Attributes</name>
<description>Manual selection of attributes.<br>Aside from a data set with the selected attributes, the widget also outputs an "Attribute List" which can be used with another "Select Attribute" widget to make the same variable selection</description>
<icon>icons/SelectAttributes.png</icon>
<priority>1100</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si) Changed by Pedro Almeida</contact>
"""
from AZutilities import dataUtilities
from OWTools import *
from OWWidget import *
import OWGUI,OWGUIEx

class OWDataDomain(OWWidget):
    contextHandlers = {"": DomainContextHandler("", [ContextField("chosenAttributes",
                                                                  DomainContextHandler.RequiredList,
                                                                  selected="selectedChosen", reservoir="inputAttributes"),
                                                     ContextField("classAttribute",
                                                                  DomainContextHandler.RequiredList,
                                                                  selected="selectedClass", reservoir="inputAttributes"),
                                                     ContextField("metaAttributes",
                                                                  DomainContextHandler.RequiredList,
                                                                  selected="selectedMeta", reservoir="inputAttributes")
                                                     ])}


    def __init__(self,parent = None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "Data Domain",wantMainArea = 0) #initialize base class

        self.inputs = [("Examples", ExampleTable, self.onDataInput), ("Attribute Subset", AttributeList, self.onAttributeList)]
        self.outputs = [("Examples", ExampleTable), ("Classified Examples", ExampleTable), ("Selected Attributes", AttributeList)]

        buttonWidth = 50
        applyButtonWidth = 101

        self.data = None
        self.receivedAttrList = []

        self.selectedInput = []
        self.inputAttributes = []
        self.selectedChosen = []
        self.chosenAttributes = []
        self.selectedClass = []
        self.classAttribute = []
        self.metaAttributes = []
        self.selectedMeta = []
        self.loadSettings()
        self.usedAttributes = {}
        self.userWarned = False
       
        import sip
        sip.delete(self.controlArea.layout()) 
        gl=QGridLayout()
        self.controlArea.setLayout(gl)
        gl.setMargin(0)

        boxAvail = OWGUI.widgetBox(self,'Available Attributes')
        gl.addWidget(boxAvail, 0,0,3,1)

        self.filterInputAttrs = OWGUIEx.lineEditFilter(boxAvail, self, None, useRE = 1, emptyText = "filter attributes...", callback = self.setInputAttributes, caseSensitive = 0)
        self.inputAttributesList = OWGUI.listBox(boxAvail, self, "selectedInput", "inputAttributes", callback = self.onSelectionChange, selectionMode = QListWidget.ExtendedSelection)
        self.filterInputAttrs.listbox = self.inputAttributesList

        vbAttr = OWGUI.widgetBox(self, addToLayout = 0)
        gl.addWidget(vbAttr, 0,1)
        self.attributesButtonUp = OWGUI.button(vbAttr, self, "Up", self.onAttributesButtonUpClick)
        self.attributesButtonUp.setMaximumWidth(buttonWidth)
        self.attributesButton = OWGUI.button(vbAttr, self, ">",self.onAttributesButtonClicked)        
        self.attributesButton.setMaximumWidth(buttonWidth)
        self.attributesButtonDown = OWGUI.button(vbAttr, self, "Down", self.onAttributesButtonDownClick)
        self.attributesButtonDown.setMaximumWidth(buttonWidth)
        
        boxAttr = OWGUI.widgetBox(self,'Attributes', addToLayout = 0)
        gl.addWidget(boxAttr, 0,2)
        self.attributesList = OWGUI.listBox(boxAttr, self, "selectedChosen", "chosenAttributes", callback = self.onSelectionChange, selectionMode = QListWidget.ExtendedSelection)

        self.classButton = OWGUI.button(self, self, ">", self.onClassButtonClicked, addToLayout = 0)
        self.classButton.setMaximumWidth(buttonWidth)
        gl.addWidget(self.classButton, 1,1)
        boxClass = OWGUI.widgetBox(self,'Class', addToLayout = 0)
        boxClass.setFixedHeight(55)
        gl.addWidget(boxClass, 1,2)
        self.classList = OWGUI.listBox(boxClass, self, "selectedClass", "classAttribute", callback = self.onSelectionChange, selectionMode = QListWidget.ExtendedSelection)
        
        vbMeta = OWGUI.widgetBox(self, addToLayout = 0)
        gl.addWidget(vbMeta, 2,1)
        self.metaButtonUp = OWGUI.button(vbMeta, self, "Up", self.onMetaButtonUpClick)
        self.metaButtonUp.setMaximumWidth(buttonWidth)
        self.metaButton = OWGUI.button(vbMeta, self, ">",self.onMetaButtonClicked)
        self.metaButton.setMaximumWidth(buttonWidth)
        self.metaButtonDown = OWGUI.button(vbMeta, self, "Down", self.onMetaButtonDownClick)
        self.metaButtonDown.setMaximumWidth(buttonWidth)
        boxMeta = OWGUI.widgetBox(self,'Meta Attributes', addToLayout = 0)
        gl.addWidget(boxMeta, 2,2)
        self.metaList = OWGUI.listBox(boxMeta, self, "selectedMeta", "metaAttributes", callback = self.onSelectionChange, selectionMode = QListWidget.ExtendedSelection)
        
        boxApply = OWGUI.widgetBox(self, addToLayout = 0, orientation = "horizontal", addSpace = 1) #QHBox(ca)
        gl.addWidget(boxApply, 3,0,1,3)
        self.applyButton = OWGUI.button(boxApply, self, "Apply", callback = self.setOutput)
        self.applyButton.setEnabled(False)
        self.applyButton.setMaximumWidth(applyButtonWidth)
        self.resetButton = OWGUI.button(boxApply, self, "Reset", callback = self.reset)
        self.resetButton.setMaximumWidth(applyButtonWidth)

        infoBox = OWGUI.widgetBox(self,"Info" ,addToLayout = 0, orientation = "horizontal", addSpace = 1)
        gl.addWidget(infoBox, 5,0,1,3)
        OWGUI.widgetLabel(infoBox, 'Aside from a data set with the selected attributes, \nthe widget also outputs an "Attribute List" which can be used with\n another "Select Attribute" widget to make the same variable selection')

        gl.setRowStretch(0, 4)
        gl.setRowStretch(1, 0)
        gl.setRowStretch(2, 2)
        

        self.icons = self.createAttributeIconDict()

        self.inChange = False
        self.resize(400,480)       

    def onAttributeList(self, attrList):
        ##scPA
        #The received attrList must be a list being the first element a list of
        #selected attributes, the second a list with one elements, the class attribute (may be an empty list) 
        #and the third element a list of metaAttributes Ex: [[a,b,c,d], [g], [e,f]]
        if attrList and type(attrList)==list:
            if type(attrList[0]) != list:
                self.receivedAttrList =[[attrList]]
            else:
                self.receivedAttrList = attrList
            self.updateSelection()
        else:
            self.receivedAttrList = []
        ##ecPA

    def onSelectionChange(self):
        if not self.inChange:
            self.inChange = True
            for lb, co in [(self.inputAttributesList, "selectedInput"),
                       (self.attributesList, "selectedChosen"),
                       (self.classList, "selectedClass"),
                       (self.metaList, "selectedMeta")]:
                if not lb.hasFocus():
                    setattr(self, co, [])
            self.inChange = False

        self.updateInterfaceState()            

    def onDataInput(self, data): ##scPA Added parameter reloadData##ecPA
        self.userWarned = False
        self.data = data
        self.updateSelection()

    def updateSelection(self):
        ##scPA - Commented because if the data changes, it will assume as the same dada
        #if self.data and data and self.data.checksum() == data.checksum() and not reloadData:
        #    return   # we received the same dataset again

        self.closeContext()
        #local variables for the list boxes. the actual list boxes variables must be only 
        #changed once. If not done this way, the GUI intergace will not be properly updated!
        locatMetaList = []
        localAttrList = []
        localClassList = []
        localInputList = []

        if self.data:
            domain = self.data.domain
            metaIds = domain.getmetas().keys()
            metaIds.sort()
            self.allAttributes = [(attr.name, attr.varType) for attr in domain] + [(domain[i].name, domain[i].varType) for i in metaIds]
            #Make attributes selection based only on name and case insensitive
            allAttributesNames = [x[0] for x in self.allAttributes]

            ##scPA
            if self.receivedAttrList:
                if len(self.receivedAttrList)>=1 and self.receivedAttrList[0]: 
                    chosenAttr = [x.name for x in self.receivedAttrList[0]]
                    for attr in chosenAttr:
                        if attr in allAttributesNames:
                            localAttrList.append(self.allAttributes[allAttributesNames.index(attr)])

                if len(self.receivedAttrList) >= 2 and self.receivedAttrList[1]:
                    classAttr = self.receivedAttrList[1][0].name
                    if classAttr in allAttributesNames:
                        localClassList = [self.allAttributes[allAttributesNames.index(classAttr)]]

                if len(self.receivedAttrList)>=3 and self.receivedAttrList[2]:
                    metaAttr=[x.name for x in self.receivedAttrList[2]] 
                    for attr in metaAttr:
                        if attr in allAttributesNames:
                            locatMetaList.append(self.allAttributes[allAttributesNames.index(attr)])
            else:
                localAttrList = [(a.name, a.varType) for a in domain.attributes]
                if domain.classVar:
                    localClassList = [(domain.classVar.name, domain.classVar.varType)]
                locatMetaList = [(a.name, a.varType) for a in domain.getmetas().values()]


        if not self.receivedAttrList:
            self.openContext("",self.data)
        else:
            self.metaAttributes = locatMetaList
            self.chosenAttributes = localAttrList
            self.classAttribute = localClassList
            self.inputAttributes = localInputList


        self.usedAttributes = dict.fromkeys(self.chosenAttributes + self.classAttribute + self.metaAttributes, 1)
        self.setInputAttributes()
        self.setOutput()
        self.updateInterfaceState()
        #self.metaList.triggerUpdate(True)
        #self.attributesList.triggerUpdate(True)
        #self.inputAttributesList.triggerUpdate(True)
        #self.classList.triggerUpdate(True)

    
    def setOutput(self):
        if self.data:
            self.applyButton.setEnabled(False)
            attributes = []
            attributes.append([self.data.domain[x[0]] for x in self.chosenAttributes])
            if self.classAttribute:
                attributes.append([self.data.domain[self.classAttribute[0][0]] ])
                classVar = attributes[1][0] 
            else:
                attributes.append([])
                classVar = None
            attributes.append([self.data.domain[x[0]] for x in self.metaAttributes])

            self.send("Selected Attributes", attributes)

            domain = orange.Domain(attributes[0], classVar)
            for meta in self.metaAttributes:
                domain.addmeta(orange.newmetaid(), self.data.domain[meta[0]])

            newdata = dataUtilities.DataTable(domain, self.data)
            newdata.name = self.data.name
            if not self.userWarned and len(newdata.domain.getmetas()) != 0:            
                QMessageBox.warning( None, "Select - Meta Attributes", "There are meta-attributes present in the dataset.\nThe presence of meta-Attributes in datasets used with Learners/Classifiers\nrequires the use of considerably more ram memory!" , QMessageBox.Ok) 
                self.userWarned = True
            if len(newdata.domain) != 0:
                self.send("Examples", newdata)
            else:
                self.send("Examples", None)
                self.send("Classified Examples", None)
                return

            if classVar:
                self.send("Classified Examples", newdata)
            else:
                self.send("Classified Examples", None)
        else:
            self.send("Examples", None)
            self.send("Classified Examples", None)

        
    def reset(self):
        data = self.data
        self.data = None
        self.onDataInput(data)

        
    def disableButtons(self, *arg):
        for b in arg:
            b.setEnabled(False)

    def setButton(self, button, dir):
        button.setText(dir)
        button.setEnabled(True)

        
    def updateInterfaceState(self):
        if self.selectedInput:
            self.setButton(self.attributesButton, ">")
            self.setButton(self.metaButton, ">")
            self.disableButtons(self.attributesButtonUp, self.attributesButtonDown, self.metaButtonUp, self.metaButtonDown)

            if len(self.selectedInput) == 1 and self.inputAttributes[self.selectedInput[0]][1] in [orange.VarTypes.Discrete, orange.VarTypes.Continuous]:
                self.setButton(self.classButton, ">")
            else:
                self.classButton.setEnabled(False)
            
        elif self.selectedChosen:
            self.setButton(self.attributesButton, "<")
            self.disableButtons(self.classButton, self.metaButton, self.metaButtonUp, self.metaButtonDown)

            mini, maxi = min(self.selectedChosen), max(self.selectedChosen)
            cons = maxi - mini == len(self.selectedChosen) - 1
            self.attributesButtonUp.setEnabled(cons and mini)
            self.attributesButtonDown.setEnabled(cons and maxi < len(self.chosenAttributes)-1)

        elif self.selectedClass:
            self.setButton(self.classButton, "<")
            self.disableButtons(self.attributesButtonUp, self.attributesButtonDown, self.metaButtonUp, self.metaButtonDown,
                                self.attributesButton, self.metaButton)

        elif self.selectedMeta:
            self.setButton(self.metaButton, "<")
            self.disableButtons(self.attributesButton, self.classButton, self.attributesButtonDown, self.attributesButtonUp)

            mini, maxi, leni = min(self.selectedMeta), max(self.selectedMeta), len(self.selectedMeta)
            cons = maxi - mini == leni - 1
            self.metaButtonUp.setEnabled(cons and mini)
            self.metaButtonDown.setEnabled(cons and maxi < len(self.metaAttributes)-1)

        else:
            self.disableButtons(self.attributesButtonUp, self.attributesButtonDown, self.metaButtonUp, self.metaButtonDown,
                                self.attributesButton, self.metaButton, self.classButton)


    def splitSelection(self, alist, selected):
        selected.sort()

        i, sele = 0, selected[0]
        selList, restList = [], []
        for j, attr in enumerate(alist):
            if j == sele:
                selList.append(attr)
                i += 1
                sele = i<len(selected) and selected[i] or None
            else:
                restList.append(attr)
        return selList, restList


    def setInputAttributes(self):
        self.selectedInput = []
        if self.data:
            self.inputAttributes = filter(lambda x:not self.usedAttributes.has_key(x), self.allAttributes)
        else:
            self.inputAttributes = []
        self.filterInputAttrs.setAllListItems()
        self.filterInputAttrs.updateListBoxItems(callCallback = 0)

    def removeFromUsed(self, attributes):
        for attr in attributes:
            del self.usedAttributes[attr]
        self.setInputAttributes()

    def addToUsed(self, attributes):
        self.usedAttributes.update(dict.fromkeys(attributes))
        self.setInputAttributes()

       
    def onAttributesButtonClicked(self):
        if self.selectedInput:
            selList, restList = self.splitSelection(self.inputAttributes, self.selectedInput)
            self.chosenAttributes = self.chosenAttributes + selList
            self.addToUsed(selList)
        else:
            selList, restList = self.splitSelection(self.chosenAttributes, self.selectedChosen)
            self.chosenAttributes = restList
            self.removeFromUsed(selList)

        self.updateInterfaceState()
        self.applyButton.setEnabled(True)


    def onClassButtonClicked(self):
        if self.selectedInput:
            selected = self.inputAttributes[self.selectedInput[0]]
            if self.classAttribute:
                self.removeFromUsed(self.classAttribute)
            self.addToUsed([selected])
            self.classAttribute = [selected]
        else:
            self.removeFromUsed(self.classAttribute)
            self.selectedClass = []
            self.classAttribute = []

        self.updateInterfaceState()
        self.applyButton.setEnabled(True)        


    def onMetaButtonClicked(self):
        if not self.userWarned:
            QMessageBox.warning( self, "Select - Meta Attributes", "The presence of meta-Attributes in datasets used with Learners/Classifiers\nrequires the use of considerably more ram memory!" , QMessageBox.Ok)
            self.userWarned = True
        if self.selectedInput:
            selList, restList = self.splitSelection(self.inputAttributes, self.selectedInput)
            self.metaAttributes = self.metaAttributes + selList
            self.addToUsed(selList)
        else:
            selList, restList = self.splitSelection(self.metaAttributes, self.selectedMeta)
            self.metaAttributes = restList
            self.removeFromUsed(selList)

        self.updateInterfaceState()
        self.applyButton.setEnabled(True)


    def moveSelection(self, labels, selection, dir):
        labs = getattr(self, labels)
        sel = getattr(self, selection)
        mini, maxi = min(sel), max(sel)+1
        if dir == -1:
            setattr(self, labels, labs[:mini-1] + labs[mini:maxi] + [labs[mini-1]] + labs[maxi:])
        else:
            setattr(self, labels, labs[:mini] + [labs[maxi]] + labs[mini:maxi] + labs[maxi+1:])
        setattr(self, selection, map(lambda x:x+dir, sel))
        self.updateInterfaceState()
        self.applyButton.setEnabled(True)

    def onMetaButtonUpClick(self):
        self.moveSelection("metaAttributes", "selectedMeta", -1)

    def onMetaButtonDownClick(self):
        self.moveSelection("metaAttributes", "selectedMeta", 1)
        
    def onAttributesButtonUpClick(self):
        self.moveSelection("chosenAttributes", "selectedChosen", -1)
        
    def onAttributesButtonDownClick(self):
        self.moveSelection("chosenAttributes", "selectedChosen", 1)
        

if __name__=="__main__":
    import sys
    data = dataUtilities.DataTable(r'..\..\doc\datasets\iris.tab')
    # add meta attribute
    data.domain.addmeta(orange.newmetaid(), orange.StringVariable("name"))
    for ex in data:
        ex["name"] = str(ex.getclass())

    a=QApplication(sys.argv)
    ow=OWDataDomain()
    a.setMainWidget(ow)
    ow.show()
    ow.onDataInput(data)
    a.exec_loop()
    ow.saveSettings()
