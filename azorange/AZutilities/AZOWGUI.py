from OWWidget import *
import OWGUI
import OWGUIEx

             
class widgetSelector(QWidget):
    """
        Class to provide a left-right select widget with an apply and reset button
        usage:
            obj = widgetSelector(self.controlArea, selectText = "Descriptors", applyText = "Retrieve", callbackOnApply=self.OnRetrieve, callbackOnReset = None)
    """
    def __init__(self,parent=None, selectText = "", applyText = "Apply", callbackOnApply=None, callbackOnReset = None):
        QWidget.__init__(self)
        if parent:
            parent.layout().addWidget(self)
        gl = QGridLayout()
        gl.setMargin(0)
        self.setLayout(gl)
        self.callbackOnApply = callbackOnApply
        self.callbackOnReset =callbackOnReset
        #Local Variables
        buttonWidth = 50
        applyButtonWidth = 101

        self.inputItems = []            # A backup of all input Items
        self.availableItems = []        # A managed list of the available items that can still be selected (inputItems that are not yet selectedItems)
        self.availableItems2 = ["a","b","c"]        # A managed list of the available items that can still be selected (inputItems that are not yet selectedItems)
        self.selectedItems = []         # A managed list of the already selected items

        
        # Left pane Box.
        boxAvail = OWGUI.widgetBox(self,'Available '+selectText,addToLayout = 0)
        gl.addWidget(boxAvail, 0,0,3,1)

        self.filterInputItems = OWGUIEx.lineEditFilter(boxAvail, self, None, useRE = 1, emptyText = "filter "+selectText+"...", callback = self.__setFilteredInput, caseSensitive = 0)

        self.inputItemsList = OWGUI.listBox(boxAvail, self, None, None, selectionMode = QListWidget.ExtendedSelection)

        self.filterInputItems.listbox = self.inputItemsList

        vbItems = OWGUI.widgetBox(self)
        gl.addWidget(vbItems, 0,1)
        self.ButtonAdd = OWGUI.button(vbItems, self, ">",self.__onButtonAddClicked)
        self.ButtonAdd.setMaximumWidth(buttonWidth)
        self.ButtonRemove = OWGUI.button(vbItems, self, "<",self.__onButtonRemoveClicked)
        self.ButtonRemove.setMaximumWidth(buttonWidth)

        # Right pane selected descriptors.
        box = OWGUI.widgetBox(self,'Selected '+selectText)
        gl.addWidget(box, 0,2)
        self.selectedItemsList = OWGUI.listBox(box,self,None, None,selectionMode = QListWidget.ExtendedSelection)

        #Apply Reset buttons
        boxApply = OWGUI.widgetBox(self, '', orientation="horizontal",addToLayout = 0)
        gl.addWidget(boxApply, 3,0,1,3)
        if self.callbackOnApply != None:
            self.applyButton = OWGUI.button(boxApply, self, applyText, callback = self.__apply)
            self.applyButton.setMaximumWidth(applyButtonWidth)
        self.resetButton = OWGUI.button(boxApply, self, "Reset", callback = self.__reset)
        self.resetButton.setMaximumWidth(applyButtonWidth)

    def getSelectedItems(self):
        return [str(x) for x in self.selectedItems]

    def __apply(self):
        if self.callbackOnApply != None:
            self.callbackOnApply([str(x) for x in self.selectedItems])
        else:
            msgItems = ""
            for x in self.selectedItems:
                msgItems += str(x)+"\n"
            QMessageBox.information("Missing Apply CallBack","There is not CallBack specifyed.\nThe Selected Items are:\n"+msgItems, QMessageBox.Ok)        

    def __updateLists(self):
        tmp = {}.fromkeys(self.availableItems)
        self.availableItems = [str(x) for x in tmp.keys()]

        tmp = {}.fromkeys(self.selectedItems)
        self.selectedItems = [str(x) for x in tmp.keys()]

        #Fill the ListBoxes
        self.inputItemsList.clear()
        self.inputItemsList.addItems(self.availableItems)
        self.selectedItemsList.clear()
        self.selectedItemsList.addItems(self.selectedItems)

        #Sort items
        self.inputItemsList.sortItems()
        self.selectedItemsList.sortItems()


    def __setFilteredInput(self):
        self.availableItems = [x for x in self.inputItems if x not in self.selectedItems]
        self.__updateLists()
        self.filterInputItems.setAllListItems()
        self.filterInputItems.updateListBoxItems(callCallback = 0)

    def __reset(self):
        self.availableItems = []
        self.selectedItems = []
        for item in self.inputItems:
            self.availableItems.append(item)
        self.__updateLists()
        if self.callbackOnReset != None:
            self.callbackOnReset()

    def __onButtonAddClicked(self):
        for item in self.inputItemsList.selectedItems():
            self.selectedItems.append(item.text())
            self.availableItems.remove(item.text())
        self.__updateLists()

    def __onButtonRemoveClicked(self):
        for item in self.selectedItemsList.selectedItems():
            self.availableItems.append(item.text())
            self.selectedItems.remove(item.text())
        self.__updateLists()

    def setInputItems(self, items):
        self.inputItems = [str(x) for x in items]
        self.__reset()


if __name__ == "__main__":
    a = QApplication(sys.argv)
    owdm = OWOptimizeMe()
    a.setMainWidget(owdm)
    owdm.show()
    a.exec_loop()
