"""
<name>File</name>
<description>Reads data from a file.</description>
<icon>icons/File.png</icon>
<contact>Janez Demsar (janez.demsar(@at@)fri.uni-lj.si) Mofifyed by Pedro Almeida</contact>
<priority>10</priority>
"""

# Don't move this - the line number of the call is important
def call(f,*args,**keyargs):
    return f(*args, **keyargs)

from OWWidget import *
import OWGUI, string, os.path, user, sys, warnings
from AZutilities import dataUtilities
import os

warnings.filterwarnings("error", ".*" , orange.KernelWarning, "OWFile", 11)


class FileNameContextHandler(ContextHandler):
    def match(self, context, imperfect, filename):
        return context.filename == filename and 2
        

class OWFile(OWWidget):
    ##scPA  Force the self.createNewOn to be by default orange.Variable.MakeStatus.OK despite what user used last time
    settingsList=["recentFiles"]#, "createNewOn"]
    contextHandlers = {"": FileNameContextHandler()}

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "File")

        self.inputs = []
        self.outputs = [("Examples", ExampleTable), ("Attribute Definitions", orange.Domain)]

        #set default settings
        self.recentFiles=["(none)"]
        self.symbolDC = "?"
        self.symbolDK = "~"
        self.createNewOn = 3 # Always create a new attribute #orange.Variable.MakeStatus.NoRecognizedValues
        self.domain = None
        self.loadedFile = ""
        #get settings from the ini file, if they exist
        self.loadSettings()

        box = OWGUI.widgetBox(self.controlArea, "Data File", addSpace = True, orientation=0)
        self.filecombo=QComboBox(box)
        self.filecombo.setMinimumWidth(150)
        box.layout().addWidget(self.filecombo)
        button = OWGUI.button(box, self, '...', callback = self.browseFile, disabled=0)
        button.setIcon(self.style().standardIcon(QStyle.SP_DirOpenIcon))
        button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed)
        self.reloadBtn = OWGUI.button(box, self, "Load", callback = self.reload)
        self.reloadBtn.setIcon(self.style().standardIcon(QStyle.SP_BrowserReload))
        self.reloadBtn.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        #button.setMaximumWidth(25)i
        #self.reloadBtn.setMaximumWidth(80)

        #OWGUI.rubber(self.controlArea)

        box = OWGUI.widgetBox(self.controlArea, "Missing Value Symbols", addSpace = True, orientation=1)#

        hbox = OWGUI.indentedBox(box)
        #box = OWGUI.widgetBox(self.controlArea, "Data File", addSpace = True, orientation=0)
        #OWGUI.widgetLabel(box, "Symbols for missing values in tab-delimited files (besides default ones)")
        #hbox = OWGUI.indentedBox(box, addSpace=True)
        le = OWGUI.lineEdit(hbox, self, "symbolDC", "Don't care:  ", labelWidth=100, orientation="horizontal", tooltip="Default values: empty fields (space), '?' or 'NA'")
        OWGUI.separator(hbox, 16, 0)
        le = OWGUI.lineEdit(hbox, self, "symbolDK", "Don't know:  ", labelWidth=100, orientation="horizontal", tooltip="Default values: '~' or '*'")
        OWGUI.radioButtonsInBox(self.controlArea, self, "createNewOn", box="Advanced", addSpace=True,
                       label = "Create a new attribute when existing attribute(s) ...",
                       btnLabels = ["Have mismatching order of values",
                                    "Have no common values with the new (recommended)",
                                    "Miss some values of the new attribute",
                                    "... Always create a new attribute"
                               ])
        OWGUI.rubber(self.controlArea)

        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.infoa = OWGUI.widgetLabel(box, 'No data loaded.')
        self.infob = OWGUI.widgetLabel(box, '')
        self.warnings = OWGUI.widgetLabel(box, '')

        # remove missing data set names
        self.recentFiles=filter(os.path.exists, self.recentFiles)
        self.setFileList()
        #Uncomment next to auto load last file oppenned
        #if len(self.recentFiles) > 0 and os.path.exists(self.recentFiles[0]):
        #   self.openFile(self.recentFiles[0], 0, self.symbolDK, self.symbolDC)
        self.connect(self.filecombo, SIGNAL('activated(int)'), self.selectFile)

        self.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding))
        #self.adjustSize()

    # set the file combo box
    def setFileList(self):
        self.filecombo.clear()
        if not self.recentFiles:
            self.filecombo.addItem("(none)")
        for file in self.recentFiles:
            if file == "(none)":
                self.filecombo.addItem("(none)")
            else:
                self.filecombo.addItem(os.path.split(file)[1])
        self.filecombo.addItem("Browse documentation data sets...")
        self.filecombo.updateGeometry()

##scPA
    #def remove(self, dw = 1, dsw = 1):
    #    self.linksOut.clear()
    #    self.signalManager.removeWidget(self)
    #    #del self.dataDomain

##ecPA

    def reload(self):
        idx = self.filecombo.currentIndex()
        file = self.recentFiles[idx]
        return self.openFile(file, 1, self.symbolDK, self.symbolDC)


    def settingsFromWidgetCallback(self, handler, context):
        context.filename = self.loadedFile
        context.symbolDC, context.symbolDK = self.symbolDC, self.symbolDK

    def settingsToWidgetCallback(self, handler, context):
        self.symbolDC, self.symbolDK = context.symbolDC, context.symbolDK

    # user selected a file from the combo box
    def selectFile(self, n):
        if n < len(self.recentFiles) :
            name = self.recentFiles[n]
            self.recentFiles.remove(name)
            self.recentFiles.insert(0, name)
        elif n:
            self.browseFile(1)

        if len(self.recentFiles) > 0:
            self.setFileList()
            self.openFile(self.recentFiles[0], 0, self.symbolDK, self.symbolDC)

    # user pressed the "..." button to manually select a file to load
    def browseFile(self, inDemos=0):
        "Display a FileDialog and select a file"
        if inDemos:
            try:
                import win32api, win32con
                t = win32api.RegOpenKey(win32con.HKEY_LOCAL_MACHINE, "SOFTWARE\\Python\\PythonCore\\%i.%i\\PythonPath\\Orange" % sys.version_info[:2], 0, win32con.KEY_READ)
                t = win32api.RegQueryValueEx(t, "")[0]
                startfile = t[:t.find("orange")] + "orange\\doc\\datasets"
            except:
                startfile = ""

            if not startfile or not os.path.exists(startfile):
                d = OWGUI.__file__
                if d[-8:] == "OWGUI.py":
                    startfile = d[:-22] + "doc/datasets"
                elif d[-9:] == "OWGUI.pyc":
                    startfile = d[:-23] + "doc/datasets"

            if not startfile or not os.path.exists(startfile):
                d = os.getcwd()
                if d[-12:] == "OrangeCanvas":
                    startfile = d[:-12]+"doc/datasets"
                else:
                    if d[-1] not in ["/", "\\"]:
                        d+= "/"
                    startfile = d+"doc/datasets"

            if not os.path.exists(startfile):
                QMessageBox.information( None, "File", "Cannot find the directory with example data sets", QMessageBox.Ok + QMessageBox.Default)
                return
        else:
            if len(self.recentFiles) == 0 or self.recentFiles[0] == "(none)":
                if sys.platform == "darwin":
                    startfile = user.home
                else:
                    startfile="."
            else:
                startfile=self.recentFiles[0]

        if not os.path.isfile(startfile) and not os.path.isdir(startfile):
            startfile = os.getcwd()
        filename = str(QFileDialog.getOpenFileName(self,"Open Data File",startfile,
        'Tab-delimited files (*.tab *.txt)\nSmiles files (*.smi)\nC4.5 files (*.data)\nAssistant files (*.dat)\nRetis files (*.rda *.rdo)\nBasket files (*.basket)\nAll files(*.*)'))

        if filename == "": return
        if filename in self.recentFiles: self.recentFiles.remove(filename)
        self.recentFiles.insert(0, filename)
        self.setFileList()

        self.openFile(self.recentFiles[0], 0, self.symbolDK, self.symbolDC)


    # Open a file, create data from it and send it over the data channel
    def openFile(self, fn, throughReload, DK=None, DC=None):
        self.error()
        self.warning()

        self.closeContext()
        self.loadedFile = ""
        
        if not fn or fn == "(none)" or not os.path.isfile(fn):
            self.send("Examples", None)
            self.send("Attribute Definitions", None)
            self.infoa.setText("No data loaded")
            self.infob.setText("")
            self.warnings.setText("")
            return
            
        self.symbolDK = self.symbolDC = ""
        self.openContext("", fn)

        self.loadedFile = ""

        argdict = {"createNewOn": 3-self.createNewOn}
        if DK:
            argdict["DK"] = str(DK)
        if DC:
            argdict["DC"] = str(DC)
        
        data = None
        ##scPA
        #Choose the appropriate dataLoader
        if os.path.splitext(fn)[1][1:].lower() == "smi":
            dataLoader =  dataUtilities.loadSMI
        else:
            dataLoader = dataUtilities.DataTable
            argdict["removeDuplicatedVars"] = True
        ##ecPA
        try:
            data = call(dataLoader, fn, **argdict)
            self.loadedFile = fn
        except Exception, (errValue):
            # If an exception occured, it might be that it was because 2 couses we assume to still load the data:
            #   1) The dataset is being loaded as a txt but it's actually a .tab
            #   2) All the values of an attribute are defined as unknown
            # in which cases it's a warining, and we want to still load the dada.
            #Since an exception was already detected in line 11, we need to be sure it's one of these cases, and reload
            #   the data not using the mehtod at line 11.
            if "is being loaded as" in str(errValue) or 'cannot determine type for attribute' in str(errValue):
                try:
                    data = dataLoader(fn, **argdict)
                    self.warning(0, str(errValue))
                except:
                    pass
            if data is None:
                self.error(str(errValue))
                self.dataDomain = None
                self.infoa.setText('No data loaded due to an error')
                self.infob.setText("")
                self.warnings.setText("")
                return

        ##scPA 
        # DataTable attributes:  attributeLoadStatus,duplicatedFixedVars,unnamedFixedVars,removedVars     
        WmsgRem = data.removedVars
        attributeLoadStatus = data.attributeLoadStatus
        Wmsg = data.duplicatedFixedVars
        WmsgUnnamed = data.unnamedFixedVars
        Wdetail = ""

        if len(WmsgUnnamed):
            Wdetail += "  -Renamed Attributes (They had no name)\n"
            for col in WmsgUnnamed:
                Wdetail += "    Column "+str(col)+" -> "+data.domain[col].name+"\n"

        if len(WmsgRem):
            Wdetail += "  -Deleted Attributes (They were the same)\n"
            for attr in WmsgRem:
                Wdetail += "    " + str(attr) + "  ("+ str(WmsgRem[attr]) +" attributes)\n"

        if len(Wmsg):
            Wdetail += "  -Renamed Attributes (They had the same name)\n"
            for attr in Wmsg:
                Wdetail += "    " + str(attr) + "  ("+ str(len(Wmsg[attr])) +" attributes)\n"
        if len(Wmsg) or len(WmsgRem) or len(WmsgUnnamed):
            QMessageBox.warning( None, "Fixed Attributes", "Check the dataset for consistency.\n\n" + Wdetail+"\n" , QMessageBox.Ok)
        ##ecPA
                
        self.dataDomain = data.domain

        self.infoa.setText('%d example(s), ' % len(data) + '%d attribute(s), ' % len(data.domain.attributes) + '%d meta attribute(s).' % len(data.domain.getmetas()))
        cl = data.domain.classVar
        if cl:
            if cl.varType == orange.VarTypes.Continuous:
                    self.infob.setText('Regression; Numerical class.')
            elif cl.varType == orange.VarTypes.Discrete:
                    self.infob.setText('Classification; Discrete class with %d value(s).' % len(cl.values))
            else:
                self.infob.setText("Class neither descrete nor continuous.")
        else:
            self.infob.setText("Data without a dependent variable.")

        warnings = ""
        metas = data.domain.getmetas()
        for status, messageUsed, messageNotUsed in [
                                (orange.Variable.MakeStatus.Incompatible,
                                 "",
                                 "The following attributes already existed but had a different order of values, so new attributes needed to be created"),
                                (orange.Variable.MakeStatus.NoRecognizedValues,
                                 "The following attributes were reused although they share no common values with the existing attribute of the same names",
                                 "The following attributes were not reused since they share no common values with the existing attribute of the same names"),
                                (orange.Variable.MakeStatus.MissingValues,
                                 "The following attribute(s) were reused although some values needed to be added",
                                 "The following attribute(s) were not reused since they miss some values")
                                ]:
            if self.createNewOn > status:
                message = messageUsed
            else:
                message = messageNotUsed
            if not message:
                continue
            attrs = [attr.name for attr, stat in zip(data.domain, data.attributeLoadStatus) if stat == status] \
                  + [attr.name for id, attr in metas.items() if data.metaAttributeLoadStatus.get(id, -99) == status]
            if attrs:
                warnings += "<li>%s: %s</li>" % (message, ", ".join(attrs))

        self.warnings.setText(warnings)
        qApp.processEvents()
        self.adjustSize()

        # make new data and send it
        fName = os.path.split(fn)[1]
        if "." in fName:
            data.name = fName[:fName.rfind('.')]
        else:
            data.name = fName

        #self.send("Examples", None)
        #self.send("Attribute Definitions", None)

        if len(data.domain.getmetas()) != 0: 
            QMessageBox.warning( None, "File - Meta Attributes", "There are meta-attributes present in the dataset.\nThe presence of meta-Attributes in datasets used with Learners/Classifiers\nrequires the use of considerably more ram memory!" , QMessageBox.Ok)   
        self.send("Examples", data)
        self.send("Attribute Definitions", data.domain)


if __name__ == "__main__":
    a=QApplication(sys.argv)
    owf=OWFile()
    owf.activateLoadedSettings()
    a.setMainWidget(owf)
    owf.show()
    a.exec_loop()
##    owf.openFile(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\Natasa Debeljak\12631098-A2flipped.txt", 1)
    owf.saveSettings()
