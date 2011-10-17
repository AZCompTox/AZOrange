"""
<name>Parameter Optimizer</name>
<description>Automatic optimization of learners parameters</description>
<icon>icons/Opt.png</icon>
<contact>Pedro Rafael Almeida</contact>
<priority>7</priority>
"""
import string
from AZutilities import dataUtilities
from OWWidget import *
import OWGUI
import orange
import AZOrangeConfig as AZOC
import os, types
from AZutilities import paramOptUtilities
from AZutilities import miscUtilities
import time
import AZLearnersParamsConfig
#from qttable import *
#from qt import QColor
#from qt import Qt
from copy import deepcopy
#import qt 
version = 9
class OWParamOpt(OWWidget):

    def __init__(self, parent=None, signalManager = None, name='ParamOptimizer'):
        OWWidget.__init__(self, parent, signalManager, name, 1)
        # Define the input and output channels
        self.inputs = [("Classified Examples", ExampleTable, self.setData), ("Learner", orange.Learner, self.setLearner)]
        self.outputs = [("Learner - Tuned", orange.Learner), ("Examples - Optimization Steps", ExampleTable)]

        self.requiredParamVer = 12
        self.name = name
	self.dataset = None
	self.learner = None
        self.optimizer = paramOptUtilities.Appspack()
        self.verbose = 0
        self.tunedPars = None
        self.intRes = None
        
        self.paramsNames=["Name","Optimize","Lower Limit","Upper Limit","Distribution","Step","Default","Actual Learner parameter"]

        self.learnerType = None
        self.parameters = None
        self.nParameters = 0

        self.nFolds = 5
        self.SMethod = 1
        self.execEnv = 0 # Serial
        self.CMethod = 0
        self.RMethod = 0  
        self.GUIparams = {}
        self.OptimizeChBox = {}
        self.DistCombo = {}
        self.DefaultCombo = {}
        self.UseGridSearch = False
        self.nInnerPoints = 5

        #Define Evaluation Methods: ["LabelName", "EvaluateFunction", "True=Best is the Min  |  False=Best is the Max"]
        self.CMethods = [("CA", "AZutilities.evalUtilities.CA",False)]

        self.RMethods = [("RMSE", "AZutilities.evalUtilities.RMSE", True)]
        #                 ("R^2", "AZutilities.evalUtilities.Rsqrt", False)]
        self.SMethods = [("Leave-One-Out", 0),
                         ("Cross Validation", 1)]

        self.execEnvs = AZOC.OWParamOptExecEnvs
       
        self.defineGUI()
        

    def defineGUI(self):

        self.sBox = OWGUI.widgetBox(self.controlArea, "Execution environment")
        itms = [e[0] for e in self.execEnvs]
        OWGUI.radioButtonsInBox(self.sBox, self, "execEnv", btnLabels=itms)

        boxGrid = OWGUI.widgetBox(self.controlArea,'Initial Point for the Optimizer')
        OWGUI.checkBox(boxGrid, self, 'UseGridSearch','Use Grid-Search',  tooltip='Use Grid-Search to find the best initial point to start the optimization.<br>If not checked, the midrange point will be used as optimization initial point.')
        OWGUI.spin(boxGrid, self, 'nInnerPoints', 1, 100, step=1, label='    Number of inner points:', tooltip='Number of points to break down each variable to evaluate the initial point.<br>It will evaluate nInnerPoints^nOptimizedVars')

        OWGUI.separator(self.controlArea)

        self.sBox = OWGUI.widgetBox(self.controlArea, "Sampling  Method")
        itms = [e[0] for e in self.SMethods]
        self.comboDistItems = ["Continuous","Power2","By Step","Specific Values"]        
        OWGUI.radioButtonsInBox(self.sBox, self, "SMethod", btnLabels=itms)

        hBox = OWGUI.widgetBox(OWGUI.indentedBox(self.sBox))
        #QWidget(hBox).setFixedSize(19, 8)
        OWGUI.spin(hBox, self, 'nFolds', 2, 100, step=1, label='Number of Folds:  ')


        OWGUI.separator(self.controlArea)

        box2 = OWGUI.widgetBox(self.controlArea,'Evaluation Method')
        width = 150
        itms = [e[0] for e in self.CMethods]
        OWGUI.comboBox(box2, self, 'CMethod', items=itms, label='For Classifiers:', labelWidth=width, orientation='horizontal',
                       tooltip='Method used for evaluation in case of classifiers.')
        itms = [e[0] for e in self.RMethods]
        OWGUI.comboBox(box2, self, 'RMethod', items=itms, label='For Regressors:', labelWidth=width, 
                       orientation='horizontal', tooltip='Method used for evaluation in case of regressors.')

        OWGUI.separator(self.controlArea)
        OWGUI.button(self.controlArea, self,"&Reload Defaults ", callback=self.reloadDefaults)
        OWGUI.separator(self.controlArea)


        #OWGUI.separator(self.controlArea, height=24)

        infoBox = OWGUI.widgetBox(self.controlArea, "Optimizer status")
        self.infoStatus = OWGUI.label(infoBox,self,'Waiting for inputs...')                    
        self.infoPars = OWGUI.label(infoBox,self,'')
        self.infoRes = OWGUI.label(infoBox,self,'')
        OWGUI.label(infoBox,self,'')
        self.infoErr = OWGUI.label(infoBox,self,'')
        OWGUI.separator(self.controlArea)
        OWGUI.button(self.controlArea, self,"&Apply Settings ", callback=self.optimizeParameters)
 

        # Main area GUI
        import sip
        sip.delete(self.mainArea.layout())
        self.mainLayout = QGridLayout(self.mainArea)

        mainRight = OWGUI.widgetBox(self.mainArea, "Parameters Configuration")
        mainRight.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding))

        self.paramsTable = OWGUI.table(self.mainArea, rows = 0, columns = 0, selectionMode = QTableWidget.MultiSelection, addToLayout = 0)
        #self.paramsTable.setLeftMargin(0)
        self.paramsTable.verticalHeader().hide()
        self.paramsTable.setSelectionMode(QTableWidget.NoSelection)
        self.paramsTable.setColumnCount(len(self.paramsNames))

        #for i, m in enumerate(self.paramsNames):
        #    self.paramsTable.setColumnStretchable(i, 0)
        #    header.setLabel(i, m)
        self.paramsTable.setHorizontalHeaderLabels(self.paramsNames)
        #self.mainLayout.setColumnStretch(1, 100)
        #self.mainLayout.setRowStretch(2, 100)
        self.mainLayout.addWidget(self.paramsTable, 0, 0, 1,2)
        self.mainLayout.addWidget(OWGUI.label(mainRight,self,'Red - Parameter selected to be optimized\r\nGreen - Parameter optimized\r\nBlack - Parameter will not be optimized'),1,0) 
        self.mainLayout.addWidget(OWGUI.label(mainRight,self,'N_EX - Number of Examples in dataset\r\nN_ATTR - Number of attributes in dataset'),1,1)
         
        self.adjustSize()
        self.create()


    def reloadDefaults(self):
        if self.learner:
            self.originalParameters = eval("AZLearnersParamsConfig." + self.learnerType)
            self.parameters = deepcopy(self.originalParameters)
            self.updateTable()


    def clearTable(self):
        self.learnerType = None
        self.parameters = None
        for row in range(self.paramsTable.rowCount()):
            self.paramsTable.removeRow(row)
        self.paramsTable.setRowCount(0)

    def getRangeParsAndDistType(self, parametersOrig ,param):

        """Get from the parametersOrig the Distribution type and the Range for the parameter named as the input var 'param'
        """
        N_EX = 5    #Just a number to test the avaluation of the expressions on atributes definition
        N_ATTR = 5  #Just a number to test the avaluation of the expressions on atributes definition

        if parametersOrig[param][1] == "interval":
                distType = 0
        elif parametersOrig[param][1] == "values":
                if "power2Range" in parametersOrig[param][2]:
                        distType = 1
                elif "Range" in parametersOrig[param][2]:
                        distType = 2
                else:
                        distType = 3
        else:
                self.setErrors("Invalid keyword in Configuration File")
                return [None,None]

        if distType in (1,2):
                txt = parametersOrig[param][2]
                RangePars = txt[txt.find("Range")+6:txt.find(")",txt.find("Range"))].split(",")
                if len(RangePars)<3:
                    RangePars.append("1")  #the default step size
                try:
                    if distType == 1 and miscUtilities.power2Range(eval(RangePars[0]),eval(RangePars[1]),eval(RangePars[2])) != eval(parametersOrig[param][2]):
                        return [[],3]
                    elif distType == 2 and miscUtilities.Range(eval(RangePars[0]),eval(RangePars[1]),eval(RangePars[2])) != eval(parametersOrig[param][2]):
                        return [[],3]
                except:
                    return [[],3]
                    
        elif distType == 0:
                txt = parametersOrig[param][2]
                RangePars = [txt[1:-1].split(",")[0],txt[1:-1].split(",")[-1]]
        else:#distType==3
                txt = parametersOrig[param][2]
                RangePars = [x.strip() for x in txt[1:-1].split(" , ")]
                if eval(parametersOrig[param][0]) in types.StringTypes or \
                         (type(eval(parametersOrig[param][0]))==types.ListType and (type(eval(parametersOrig[param][0])[0]) in types.StringTypes)):
                    newPars = []
                    for par in RangePars:
                        if par[0]==par[-1] and par[0] in ("'",'"'):
                            newPars.append(par[1:-1])
                        else:
                            newPars.append(par)
                    RangePars = newPars
        RangePars=[x.strip() for x in RangePars]
        return [RangePars,distType]


    def setCellText(self,table,row,col,text):
        table.removeCellWidget( row, col)
        it = QTableWidgetItem()
        it.setFlags(Qt.ItemIsEnabled | (Qt.ItemIsSelectable or Qt.NoItemFlags))
        it.setTextAlignment(Qt.AlignRight)
        it.setText(text)
        table.setItem(row, col, it)
        return it

    def setCellComboBox(self,table,row,col,items):
        table.removeCellWidget( row, col)
        #it = OWGUI.comboBox(None, self, None, items = items, tooltip = "")
        combo = QComboBox()
        combo.addItems([unicode(i) for i in items])
        table.setCellWidget(row,col,combo)
        return combo
 
    def updateTable(self):
        """Updates the GUI table parameters with the ones specified on self.parameters variable
           Also updates the optimized parameters column if the self.optimizer is optimized with success
        """
        if not self.parameters:
            return   

        #self.parameters["NAME"]: [ParameterType, valuesRangeType,valuesRange,valuesAlias,Default,Optimize, EditableRange] 
        #self.paramsNames: ["Name","Optimize","Lower Limit","Upper Limit","Distribution","Step","Default"]
        #self.comboDistItems: ["Continuous","Power2","By Step","Specific Values"]
        self.nParameters = len(self.parameters)
        self.paramsTable.setRowCount(len(self.parameters))      
        self.OptimizeChBox = {}
        self.DistCombo = {}
        self.DefaultCombo = {}
        #QcomboDistItems = QStringList()
        #[QcomboDistItems.append(str(x)) for x in self.comboDistItems]
        for row, param in enumerate(self.parameters):
            origRangePars,origDistType = self.getRangeParsAndDistType(self.originalParameters,param)
            RangePars,distType = self.getRangeParsAndDistType(self.parameters,param)
            if RangePars==None or distType==None:
                self.setErrors("Parameter "+str(param)+" has undefined distribution!","WARNING")
                return
            for col, m in enumerate(self.paramsNames):
                if col==0:#Name
                    self.setCellText(self.paramsTable, row, col, str(param))
                elif col==1:#Optimize Flag
                    self.paramsTable.removeCellWidget( row, col)
                    self.OptimizeChBox[str(param)] = QCheckBox()
                    self.OptimizeChBox[str(param)].setTristate(False)
                    self.OptimizeChBox[str(param)].setCheckState(int(self.parameters[param][5])*2)   #0-Unchecked    1-Partially checked    2-Checked
                    self.paramsTable.setCellWidget(row,col,self.OptimizeChBox[str(param)])                    
                    if not self.parameters[param][5]: self.nParameters -= 1
                elif col==2:#Lower
                    if distType != 3:
                        self.setCellText(self.paramsTable, row, col, str(RangePars[0]))
                    else:
                        self.setCellText(self.paramsTable, row, col, "")
                elif col==3:#Upper
                    if distType != 3:
                        self.setCellText(self.paramsTable, row, col, str(RangePars[1]))
                    else:
                        self.setCellText(self.paramsTable, row, col, "")
                elif col==4:#Distribution Type
                    self.paramsTable.removeCellWidget( row, col)
                    if "Combo" in str(type(self.paramsTable.cellWidget(row,6))) and not self.parameters[param][6]:
                        self.setCellText(self.paramsTable, row, col, self.comboDistItems[-1]) 
                    else:
                        if distType == len(self.comboDistItems)-1:
                            self.DistCombo[str(param)] = self.setCellComboBox(self.paramsTable,row,col,[str(x) for x in self.comboDistItems])
                        else:
                            self.DistCombo[str(param)] = self.setCellComboBox(self.paramsTable,row,col,[str(x) for x in self.comboDistItems][:-1])                            
                        self.DistCombo[str(param)].setCurrentIndex(distType)
                    #if distType==3:
                     #   if self.paramsTable.item(row,7) != None and hasattr(self.paramsTable.item(row,7),"setToolTip"):
                      #      print "Seted"
                       #     self.paramsTable.item(row,7).setToolTip("is This OK?")
                        #    self.paramsTable.updateCell(row, col)
                        #QToolTip.add(self.paramsTable,self.paramsTable.cellGeometry(row+1, col-2),str(self.originalParameters[param][2]))
                    
                elif col==5: #Step
                    if distType in (0,3):
                        cellItem = ""
                    elif distType in (1,2):
                        cellItem = str(RangePars[2])
                    else:
                        cellItem = str(self.parameters[param][2])
                    self.setCellText(self.paramsTable, row, col, cellItem)
                    if distType==3 and self.originalParameters[param][7]:
                        if self.paramsTable.item(row,col) != None and hasattr(self.paramsTable.item(row,col),"setToolTip"):
                            toolTipTXT = self.originalParameters[param][7]
                            self.paramsTable.item(row,col).setToolTip("Specific Values from Config File:\n"+toolTipTXT)
                    #self.paramsTable.updateCell(row, col)
                    self.paramsTable.update()

                elif col ==6:#Default
                    #QComboDef = QStringList()
                    self.paramsTable.removeCellWidget( row, col)
                    if distType == 3 and len(origRangePars) > 0 and len(self.parameters[param][3])==len(origRangePars) and not self.parameters[param][6]:
                        if self.parameters[param][4].strip() in origRangePars:
                            #[QComboDef.append(str(x)) for x in self.parameters[param][3]]
                            self.DefaultCombo[str(param)] = self.setCellComboBox(self.paramsTable,row,col,[str(x) for x in self.parameters[param][3]])
                            self.DefaultCombo[str(param)].setCurrentIndex(origRangePars.index(self.parameters[param][4].strip()))
                        else:#The alias are present, but could not find a match
                            if self.verbose > 0: self.setErrors("No match found for the Default parameter of " + str(param) + " to resolve the alias","WARNING")
                            self.DefaultCombo[str(param)] = self.setCellComboBox(self.paramsTable,row,col,[str(x) for x in self.parameters[param][3]]+[str(self.parameters[param][4].strip())])
                            self.DefaultCombo[str(param)].setCurrentIndex(self.DefaultCombo[str(param)].count-1)
                            #self.setCellText(self.paramsTable, row, col, str(self.parameters[param][4]).strip())
                    elif len(origRangePars)>1 and distType == 3 and not self.parameters[param][6]: 
                        if self.parameters[param][4].strip() in origRangePars: 
                            self.DefaultCombo[str(param)] = self.setCellComboBox(self.paramsTable,row,col,[str(x) for x in origRangePars])
                            self.DefaultCombo[str(param)].setCurrentIndex(origRangePars.index(self.parameters[param][4].strip()))
                        else:   
                            if self.verbose > 0: self.setErrors("The default parameter for " + str(param) + " is not listed in the values list","WARNING")
                            self.DefaultCombo[str(param)] = self.setCellComboBox(self.paramsTable,row,col,[str(x) for x in origRangePars]+[str(self.parameters[param][4].strip())])
                            self.DefaultCombo[str(param)].setCurrentItem(self.DefaultCombo[str(param)].count-1)
                    else:
                        self.setCellText(self.paramsTable, row, col, str(self.parameters[param][4]).strip())
                elif col==7: #Optimization
                   self.paramsTable.removeCellWidget( row, col)
                   if self.learner:
                        if hasattr(self.learner,str(param)):
                            aParam = str(eval("self.learner."+str(param))) 
                            # Check if an alias for this value exists
                            if aParam in origRangePars:
                                aPidx = origRangePars.index(aParam)
                                if len(self.parameters[param][3]) > aPidx:
                                    aParam = self.parameters[param][3][aPidx]
                            #cellItem = ColorTableItem(self.paramsTable,aParam ,QTableWidgetItem.Never)
                            cellItem = self.setCellText(self.paramsTable, row, col, aParam) 
                            if hasattr(self.learner,"optimized") and self.learner.optimized and  self.parameters[param][5]:
                                cellItem.setForeground(QBrush(QColor("Green")))
                            elif not self.parameters[param][5]:
                                cellItem.setForeground(QBrush(QColor("black")))
                            else:
                                cellItem.setForeground(QBrush(QColor("red")))    
                        else:
                            self.setCellText(self.paramsTable, row, col, "N/A")
                        #self.paramsTable.updateCell(row, col)
                        self.paramsTable.update()
                else:
                    pass
        #for i in range(len(self.paramsNames)):
        #    self.paramsTable.adjustColumn(i)
        self.paramsTable.resizeColumnsToContents()
        self.paramsTable.setColumnWidth(5,30)

    def updateParametersFromTable(self):
        """Updates the parameters of the optimizer with the ones present on GUI table
           Returns True if all OK
           Returns False if Errors occurred
        """
        #self.paramsNames: ["Name","Optimize","Lower Limit","Upper Limit","Distribution","Step","Default"]
        #self.comboDistItems: ["Continuous","Power2","By Step","Specific Values"]
        if self.paramsTable.columnCount() < 7:
            self.setErrors("Wrong number of columns in table!")
            return False
        RangePars,distType = [None,None]
        self.nParameters = len(self.parameters)
        for row in range(self.paramsTable.rowCount()):
            for col in range(7):
                if col == 0:#Name
                    name=str(self.paramsTable.item(row,col).text()).strip()
                    RangePars,distType = self.getRangeParsAndDistType(self.originalParameters,name)  
                    if RangePars==None or distType==None:
                        self.setErrors("It was not possible to identify the range parameters")
                        return False
                elif col == 1:#Optimize
                    if self.paramsTable.cellWidget(row,col).checkState()==2:
                        optimize=True
                    else:
                        optimize = False
                        self.nParameters -= 1
                elif col == 2:#Lower Limit
                    Llimit = str(self.paramsTable.item(row,col).text()).strip()
                elif col == 3:#Upper Limit
                    Ulimit = str(self.paramsTable.item(row,col).text()).strip()
                elif col == 4:#Distribution
                    if "Combo" in str(type(self.paramsTable.cellWidget(row,col))):
                        dist =  str(self.paramsTable.cellWidget(row,col).currentIndex())
                    else:
                        dist = self.comboDistItems.index(self.paramsTable.item(row,col).text())
                elif col == 5:#Step
                    step = str(self.paramsTable.item(row,col).text()).strip()
                elif col == 6:#Default
                    if "Combo" in str(type(self.paramsTable.cellWidget(row,col))):
                        default =  self.paramsTable.cellWidget(row,col).currentIndex()
                        if len(self.parameters[name][3]) <= default:  #the parameter did not have an alias, the text is the value
                            default = str(self.paramsTable.cellWidget(row,col).currentText()).strip()
                        else:
                            default = RangePars[default]
                    else:#is string
                        default = str(self.paramsTable.item(row,col).text()).strip()

#['types.StringType', 'values', "['kernel' , 'pls1' , 'simpls']", ['Kernel', 'PLS1', 'SimPLS'], "'simpls'", True, False]

            if optimize:
                self.parameters[name][5] = True
                if "Combo" in str(type(self.paramsTable.cellWidget(row,4))):
                    comboDistType = self.paramsTable.cellWidget(row,4).currentIndex()
                else:
                    comboDistType = self.comboDistItems.index(self.paramsTable.item(row,4).text())
                if comboDistType==0:     #Continuous
                    if not miscUtilities.isNumber(Llimit):
                        if "N_EX" not in Llimit and "N_ATTR" not in Llimit:
                            QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has invalid Lower limit",QMessageBox.Ok)               
                            return False
                    if not miscUtilities.isNumber(Ulimit):
                        if "N_EX" not in Ulimit and "N_ATTR" not in Ulimit:
                            QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has invalid Upper limit",QMessageBox.Ok)
                            return False
                    self.parameters[name][1] = "interval"
                    self.parameters[name][2] = "[" + Llimit + " , " + Ulimit + "]"
                    self.parameters[name][3] = ""
                    self.parameters[name][4] = default
                elif comboDistType==1:   #Power2
                    if not miscUtilities.isNumber(Llimit):
                        if "N_EX" not in Llimit and "N_ATTR" not in Llimit:
                            QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has invalid Lower limit",QMessageBox.Ok)
                            return False
                    if not miscUtilities.isNumber(Ulimit):
                        if "N_EX" not in Ulimit and "N_ATTR" not in Ulimit:
                            QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has invalid Upper limit",QMessageBox.Ok)
                            return False
                    if not miscUtilities.isNumber(step):
                        QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has an invalid step value",QMessageBox.Ok)
                        return False
                    self.parameters[name][1] = "values"
                    self.parameters[name][2] = "miscUtilities.power2Range(" + Llimit + "," + Ulimit + "," + step + ")"
                    self.parameters[name][3] = ""
                    self.parameters[name][4] = default
                elif comboDistType==2:   #By Step
                    if not miscUtilities.isNumber(Llimit):
                        if "N_EX" not in Llimit and "N_ATTR" not in Llimit:
                            QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has invalid Lower limit",QMessageBox.Ok)
                            return False
                    if not miscUtilities.isNumber(Ulimit):
                        if "N_EX" not in Ulimit and "N_ATTR" not in Ulimit:
                            QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has invalid Upper limit",QMessageBox.Ok)
                            return False
                    if not miscUtilities.isNumber(step):
                        QMessageBox.warning(self,"Invalid parameter","Parameter "+name+" has an invalid step value",QMessageBox.Ok)
                        return False
                    self.parameters[name][1] = "values"
                    self.parameters[name][2] = "miscUtilities.Range(" + Llimit + "," + Ulimit + "," + step + ")"
                    self.parameters[name][3] = ""
                    self.parameters[name][4] = default
                else:                                                 #Specific Values
                    #The 'Specific Values' refere to the parameters specified in the original AZLearnersParamsConfig.py file
                    self.parameters[name][1] = self.originalParameters[name][1]
                    self.parameters[name][2] = self.originalParameters[name][2]
                    self.parameters[name][3] = self.originalParameters[name][3]   
                    self.parameters[name][4] = default 
            else:
                self.parameters[name][5] = False
                #self.parameters[name][1]='values'
                if "Combo" in str(type(self.paramsTable.cellWidget(row,6))):
                    default =  self.paramsTable.cellWidget(row,6).currentIndex()
                    if len(RangePars) <= default:  #the parameter did not have an alias, the text is the value
                        #self.parameters[name][2]='[' +  str(self.paramsTable.cellWidget(row,6).currentText()).strip() + ']'
                        self.parameters[name][4]=str(self.paramsTable.cellWidget(row,6).currentText()).strip()
                    else:
                        #self.parameters[name][2]='[' + RangePars[default]  + ']'
                        self.parameters[name][4]= RangePars[default]
                else:#is string
                    #self.parameters[name][2]='[' +  str(self.paramsTable.item(row,6).text()).strip() + ']'
                    self.parameters[name][4]=str(self.paramsTable.item(row,6).text()).strip()
        return True
            

    def showEvent(self, ev):
        self.updateInfo()

    def onDeleteWidget(self):
	#self.linksOut.clear()
	if self.dataset:
	    del self.dataset
        if self.intRes:
            del self.intRes

    def updateInfo(self):
        self.infoStatus.setText('')
        self.infoPars.setText('')
        self.infoRes.setText('')
        
        
        if not self.learner and not self.dataset:
            self.infoStatus.setText('Waiting for inputs...')
            self.updateTable()
            self.adjustSize()
            return
        elif not self.learner:
            self.infoStatus.setText('Waiting for learner...')
            self.updateTable()
            self.adjustSize()
            return
        elif not self.dataset:
            self.infoStatus.setText('Waiting for data...')
            self.updateTable()
            self.adjustSize()
            return
        

        if type(self.tunedPars) != types.ListType or (not hasattr(self.learner,"optimized")):
            self.setErrors("Some error occurred!\nTuned parameters: \r\n  "+str(self.tunedPars), "WARNING")
            self.infoPars.setText("Check the returned Parametrs\nin the output window")
            if self.verbose > 0: self.setErrors("Tuned parameters:" + str(self.tunedPars),"WARNING")
        else:
            if self.learner.optimized == True:
                if len(self.intRes) <= 2:
                    self.setErrors("It was detected that only the default and initial points were evaluated.\nIt seems that there was no optimization!\nPlease check for any other printed errors/warnings.","WARNING")
                else:
                    if self.verbose > 0: print "Learner Optimized!"
                    self.infoStatus.setText('== Learner Optimized! ==')
            else:
                if self.verbose > 0: self.setErrors("optimization Flag = "+str(self.learner.optimized),"WARNING")
                self.setErrors('Optimization seems to be done, but repective flag is not checked ('+str(self.learner.optimized)+')',"WARNING")
            LPars = ""
            self.infoPars.setText("")
            if self.verbose > 0:
                for parName in self.tunedPars[1]:
                    if hasattr(self.learner,parName):
                        LPars += parName + " = "
                        if self.parameters != None and self.parameters[str(parName)][1]=="values" and len(self.parameters[str(parName)])>=3 and self.parameters[str(parName)][3]!="":
                            try:
                                LPars += str(self.parameters[str(parName)][3][  eval(self.parameters[str(parName)][2]).index(eval("self.learner."+parName)) ])
                            except:
                                LPars += str(eval("self.learner."+parName)) +" (Could not find parameter value alias)"
                                print "Missing alias for: ",eval("self.learner."+parName)
                        else:
                            LPars += str(eval("self.learner."+parName))
                        LPars += " ("+str(self.tunedPars[1][parName])+")\r\n\t"
                    else:
                        LPars += parName+ " = NoParameter" + " ("+str(self.tunedPars[1][parName])+")\r\n\t"        
                self.infoPars.setText("Check the returned Parametrs\nin the output window")
                print "Learner parameters (tuned):\r\n\t" + LPars
            
            self.infoRes.setText("Best optimization result = " + str(self.tunedPars[0]))
        self.adjustSize()
        self.updateTable()

    def setData(self, dataset):
        if dataset:
            self.dataset = dataset
        else:
            self.dataset = None
        self.optimizeParameters()


    def setLearner(self, learner):
        self.infoStatus.setText('')
        self.infoPars.setText('')
        self.infoRes.setText('')
        self.clearTable()
 
        if learner:
            self.learner = learner
        else:
            self.learner = None
            return

        # check if learner is defined in AZLearnersParamsConfig
        self.learnerType = str(self.learner).split()[0]
        if self.learnerType[0] == "<":
            self.learnerType = self.learnerType[1:]
        if self.learnerType.rfind('.') >=0:
            self.learnerType = self.learnerType[self.learnerType.rfind('.')+1:]
        if AZLearnersParamsConfig.version != self.requiredParamVer:
            self.setErrors("The version of AZLearnersParamsConfig.py is not the correct for this widget. Please use the version "+str(self.requiredParamVer))
            self.learnerType = None
            self.parameters = None
            return

        if not hasattr(AZLearnersParamsConfig, self.learnerType):
            self.setErrors("The learner " + str(self.learnerType) +" is not compatible with the Optimizer!")
            self.learnerType = None
            self.parameters = None
            return
        else:
            self.originalParameters = eval("AZLearnersParamsConfig." + self.learnerType)
            self.parameters = deepcopy(self.originalParameters)
        for parameter in self.parameters:
            if len(self.parameters[parameter]) != 8:
                self.setErrors("The version of AZLearnersParamsConfig.py is not the correct for this widget")
                self.parameters = None
                return

        self.updateTable()  
        self.optimizeParameters()

    def clearErrors(self):
        self.warning()
        self.infoErr.setText('')

    def setErrors(self, msg, errType = "ERROR"):
        self.warning(0,msg)
        self.infoErr.setText("=== "+errType+" ===\n Please Check the output window\nfor more details.")
        print errType," (",time.asctime(),"):"
        print "  ",msg

    
    def optimizeParameters(self): 
        """ Sets up the input learner with tuned parameters  """

        self.clearErrors()
        self.tunedPars = None
        if hasattr(self.learner,"optimized"):
            self.learner.optimized = False

        if not self.learner:
            self.send("Learner - Tuned", None)
            self.send("Examples - Optimization Steps", None)
            self.updateInfo()
            return

        # Apply the parameters var with values  on configuration table of GUI (user could have changed them!)
        if not self.updateParametersFromTable():
            return
   
        if not self.dataset:
            self.dataset = None
            self.send("Learner - Tuned", None)
            self.send("Examples - Optimization Steps", None)
            self.updateInfo()
            return

        # Progess Bar 1
        optSteps = 3
        progress1 = QProgressDialog("Gathering data and configuring the optimizer...", "Cancel", 0, optSteps, self,Qt.Dialog)#, "progress", True )
        progress1.setWindowModality(Qt.WindowModal)
        bar1 = QProgressBar(progress1)
        bar1.show()
        progress1.setBar(bar1)
        #progress1.setTotalSteps(optSteps)
        progress1.setMinimumDuration(0)
        progress1.forceShow()
        progress1.setValue(0)
        time.sleep(0.1)
        progress1.setValue(0)

        # Create path for running the optimizer
        randNr = random.randint(0,10000)
        if self.execEnv == 0:
            scratchdir = miscUtilities.createScratchDir(desc = "OWParamOpt_Serial")
        else:
            scratchdir = miscUtilities.createScratchDir(desc ="OWParamOpt_MPI", baseDir = AZOC.NFS_SCRATCHDIR)
        # Save the dataset to the optimizer running path
        OrngFile = os.path.join(scratchdir,"OrngData.tab")
        orange.saveTabDelimited(OrngFile,self.dataset)
        # Advance Progress Bar
        progress1.setValue(1)
        # Define the evaluation method to use
        if self.dataset.domain.classVar.varType == orange.VarTypes.Continuous:
            fMin = self.RMethods[self.RMethod][2]
            evalM = self.RMethods[self.RMethod][1]
        else:
            fMin = self.CMethods[self.CMethod][2]
            evalM= self.CMethods[self.CMethod][1]
        try:
            if os.path.exists(os.path.join(scratchdir,"AZLearnersParamsConfig.py")):
                os.system("rm "+str(os.path.join(scratchdir,"AZLearnersParamsConfig.py"))) 
            paramFile=file(os.path.join(scratchdir,"AZLearnersParamsConfig.py"),"w")
            paramFile.write(self.learnerType + "= " + str(self.parameters)+"\r\n")
            paramFile.close()

            progress1.setValue(2)
            # Run the optimizer which will configure the input learner and aditionaly return [<minimum of objective function found>, <optimized parameters>]
            # Serial
            print "ENV:",self.execEnv
            if self.execEnv == 0:
                print "Executing the optimizer in serial mode on local machine"
                optPID = self.optimizer(learner=self.learner, dataSet=OrngFile, evaluateMethod = evalM , findMin=fMin, nFolds = self.nFolds, samplingMethod = self.SMethods[self.SMethod][1], runPath = scratchdir, verbose = self.verbose, externalControl = 1,useParameters = self.parameters, useGridSearchFirst = self.UseGridSearch,  gridSearchInnerPoints=self.nInnerPoints, np = None, machinefile = None, advancedMPIoptions = "",)
            # Local mpi
            elif self.execEnv == 1:
                print "Executing the optimizer in parallel mode on local machine"
                optPID = self.optimizer(learner=self.learner, dataSet=OrngFile, evaluateMethod = evalM , findMin=fMin, nFolds = self.nFolds, samplingMethod = self.SMethods[self.SMethod][1], runPath = scratchdir, verbose = self.verbose, externalControl = 1,useParameters = self.parameters, useGridSearchFirst = self.UseGridSearch,  gridSearchInnerPoints=self.nInnerPoints, machinefile = 0)
            # Sge Molndal
            elif self.execEnv == 2:
                print "Executing the optimizer in parallel mode in the batch queue on the sge"
                print "*****************runPath*****************"
                optPID = self.optimizer(learner=self.learner, dataSet=OrngFile, evaluateMethod = evalM , findMin=fMin, nFolds = self.nFolds, samplingMethod = self.SMethods[self.SMethod][1], runPath = scratchdir, verbose = self.verbose, externalControl = 1,useParameters = self.parameters, useGridSearchFirst = self.UseGridSearch,  gridSearchInnerPoints=self.nInnerPoints, np = 8,machinefile = "qsub")#, sgeEnv = "sge_seml")
            elif self.execEnv == 3:
                print "Executing the optimizer in parallel mode in the quick queue on the sge"
                print "*****************runPath*****************"
                optPID = self.optimizer(learner=self.learner, dataSet=OrngFile, evaluateMethod = evalM , findMin=fMin, nFolds = self.nFolds, samplingMethod = self.SMethods[self.SMethod][1], runPath = scratchdir, verbose = self.verbose, externalControl = 1,useParameters = self.parameters, useGridSearchFirst = self.UseGridSearch,  gridSearchInnerPoints=self.nInnerPoints, np = 8,machinefile = "qsub",queueType = "quick.q")#, sgeEnv = "sge_seml")
            else:
                print "No SGE Env. selected. Nothing will happen."
        except:
            progress1.close()
            self.updateInfo()
            self.setErrors("Some error(s) occurred during the optimization.\nCheck the "+str(scratchdir)+" and the output terminal for more information")
            self.send("Learner - Tuned", None)
            self.send("Examples - Optimization Steps", None)
            return

        progress1.setValue(3)

        if type(optPID)!=types.IntType:
            progress1.close()
            self.updateInfo()
            self.setErrors("Some error(s) occurred during optimization:\n"+str(optPID))
            self.send("Learner - Tuned", None)
            self.send("Examples - Optimization Steps", None)
            return


        progress1.close()

        # Progess Bar
        optSteps = (1+round((len(self.dataset)*len(self.dataset.domain.attributes)*self.nParameters)/1000))*8
        print "Learner optimization started at "+time.asctime()
        print "Optimization steps = ",int(optSteps)," (estimated to aprox. ",optSteps/2," seconds)"
        progress = QProgressDialog("Learner optimization started at "+time.asctime()+" ,please wait...", "Abort Optimization", 0,optSteps ,self,Qt.Dialog)#, "progress", True )
        progress.setWindowModality(Qt.WindowModal)
        bar = QProgressBar(progress)
        bar.show()
        progress.setBar(bar)
        #progress.setTotalSteps(optSteps)
        progress.setMinimumDuration(0)
        stepsDone = 0
        progress.setValue(stepsDone)
        progress.forceShow()
        #Loop waiting for the optimizer to finish
        while 1:
            if stepsDone < (progress.maximum()-1):
                progress.setValue(stepsDone)
                stepsDone+=1
                time.sleep(0.5)
            else:
                bar.setTextVisible(False)
                progress.setLabelText("The optimizer is taking longer than expected, please wait some more time...")
                stepsDone = 0
                progress.setValue(stepsDone)
                time.sleep(0.5)
            if progress.wasCanceled():
                if not self.optimizer.stop():
                    progress.setLabelText("Could not stop the optimizer! Please wait until it finish...")
                else:
                    self.setErrors("Learner optimization stopped by user at "+time.asctime(),"WARNING")
                    break
            if self.optimizer.isFinished():
                print "Learner optimization finished at "+time.asctime()
                break
        progress.setValue(progress.maximum()-1)
        time.sleep(0.5)
        progress.setValue(progress.maximum())   
        self.tunedPars = self.optimizer.tunedParameters 
        if self.verbose > 0:
            if self.optimizer.usedMPI:
                print "appspack version used in fact: MPI"
            else:
                print "appspack version used in fact: SERIAL"
        if type(self.tunedPars) != types.ListType or self.learner.optimized == False:
            self.send("Learner - Tuned", None)
            self.send("Examples - Optimization Steps", None)
        else:
            self.send("Learner - Tuned", self.learner)
            self.intRes = dataUtilities.DataTable(scratchdir+"/optimizationLog.txt") 
            self.send("Examples - Optimization Steps", self.intRes)
        self.updateInfo()

        if self.verbose == 0:
            miscUtilities.removeDir(scratchdir)
        else:
            self.setErrors("The directory " + str(scratchdir) + " was not deleted because verbose flag is ON","DEBUG")

class ProgressBar:
    def __init__(self, widget, iterations):
        self.iter = iterations
        self.widget = widget
        self.count = 0
        self.widget.progressBarInit()
    def advance(self):
        self.count += 1
        self.widget.progressBarSet(int(self.count*100/self.iter))
    def finish(self):
        self.widget.progressBarFinished()



if __name__ == "__main__": 
    appl = QApplication(sys.argv) 
    ow = OWParamOpt() 
    appl.setMainWidget(ow) 
    ow.show() 
    dataset = dataUtilities.DataTable('iris.tab') 
    ow.data(dataset) 
    appl.exec_loop()

