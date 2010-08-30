"""
<name>Rank</name>
<description>Ranks and filters attributes by their relevance.</description>
<icon>icons/Rank.png</icon>
<contact>Janez Demsar (janez.demsar(@at@)fri.uni-lj.si), changed By Pedro Almeida</contact>
<priority>1102</priority>
"""

from AZutilities import dataUtilities
from OWWidget import *
#from qttable import *
import OWGUI
from trainingMethods import AZorngRF


class OWRank(OWWidget):
    ##scPA       Added the MSE estimator that handles continuous class data              ##ecPA
    #             Added the RF estimator
    settingsList =  ["showDistributions","nDecimals", "reliefK", "reliefN", "nIntervals", "sortBy", "nSelected", "selectMethod", "distColorRgb","autoApply"]
    measures          = ["ReliefF", "Information Gain", "Gain Ratio", "Gini Gain", "Mean Squares Errors","RandomForest"]
    measuresShort     = ["ReliefF", "Inf. gain", "Gain ratio", "Gini", "MSE", "RF"]
    measuresAttrs     = ["computeReliefF", "computeInfoGain", "computeGainRatio", "computeGini", "computeMSE", "computeRF"]
    estimators        = [orange.MeasureAttribute_relief, orange.MeasureAttribute_info, orange.MeasureAttribute_gainRatio, orange.MeasureAttribute_gini, orange.MeasureAttribute_MSE, AZorngRF.RFLearner]
    handlesContinuous = [True, False, False, False,False,True]

    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "Rank")#,wantMainArea = 0)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Reduced Example Table", ExampleTable, Default + Single), ("ExampleTable Attributes", ExampleTable, NonDefault)]

        self.settingsList += self.measuresAttrs

        self.nDecimals = 3
        self.reliefK = 10
        self.reliefN = 20
        self.nIntervals = 4
        self.sortBy = 0
        self.selectMethod = 3
        self.nSelected = 5
        self.showDistributions = 1
        self.distColorRgb = (220,220,220, 255)
        self.distColor = QColor(*self.distColorRgb)
        self.autoApply = True

        for meas in self.measuresAttrs:
            setattr(self, meas, True)

        self.loadSettings()

        labelWidth = 80

        box = OWGUI.widgetBox(self.controlArea, "Measures", addSpace=True)
        for meas, valueName in zip(self.measures, self.measuresAttrs):
            OWGUI.checkBox(box, self, valueName, meas, callback=self.measuresChanged)
            if valueName == "computeReliefF":
                ibox = OWGUI.indentedBox(box)
                OWGUI.spin(ibox, self, "reliefK", 1, 20, label="Neighbours", labelWidth=labelWidth, orientation=0, callback=self.reliefChanged, callbackOnReturn = True)
                OWGUI.spin(ibox, self, "reliefN", 20, 100, label="Examples", labelWidth=labelWidth, orientation=0, callback=self.reliefChanged, callbackOnReturn = True)
        OWGUI.separator(box)

        OWGUI.comboBox(box, self, "sortBy", label = "Sort by"+"  ", items = ["No Sorting", "Attribute Name", "Number of Values"] + self.measures, orientation=0, valueType = int, callback=self.sortingChanged)


        box = OWGUI.widgetBox(self.controlArea, "Discretization", addSpace=True)
        OWGUI.spin(box, self, "nIntervals", 2, 20, label="Intervals", labelWidth=labelWidth, orientation=0, callback=self.discretizationChanged, callbackOnReturn = True)

        box = OWGUI.widgetBox(self.controlArea, "Precision", addSpace=True)
        OWGUI.spin(box, self, "nDecimals", 1, 6, label="No. of decimals", labelWidth=labelWidth, orientation=0, callback=self.decimalsChanged)

        box = OWGUI.widgetBox(self.controlArea, "Distributions", addSpace=True)
        self.cbShowDistributions = OWGUI.checkBox(box, self, "showDistributions", 'Visualize values', callback = self.cbShowDistributions)
        colBox = OWGUI.indentedBox(box, orientation = "horizontal")
        OWGUI.widgetLabel(colBox, "Color: ")
        self.colButton = OWGUI.toolButton(colBox, self, self.changeColor, width=20, height=20, debuggingEnabled = 0)
        OWGUI.rubber(colBox)

        OWGUI.rubber(self.controlArea)
        selMethBox = OWGUI.radioButtonsInBox(self.controlArea, self, "selectMethod", ["None", "All", "Manual", "Best ranked"], box="Select attributes", callback=self.selectMethodChanged)
        self.selNAttr = OWGUI.spin(OWGUI.indentedBox(selMethBox), self, "nSelected", 1, 100, label="No. selected"+"  ", orientation=0, callback=self.nSelectedChanged)

        OWGUI.separator(selMethBox)

        applyButton = OWGUI.button(selMethBox, self, "Commit", callback = self.apply)
        autoApplyCB = OWGUI.checkBox(selMethBox, self, "autoApply", "Commit automatically")

        OWGUI.setStopper(self, applyButton, autoApplyCB, "dataChanged", self.apply)

        self.table = QTableWidget()
        self.mainArea.layout().addWidget(self.table)

        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.MultiSelection)
        self.table.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
        self.table.setItemDelegate(RankItemDelegate(self, self.table))

        self.topheader = self.table.horizontalHeader()
        self.topheader.setSortIndicatorShown(1)
        self.topheader.setHighlightSections(0)

        self.resetInternals()
        self.apply()

        self.connect(self.table.horizontalHeader(), SIGNAL("sectionClicked(int)"), self.headerClick)
        self.connect(self.table, SIGNAL("clicked (const QModelIndex&)"), self.selectItem)
        self.resize(690,500)
        self.updateColor()


    def cbShowDistributions(self):
        self.table.reset()

    def changeColor(self):
        color = QColorDialog.getColor(self.distColor, self)
        if color.isValid():
            self.distColorRgb = color.getRgb()
            self.updateColor()

    def updateColor(self):
        self.distColor = QColor(*self.distColorRgb)
        w = self.colButton.width()-8
        h = self.colButton.height()-8
        pixmap = QPixmap(w, h)
        painter = QPainter()
        painter.begin(pixmap)
        painter.fillRect(0,0,w,h, QBrush(self.distColor))
        painter.end()
        self.colButton.setIcon(QIcon(pixmap))


    def resetInternals(self):
        self.data = None
        self.discretizedData = None
        self.attributeOrder = []
        self.selected = []
        self.measured = {}
        self.usefulAttributes = []
        self.dataChanged = False
        self.adjustCol0 = True
        self.lastSentAttrs = None

        self.table.setColumnCount(0)


    def selectMethodChanged(self):
        if self.selectMethod == 0:
            self.selected = []
            self.reselect()
        elif self.selectMethod == 1:
            self.selected = self.attributeOrder[:]
            self.reselect()
        elif self.selectMethod == 3:
            self.selected = self.attributeOrder[:self.nSelected]
            self.reselect()
        self.applyIf()


    def nSelectedChanged(self):
        self.selectMethod = 3
        self.selectMethodChanged()


    def sendSelected(self):
        attrs = self.data and [attr for i, attr in enumerate(self.attributeOrder) if self.table.isRowSelected(i)]
        if not attrs:
            self.send("ExampleTable Attributes", None)
            return

        self.send("ExampleTable Attributes", dataUtilities.DataTable(orange.Domain(attrs, self.data.domain.classVar), self.data))


    def setData(self,data):
        self.resetInternals()
        ##scPA       Why they only alowed data with Discrete class?  
        #self.data = self.isDataWithClass(data, orange.VarTypes.Discrete) and data or None
        self.data = data or None
        ##ecPA
        if self.data:
            ##scPA
            self.selNAttr.control.setMaximum(len(self.data.domain.attributes))
            self.selNAttr.setToolTip("Available attributes: " + str(len(self.data.domain.attributes)))
            #QToolTip.add(self.selNAttr, "Available attributes: " + str(len(self.data.domain.attributes)))
            ##ecPA
            self.adjustCol0 = True
            self.usefulAttributes = filter(lambda x:x.varType in [orange.VarTypes.Discrete, orange.VarTypes.Continuous], self.data.domain.attributes)
            self.table.setRowCount(len(self.data.domain.attributes))
            #scPA    Show the results in the GUI for the case that AutoCommit is disabled
            self.setMeasures()
            #ecPA
            self.reprint()
        ##scPA
        else:
            self.selNAttr.control.setMaximum(100)
            self.selNAttr.setToolTip("")
            #QToolTip.add(self.selNAttr, "")
        ##ecPA

        self.resendAttributes()
        self.applyIf()

    def discretizationChanged(self):
        self.discretizedData = None

        removed = False
        for meas, cont in zip(self.measuresAttrs, self.handlesContinuous):
            if not cont and self.measured.has_key(meas):
                del self.measured[meas]
                removed = True

        if self.data and self.data.domain.hasContinuousAttributes(False):
            sortedByThis = self.sortBy>=3 and not self.handlesContinuous[self.sortBy-3]
            if removed or sortedByThis:
                self.reprint()
                self.resendAttributes()
                if sortedByThis and self.selectMethod == 3:
                    self.applyIf()


    def reliefChanged(self):
        removed = False
        if self.measured.has_key("computeReliefF"):
            del self.measured["computeReliefF"]
            removed = True

        if self.data:
            sortedByReliefF = self.sortBy-3 == self.measuresAttrs.index("computeReliefF")
            if removed or sortedByReliefF:
                self.reprint()
                self.resendAttributes()
                if sortedByReliefF and self.selectMethod == 3:
                    self.applyIf()

    def selectItem(self, index):
        row = index.row()
        attr = self.attributeOrder[row]
        if attr in self.selected:
            self.selected.remove(attr)
        else:
            self.selected.append(attr)
        self.applyIf()

    def headerClick(self, index):
        if index < 0: return
        if index < 2:
            self.sortBy = 1 + index
        else:
            self.sortBy = 3 + self.measuresShort.index(str(self.table.horizontalHeader().model().headerData(index, Qt.Horizontal).toString()))
        self.sortingChanged()



    def sortingChanged(self):
        self.reprint()
        self.resendAttributes()
        if self.selectMethod == 3:
            self.applyIf()


    def setMeasures(self):
        self.selectedMeasures = [i for i, ma in enumerate(self.measuresAttrs) if getattr(self, ma)]
        self.table.setColumnCount(2 + len(self.selectedMeasures))
        self.table.setHorizontalHeaderLabels(["#", "Attributes"]+[self.measuresShort[meas_idx] for meas_idx in self.selectedMeasures])


    def measuresChanged(self):
        self.setMeasures()
        if self.data:
            self.reprint(True)
            self.resendAttributes()


    def sortByColumn(self, col):
        if col < 2:
            self.sortBy = 1 + col
        else:
            self.sortBy = 3 + self.selectedMeasures[col-2]
        self.sortingChanged()


    def decimalsChanged(self):
        self.reprint(True)


    def getMeasure(self, meas_idx):
        measAttr = self.measuresAttrs[meas_idx]
        mdict = self.measured.get(measAttr, False)
        if mdict:
            return mdict
        estimator = self.estimators[meas_idx]()
        if measAttr == "computeReliefF":
            estimator.k, estimator.m = self.reliefK, self.reliefN
        elif measAttr == "computeRF":
            #Create the RF Learner and train it
            RF = AZorngRF.RFLearner(trainingData = self.data, getVarVariance = "true")
            estimator = RF.getVarImportance 

        handlesContinuous = self.handlesContinuous[meas_idx]
        mdict = {}
        for attr in self.data.domain.attributes:
            if handlesContinuous or attr.varType == orange.VarTypes.Discrete:
                act_attr, act_data = attr, self.data
            else:
                if not self.discretizedData:
                    discretizer = orange.EquiNDiscretization(numberOfIntervals=self.nIntervals)
                    contAttrs = filter(lambda attr: attr.varType == orange.VarTypes.Continuous, self.data.domain.attributes)
                    at = []
                    attrDict = {}
                    for attri in contAttrs:
                        try:
                            nattr = discretizer(attri, self.data)
                            at.append(nattr)
                            attrDict[attri] = nattr
                        except:
                            pass
                    self.discretizedData = self.data.select(orange.Domain(at, self.data.domain.classVar))
                    self.discretizedData.setattr("attrDict", attrDict)

                act_attr, act_data = self.discretizedData.attrDict.get(attr, None), self.discretizedData

            try:
                if act_attr:
                    mdict[attr] = act_attr and estimator(act_attr, act_data)
                else:
                    mdict[attr] = None
            except:
                mdict[attr] = None

        self.measured[measAttr] = mdict
        return mdict


    def reprint(self, noSort = False):
        if not self.data:
            return

        prec = " %%.%df" % self.nDecimals

        if not noSort:
            self.resort()

        for row, attr in enumerate(self.attributeOrder):
            OWGUI.tableItem(self.table, row, 0, attr.name)
            OWGUI.tableItem(self.table, row, 1, attr.varType==orange.VarTypes.Continuous and "C" or str(len(attr.values)))

        self.minmax = {}

        for col, meas_idx in enumerate(self.selectedMeasures):
            mdict = self.getMeasure(meas_idx)
            values = filter(lambda val: val != None, mdict.values())
            if values != []:
                self.minmax[col] = (min(values), max(values))
            else:
                self.minmax[col] = (0,1)
            for row, attr in enumerate(self.attributeOrder):
                if mdict[attr] is None:
                    mattr = "NA"
                elif isinstance(mdict[attr], (int, float)):
                    mattr = prec % mdict[attr]
                else:
                    mattr = mdict[attr]
                OWGUI.tableItem(self.table, row, col+2, mattr)
        self.reselect()

        if self.sortBy < 3:
            self.topheader.setSortIndicator(self.sortBy-1, Qt.DescendingOrder)
        elif self.sortBy-3 in self.selectedMeasures:
            self.topheader.setSortIndicator(2 + self.selectedMeasures.index(self.sortBy-3), Qt.DescendingOrder)
        else:
            self.topheader.setSortIndicator(-1, Qt.DescendingOrder)

        #self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        self.table.setColumnWidth(0, 100)
        self.table.setColumnWidth(1, 20)
        for col in range(len(self.selectedMeasures)):
            self.table.setColumnWidth(col+2, 80)



    def resendAttributes(self):
        if not self.data:
            self.send("ExampleTable Attributes", None)
            return

        attrDomain = orange.Domain(  [orange.StringVariable("attributes"), orange.EnumVariable("D/C", values = "DC"), orange.FloatVariable("#")]
                                   + [orange.FloatVariable(self.measuresShort[meas_idx]) for meas_idx in self.selectedMeasures],
                                     None)
        attrData = dataUtilities.DataTable(attrDomain)
        measDicts = [self.measured[self.measuresAttrs[meas_idx]] for meas_idx in self.selectedMeasures]
        for attr in self.attributeOrder:
            cont = attr.varType == orange.VarTypes.Continuous
            attrData.append([attr.name, cont, cont and "?" or len(attr.values)] + [meas[attr] or "?" for meas in measDicts])

        self.send("ExampleTable Attributes", attrData)


    def selectRow(self, i, *foo):
        if i < 0:
            return

        attr = self.attributeOrder[i]
        if attr in self.selected:
            self.selected.remove(attr)
        else:
            self.selected.append(attr)
        self.reselect()
        self.applyIf()


    def reselect(self):
        self.table.clearSelection()

        if not self.data:
            return

        for attr in self.selected:
            self.table.selectRow(self.attributeOrder.index(attr))

        if self.selectMethod == 2 or self.selectMethod == 3 and self.selected == self.attributeOrder[:self.nSelected]:
            pass
        elif self.selected == self.attributeOrder:
            self.selectMethod = 1
        else:
            self.selectMethod = 2


    def resort(self):
        self.attributeOrder = self.usefulAttributes

        if self.sortBy:
            if self.sortBy == 1:
                st = [(attr, attr.name) for attr in self.attributeOrder]
                st.sort(lambda x,y: cmp(x[1], y[1]))
            elif self.sortBy == 2:
                st = [(attr, attr.varType == orange.VarTypes.Continuous and 1e30 or len(attr.values)) for attr in self.attributeOrder]
                st.sort(lambda x,y: cmp(x[1], y[1]))
                self.topheader.setSortIndicator(1, Qt.DescendingOrder)
            else:
                st = [(m, a == None and -1e20 or a) for m, a in self.getMeasure(self.sortBy-3).items()]
                st.sort(lambda x,y: -cmp(x[1], y[1]) or cmp(x[0], y[0]))

            self.attributeOrder = [attr for attr, meas in st]

        if self.selectMethod == 3:
            self.selected = self.attributeOrder[:self.nSelected]


    def applyIf(self):
        if self.autoApply:
            self.apply()
        else:
            self.dataChanged = True

    def apply(self):
        #self.measuresChanged()
        if not self.data or not self.selected:
            self.send("Reduced Example Table", None)
            self.lastSentAttrs = []
        else:
            if self.lastSentAttrs != self.selected:
                self.send("Reduced Example Table", dataUtilities.DataTable(orange.Domain(self.selected, self.data.domain.classVar), self.data))
                self.lastSentAttrs = self.selected[:]

        self.dataChanged = False

class RankItemDelegate(QItemDelegate):
    def __init__(self, widget = None, table = None):
        QItemDelegate.__init__(self, widget)
        self.table = table
        self.widget = widget

    def paint(self, painter, option, index):
        if not self.widget.showDistributions:
            QItemDelegate.paint(self, painter, option, index)
            return

        col = index.column()
        row = index.row()

        if col < 2 or not self.widget.minmax.has_key(col-2):        # we don't paint first two columns
            QItemDelegate.paint(self, painter, option, index)
            return

        min, max = self.widget.minmax[col-2]

        painter.save()
        self.drawBackground(painter, option, index)
        value, ok = index.data(Qt.DisplayRole).toDouble()

        if ok:        # in case we get "?" it is not ok
            smallerWidth = option.rect.width() * (max - value) / (max-min or 1)
            painter.fillRect(option.rect.adjusted(0,0,-smallerWidth,0), self.widget.distColor)

        self.drawDisplay(painter, option, option.rect, index.data(Qt.DisplayRole).toString())
        painter.restore()



if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWRank()
    a.setMainWidget(ow)
    ow.setData(dataUtilities.DataTable("../../doc/datasets/iris.tab"))
    ow.show()
    a.exec_loop()

    ow.saveSettings()

