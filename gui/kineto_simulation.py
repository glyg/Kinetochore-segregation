#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Graphical User Interface for the kinetochore dynamics simulation

'''

import sys, os, random
from PyQt4 import QtCore, QtGui
from numpy import arange, sin, pi, array, linalg, mod
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

import pyximport
pyximport.install()


sys.path.append('/home/guillaume/Python/')


from kt_simul import simul_spindle as Sim
from param_seter import *
from game import *

from kt_simul.eval_simul import metaph_kineto_dist 
from kt_simul.xml_handler import ParamTree



paramfile = Sim.paramfile
measurefile = Sim.measurefile


__all__ = ['MainWindow', paramfile, measurefile]

#matplotlib.use('Qt4Agg') #Useless

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, span = 800., parent=None, width=5, height=4, dpi=100):
        self.fig = Figure( figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure(span)

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

      

    def compute_initial_figure(self, span):
        self.axes.axis([0, span , -10, 10 ])

        self.axes.set_xlabel('Time (seconds)', fontsize = 12)
        self.axes.set_ylabel(u'Distance from center (µm)', fontsize = 12)
        
    def update_figure(self, mt, n = None):

        ''' Plot the different trajectories
        '''
        
        self.axes.clear()
        self.axes.hold(True)
        if n == None:
            mt.show_trajs(self.axes)
        else:
            mt.show_one(n = n, fig = self.fig)

        self.draw()


class SigMetaphase(Sim.Metaphase, QtGui.QWidget):

    '''
    See if we can retrieve signals from this hybrid

    overrides _one_step method of the Metaphase class

    '''
    
    def __init__(self, paramtree, measuretree, parent = None):
        paramtree.create_dic(adimentionalized = True)
        duration = paramtree.dic["span"]
        Sim.Metaphase.__init__(self, duration, paramtree, measuretree)
        QtGui.QWidget.__init__(self, None)
        self.date = 0
        
    def _one_step(self):
            
        if not self.KD.anaphase:
            self.KD.plug_unplug()
            self.emit(QtCore.SIGNAL('inMetaphase'))
            
        A = self.KD.calcA()
        b = - self.KD.calcb()
        speeds = linalg.solve(A,b)
        self.KD.position_update(speeds)

        nb_mero = self._mero_checkpoint()
        if nb_mero > 0:
            self.emit(QtCore.SIGNAL('meroCheckPoint'), nb_mero)

        self.emit(QtCore.SIGNAL('plugCheckPoint'), self._plug_checkpoint())
        self.date += 1
        self.emit(QtCore.SIGNAL('stepDone'), self.date)
        self.emit(QtCore.SIGNAL('stepDone_nop'))
        
    def sig_simul(self):
        
        self.simul()
        self.emit(QtCore.SIGNAL('simulDone'),self.report)
        

class MainWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):

        self.paramtree = ParamTree(paramfile)
        self.paramtree.create_dic(adimentionalized = False)

        self.measuretree = ParamTree(measurefile)
        self.measuretree.create_dic(adimentionalized = False)
        self.measures = self.measuretree.dic

        QtGui.QMainWindow.__init__(self, parent)
        w = QtGui.QWidget()
        self.setCentralWidget(w)

        #Parameter Setting in a Dock Widget
        self.dock = QtGui.QDockWidget('Parameters Setting')
        self.setParameters = SetParameters(self.paramtree)
        scrollArea = QtGui.QScrollArea()
        scrollArea.setWidget(self.setParameters)
        self.dock.setWidget(scrollArea)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.dock);
        #

        #Measures Setting in another Dock Widget
        self.dock = QtGui.QDockWidget('Measures Setting')
        self.setMeasures = SetMeasures(self.measuretree)
        scrollArea = QtGui.QScrollArea()
        scrollArea.setWidget(self.setMeasures)
        self.dock.setWidget(scrollArea)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.dock);
        #

        #Text Area
        s = ("Welcome to the S.Pombe kinetochore motion "
             "Simulator")
        self.simLog = QtGui.QTextEdit(s)
        self.simLog.setReadOnly(True)

        #Plotting Areas
        span = self.paramtree.dic['span']
        self.plotarea1 = MyMplCanvas(span)
        mpl_toolbar = NavigationToolbar(self.plotarea1, w)
        vbox1 = QtGui.QVBoxLayout()
        vbox1.setMargin(5)
        vbox1.addWidget(self.plotarea1)
        vbox1.addWidget(mpl_toolbar)
        self.w1 = QtGui.QWidget()
        self.w1.setLayout(vbox1)

        
        self.plotarea2 = MyMplCanvas(span)
        mpl_toolbar = NavigationToolbar(self.plotarea2, w)
        vbox2 = QtGui.QVBoxLayout()
        vbox2.setMargin(5)
        vbox2.addWidget(self.plotarea2)
        vbox2.addWidget(mpl_toolbar)
        self.w2 = QtGui.QWidget()
        self.w2.setLayout(vbox2)


        #All this goes in a tab widget
        self.tabWidget = QtGui.QTabWidget()
        self.tabWidget.addTab(self.simLog, "Log")
        #self.tabWidget.addTab(w1, "All trajectories")
        #self.tabWidget.addTab(w2, "One trajectory")
        #Buttons
        #self.buttonGroup = QtGui.QButtonGroup()
        runButton = QtGui.QPushButton('Run the simulation')
        exec_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "exec.svg")
        runButton.setIcon(QtGui.QIcon(exec_icon))
        self.connect(runButton, QtCore.SIGNAL('clicked()'), self.run_simulation)

        showTrajButton = QtGui.QPushButton('Show trajectories')
        self.connect(showTrajButton, QtCore.SIGNAL('clicked()'), self.show_trajs)
       
        self.traj_num = 0
        showOneButton = QtGui.QPushButton('Show one trajectory')
        self.connect(showOneButton, QtCore.SIGNAL('clicked()'), self.show_one)
        
        self.interactiveButton = QtGui.QRadioButton('Interactive Simulation')
        self.interactiveButton.setChecked(True)
        #self.buttonGroup.addButton(runButton)

        #Progress Bar
        self.progressBar = QtGui.QProgressBar()

        hbox = QtGui.QHBoxLayout()
        hbox.setMargin(5)

        hbox.addWidget(runButton)
        hbox.addWidget(showTrajButton)
        hbox.addWidget(showOneButton)
        hbox.addWidget(self.progressBar)
        hbox.addWidget(self.interactiveButton)

        vbox = QtGui.QVBoxLayout()
        vbox.setMargin(5)
        vbox.addLayout(hbox)#self.buttonGroup)
        vbox.addWidget(self.tabWidget)

        w.setLayout(vbox)

        self.createActions()
        self.createMenus()
        
        self.createToolBars()
        self.createStatusBar()
        self.setCurrentFile(paramfile)

        self.setWindowTitle(self.tr("Kinetochore Dynamics Simulation"))
        self.setMinimumSize(160,160)
        self.resize(1000,600)


    def run_simulation(self):

        QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.mt = SigMetaphase(self.paramtree, self.measuretree)
        self.progressBar.setMaximum(int(self.paramtree.dic['span']/self.paramtree.dic['dt']))
        self.progressBar.setMinimum(0)
        self.connect(self.mt, QtCore.SIGNAL('plugCheckPoint'), self.active_checkpoint)
        self.connect(self.mt,  QtCore.SIGNAL('stepDone'), self.progressBar.setValue)
        self.connect(self.mt,  QtCore.SIGNAL('simulDone'), self.print_report)

        run_interactively = self.interactiveButton.isChecked()
        if run_interactively:
            self.iw = InteractiveCellWidget(self.mt)
            self.iw.setRenderHint(QtGui.QPainter.Antialiasing)

            idx = self.tabWidget.currentIndex()+1
            self.tabWidget.insertTab(idx, self.iw, "Interactive Simulation")
            self.tabWidget.setCurrentIndex(idx)

            self.connect(self.mt,  QtCore.SIGNAL('simulDone'),
                         self.iw.startAnim)
            # self.connect(self.mt,  QtCore.SIGNAL('stepDone_nop'),
            #              self.iw.scene().advance)
            self.mt.sig_simul()

            QtGui.QApplication.restoreOverrideCursor()
            
        else:
            self.mt.sig_simul()
            #self.plotarea.update_figure(self.mt)
            QtGui.QApplication.restoreOverrideCursor()

    def show_trajs(self):

        self.plotarea1.update_figure(self.mt)
        self.tabWidget.insertTab(1,self.w1, "All trajectories")
        self.tabWidget.setCurrentIndex(1)


    def show_one(self):

        #self.tabWidget.removeTab(2)
        t_num = mod(self.traj_num, 3)
        self.plotarea2.update_figure(self.mt, t_num)
        self.traj_num += 1
        self.tabWidget.insertTab(2, self.w2, "Trajectory of chromosome # %i" %(t_num + 1) )
        self.tabWidget.setCurrentIndex(2)
        
    def update_progressBar(self, val):
        self.progressBar.setValue(val)

    
    def print_report(self, report):
        self.simLog.append(QtCore.QString("Simulation's done!"))
        for l in report:
            ls = QtCore.QString(l)
            self.simLog.append(QtCore.QString(ls))

        self.progressBar.setValue(0)
            
    def active_checkpoint(self, cp):
        if cp:
            self.statusBar().showMessage(self.tr("Active Plug/unplug checkpoint"), 2000)
        else:
            self.statusBar().showMessage(self.tr(" "), 2000)

        
    def closeEvent(self, event):
        if self.maybeSave():
            event.accept()
        else:
            event.ignore()

    def newFile(self):
        if self.maybeSave():
            self.textEdit.clear()
            self.setCurrentFile(QtCore.QString())

    def open(self):
        if self.maybeSave():
            fileName = QtGui.QFileDialog.getOpenFileName(self)
            if not fileName.isEmpty():
                self.loadFile(fileName)

    def save(self):
        #         if self.curFile.isEmpty():
        #             return self.saveAs()
        #         else:
        return self.saveFile(self.curFile)

    def saveAs(self):
        fileName = QtGui.QFileDialog.getSaveFileName(self)
        if fileName.isEmpty():
            return False

        return self.saveFile(fileName)

    def about(self):
        QtGui.QMessageBox.about(self, self.tr("About Application"),
            self.tr("This little stuff allows to simulate "
                    "kinetochore dynamics in S.Pombe "))

    def documentWasModified(self):
        self.setWindowModified(self.setParameters.isModified)

    def createActions(self):

        open_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "open.png")
        self.openAct = QtGui.QAction(QtGui.QIcon(open_icon),
                                     self.tr("&Open Parameters File"), self)
        self.openAct.setShortcut(self.tr("Ctrl+O"))
        self.openAct.setStatusTip(self.tr("Open an existing parameters file"))
        self.connect(self.openAct, QtCore.SIGNAL("triggered()"), self.open)
        save_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "save.png")
        self.saveAct = QtGui.QAction(QtGui.QIcon(save_icon),
                                     self.tr("&Save"), self)
        self.saveAct.setShortcut(self.tr("Ctrl+S"))
        self.saveAct.setStatusTip(self.tr("Save the paramters to disk"))
        self.connect(self.saveAct, QtCore.SIGNAL("triggered()"), self.save)

        self.saveAsAct = QtGui.QAction(self.tr("Save &As..."), self)
        self.saveAsAct.setStatusTip(self.tr("Save the document under a new name"))
        self.connect(self.saveAsAct, QtCore.SIGNAL("triggered()"), self.saveAs)

        self.exitAct = QtGui.QAction(self.tr("E&xit"), self)
        self.exitAct.setShortcut(self.tr("Ctrl+Q"))
        self.exitAct.setStatusTip(self.tr("Exit the application"))
        self.connect(self.exitAct, QtCore.SIGNAL("triggered()"), self, QtCore.SLOT("close()"))

        self.aboutAct = QtGui.QAction(self.tr("&About"), self)
        self.aboutAct.setStatusTip(self.tr("Show the application's About box"))
        self.connect(self.aboutAct, QtCore.SIGNAL("triggered()"), self.about)

        self.aboutQtAct = QtGui.QAction(self.tr("About &Qt"), self)
        self.aboutQtAct.setStatusTip(self.tr("Show the Qt library's About box"))
        self.connect(self.aboutQtAct, QtCore.SIGNAL("triggered()"), QtGui.qApp, QtCore.SLOT("aboutQt()"))


    def createMenus(self):
        self.fileMenu = self.menuBar().addMenu(self.tr("&File"))
        self.fileMenu.addAction(self.openAct)
        self.fileMenu.addAction(self.saveAct)
        self.fileMenu.addAction(self.saveAsAct)
        self.fileMenu.addSeparator();
        self.fileMenu.addAction(self.exitAct)


        self.menuBar().addSeparator()

        self.helpMenu = self.menuBar().addMenu(self.tr("&Help"))
        self.helpMenu.addAction(self.aboutAct)
        self.helpMenu.addAction(self.aboutQtAct)

    def createToolBars(self):
        self.fileToolBar = self.addToolBar(self.tr("File"))
        self.fileToolBar.addAction(self.openAct)
        self.fileToolBar.addAction(self.saveAct)

    def createStatusBar(self):
        self.statusBar().showMessage(self.tr("Ready"))


    def maybeSave(self):
        if self.setParameters.isModified:
            ret = QtGui.QMessageBox.warning(self, self.tr("Application"),
                        self.tr("The document has been modified.\n"
                                "Do you want to save your changes?"),
                        QtGui.QMessageBox.Yes | QtGui.QMessageBox.Default,
                        QtGui.QMessageBox.No,
                        QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Escape)
            if ret == QtGui.QMessageBox.Yes:
                return self.save()
            elif ret == QtGui.QMessageBox.Cancel:
                return False
        return True

    def loadFile(self, fileName):
        
        QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.paramtree = ParamTree(fileName)

        self.removeDockWidget (self.dock)
        del self.dock
        self.dock = QtGui.QDockWidget('Parameters Setting')
        self.setParameters = SetParameters(self.paramtree)
        scrollArea = QtGui.QScrollArea()
        scrollArea.setWidget(self.setParameters)
        self.dock.setWidget(scrollArea)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.dock);

        QtGui.QApplication.restoreOverrideCursor()

        self.setCurrentFile(fileName)
        self.statusBar().showMessage(self.tr("Parameters loaded"), 2000)


    def saveFile(self, fileName):

        outf = open(fileName, 'w+')
        print fileName
        indent(self.paramtree.root)
        outf.write(tostring(self.paramtree.root))
        outf.close()

        self.setCurrentFile(fileName);
        self.statusBar().showMessage(self.tr("File saved"), 2000)
        return True

    def setCurrentFile(self, fileName):

        self.curFile = fileName
        self.setParameters.setModified(False)
        self.setWindowModified(False)


    def strippedName(self, fullFileName):
        return QtCore.QFileInfo(fullFileName).fileName()




if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    mainwindow = MainWindow()
    mainwindow.show()
    sys.exit(app.exec_())
    
