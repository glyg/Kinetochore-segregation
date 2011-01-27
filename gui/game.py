#!/usr/bin/env python -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:uillaume Gay<elagachado AT  gmail DOT com>
## Commentary:

from PyQt4 import QtCore, QtGui

import math
import sys, os
sys.path.append('/home/guillaume/Python')

from kt_simul import simul_spindle as Sim
from param_seter import *
from game import *

from kt_simul.eval_simul import metaph_kineto_dist 
from kt_simul.xml_handler import ParamTree
from graphcell import *

FRAME_RATE = 25 #image/seconds        
        
class InteractiveCellWidget(QtGui.QGraphicsView):

    def __init__(self, mt):
        
        super(InteractiveCellWidget, self).__init__()
        vbox = QtGui.QVBoxLayout()

        self.timerId = 0

        self.view = ViewCellWidget(mt)
        self.ccw = ControlCellWidget(mt)
        self.connect(self.ccw.playButton, QtCore.SIGNAL('clicked()'),
                     self.play)
        self.connect(self.ccw.pauseButton, QtCore.SIGNAL('clicked()'),
                     self.pause)
        # self.connect(self.ccw.slider, QtCore.SIGNAL('sliderMoved()'),
        #              self.gotoTime)
        self.ccw.slider.valueChanged.connect(self.gotoTime)

        vbox.setMargin(5)
        vbox.addWidget(self.view)
        vbox.addWidget(self.ccw)
        self.setLayout(vbox)
        

    def pause(self):

        if not self.timerId:
            pass
        else:
            self.killTimer(self.timerId)
            self.timerId = 0

    def gotoTime(self, time):
        if self.ccw.slider.hasFocus():
            self.pause()
            rects = [self.view.cell.boundingRect()]
            self.view.updateScene(rects)
            self.view.cell.gotoTime(time)
            self.timerId = 0

    def play(self):
        ### I'll re-implement this when user can move things in graphcell
        if not self.timerId:
            self.timerId = self.startTimer(1000 / FRAME_RATE)
        

    def startAnim(self):
        if not self.timerId:
            self.timerId = self.startTimer(1000 / FRAME_RATE)
    

    def timerEvent(self, event):

        itemsMoved = False
        rects = [self.view.cell.boundingRect()]
        self.view.updateScene(rects)
        if self.view.cell.advance():
            itemsMoved = True
            self.ccw.slider.setValue(self.view.cell.time_point)

        if not itemsMoved:
            self.killTimer(self.timerId)
            self.timerId = 0
            





class ControlCellWidget(QtGui.QWidget):    
    
    def __init__(self, mt, parent = None):
        
        super(ControlCellWidget, self).__init__(parent = parent)
        numsteps = int(mt.KD.params['span'] / mt.KD.params['dt'])

        orientation = QtCore.Qt.Horizontal
        self.slider = QtGui.QSlider(orientation)
        self.slider.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider.setTickPosition(QtGui.QSlider.TicksBothSides)
        self.slider.setTickInterval(10)
        self.slider.setSingleStep(1)
    
        self.slider.setMinimum(0)
        self.slider.setMaximum(numsteps)
        self.slider.setValue(0)

        self.playButton = QtGui.QPushButton(self)
        play_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "player_play.svg")
        self.playButton.setIcon(QtGui.QIcon(play_icon))

        self.pauseButton = QtGui.QPushButton(self)
        pause_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "player_pause.svg")
        self.pauseButton.setIcon(QtGui.QIcon(pause_icon))

        hbox = QtGui.QHBoxLayout()
        hbox.setMargin(5)
        hbox.addWidget(self.playButton)
        hbox.addWidget(self.pauseButton)
        hbox.addWidget(self.slider)
        self.setLayout(hbox)
        

class ViewCellWidget(QtGui.QGraphicsView):

    def __init__(self, mt):

        super(ViewCellWidget, self).__init__()
        self.timerId = 0
        self.cell = GraphCell(mt)
        scene = QtGui.QGraphicsScene(self)

        self.setCacheMode(QtGui.QGraphicsView.CacheBackground)
        #self.setViewportUpdateMode(QtGui.QGraphicsView.BoundingRectViewportUpdate)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.setTransformationAnchor(QtGui.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtGui.QGraphicsView.AnchorViewCenter)

        scene.setSceneRect(-9, -4, 18, 8)
        self.setScene(scene)
        
        scene.addItem(self.cell)
        
        self.scale(40, 40)


        

    def wheelEvent(self, event):
        self.scaleView(math.pow(2.0, -event.delta() / 240.0))

    def drawBackground(self, painter, rect):
        # Shadow.
        sceneRect = self.sceneRect()
        rightShadow = QtCore.QRectF(sceneRect.right(), sceneRect.top() + 0.1, 0.1,
                sceneRect.height())
        bottomShadow = QtCore.QRectF(sceneRect.left() + 0.1, sceneRect.bottom(),
                sceneRect.width(), 0.1)
        if rightShadow.intersects(rect) or rightShadow.contains(rect):
	        painter.fillRect(rightShadow, QtCore.Qt.darkGray)
        if bottomShadow.intersects(rect) or bottomShadow.contains(rect):
	        painter.fillRect(bottomShadow, QtCore.Qt.darkGray)

        # Fill.
        gradient = QtGui.QLinearGradient(sceneRect.topLeft(),
                sceneRect.bottomRight())
        gradient.setColorAt(0, QtCore.Qt.white)
        gradient.setColorAt(1, QtCore.Qt.lightGray)
        painter.fillRect(rect.intersect(sceneRect), QtGui.QBrush(gradient))
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawRect(sceneRect)


    def scaleView(self, scaleFactor):
        factor = self.matrix().scale(scaleFactor, scaleFactor).mapRect(
            QtCore.QRectF(0, 0, 1, 1)).width()

        if factor < 0.1 or factor > 200:
            return

        self.scale(scaleFactor, scaleFactor)




if __name__ == '__main__':

    from kineto_simulation import SigMetaphase

    app = QtGui.QApplication(sys.argv)
    QtCore.qsrand(QtCore.QTime(0,0,0).secsTo(QtCore.QTime.currentTime()))

    paramtree = ParamTree(paramfile)
    paramtree.create_dic(adimentionalized = False)
    
    measuretree = ParamTree(measurefile)
    measuretree.create_dic(adimentionalized = False)
    measures = measuretree.dic


    mt = SigMetaphase(paramtree, measuretree)
    widget = InteractiveCellWidget(mt)
    widget.setRenderHint(QtGui.QPainter.Antialiasing)
    widget.show()
    
    
    
    sys.exit(app.exec_())












# class PlayGround(QtGui.QGraphicsScene):

#     def __init__(self, parent = None):
#         QtGui.QGraphiscScene.__init__(-300, -100, 600,
#                                       200, parent = parent)
        
