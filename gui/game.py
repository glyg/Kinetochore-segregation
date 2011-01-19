#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:uillaume Gay<elagachado AT  gmail DOT com>
## Commentary:

from PyQt4 import QtCore, QtGui

import math
import sys
sys.path.append('/home/guillaume/Python')

from kt_simul import simul_spindle as Sim
from param_seter import *
from game import *

from kt_simul.eval_simul import metaph_kineto_dist 
from kt_simul.xml_handler import ParamTree



paramfile = Sim.paramfile
measurefile = Sim.measurefile

SITE_OFFSET = 0.2 #Vertical distance between attachment sites
CH_OFFSET = 0.4 # Vertical distance between chromosomes

SPB_COLOR = QtGui.QColor(0,0,255)
CH_COLOR = QtGui.QColor(0.2,0.2,0.2, alpha = 200)
GOOD_PLUGSITE_COLOR = QtGui.QColor(0,1.,0, alpha = 200)
BAD_PLUGSITE_COLOR = QtGui.QColor(1.,0,0, alpha = 255)

        
        


class GraphCell(QtGui.QGraphicsItem):
    '''
    This is the parent item containing all the objects within the cell
    '''
    
    def __init__(self, mt, parent = None):

        QtGui.QGraphicsItem.__init__(self, parent = parent)

        self.mt = mt # A SimMetaphase instance
        self.items = []
        self.spbR = GraphSPB(0, parent = self)
        self.items.append(self.spbR)
        self.spbL = GraphSPB(-1, parent = self)
        self.items.append(self.spbL)

        for item in self.items:
            item.setPos(item.newPos)
    

    def boundingRect(self):
        '''
        As this is not easily changed (or not supposed to, we give a
        large box
        '''
        N = self.mt.KD.params['N']
        Mk = self.mt.KD.params['Mk']
        
        height = N * ( Mk * SITE_OFFSET + CH_OFFSET)
        width = 14 # microns, shall be enough
        
        return QtCore.QRectF(-width/2, -height/2, width, height)

    def paint(self, painter, option, widget):

        painter.setBrush(QtCore.Qt.white)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.1))
        painter.drawRoundedRect(self.boundingRect(), 30, 100, QtCore.Qt.RelativeSize)

class GraphChromosome(QtGui.QGraphicsItem):
    def __init__(self, n, parent = None):
        
        QtGui.QGraphicsItem.__init__(self, parent = parent)
        self.graphcell = parent
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
        
        if side == 0:
            self.sim = self.graphcell.mt.KD.chromosomes[n].pos
        else:
            self.sim = self.graphcell.mt.KD.spbL
                
        self.newPos = QtCore.QPointF(self.sim.pos, 0.)

     
class GraphSPB(QtGui.QGraphicsItem):

    def __init__(self, side, parent = None):

        QtGui.QGraphicsItem.__init__(self, parent = parent)
        self.graphcell = parent
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)

        if side == 0:
            self.sim = self.graphcell.mt.KD.spbR
        else :
            self.sim = self.graphcell.mt.KD.spbL
        self.newPos = QtCore.QPointF(self.sim.pos, 0.)

    def advance(self):
        if self.newPos == self.pos():
            return False

        self.setPos(self.newPos)
        return True

    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(-0.1, -0.3, 0.2, 0.6)
            
    def paint(self, painter, option, widget):
        brush = QtGui.QBrush(SPB_COLOR)
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.newPos
        painter.drawEllipse(center, 0.2, 0.6)

    def boundingRect(self):
        adjust = 0.1
        return QtCore.QRectF(-0.1 - adjust, -0.3  - adjust,
                             0.2  + adjust, 0.6 + adjust)

    
        
class InteractiveWidget(QtGui.QGraphicsView):

    def __init__(self, mt):

        #QtGui.QGraphicsView.__init__()
        super(InteractiveWidget, self).__init__()
        self.timerId = 0

        cell = GraphCell(mt)

        scene = QtGui.QGraphicsScene(self)
        scene.setSceneRect(-9, -4, 18, 8)
        self.setScene(scene)
        
        scene.addItem(cell)
        
        self.scale(30, 30)



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

        # # Text.
        # textRect = QtCore.QRectF(sceneRect.left() + 4, sceneRect.top() + 4,
        #         sceneRect.width() - 4, sceneRect.height() - 4)
        # message = "Click and drag the nodes around, and zoom with the " \
        #         "mouse wheel or the '+' and '-' keys"

        

    def scaleView(self, scaleFactor):
        factor = self.matrix().scale(scaleFactor, scaleFactor).mapRect(
            QtCore.QRectF(0, 0, 1, 1)).width()

        if factor < 0.07 or factor > 100:
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
    widget = InteractiveWidget(mt)
    widget.setRenderHint(QtGui.QPainter.Antialiasing)
    widget.show()
    
    sys.exit(app.exec_())












# class PlayGround(QtGui.QGraphicsScene):

#     def __init__(self, parent = None):
#         QtGui.QGraphiscScene.__init__(-300, -100, 600,
#                                       200, parent = parent)
        
