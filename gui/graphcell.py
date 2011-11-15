#!/usr/bin/env python -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:uillaume Gay<elagachado AT  gmail DOT com>
## Commentary:

from PyQt4 import QtCore, QtGui

import math
import sys

from kt_simul import simul_spindle as Sim
from param_seter import *

from kt_simul.eval_simul import metaph_kineto_dist 
from kt_simul.xml_handler import ParamTree



paramfile = Sim.paramfile
measurefile = Sim.measurefile

SITE_OFFSET = 0.2 #Vertical distance between attachment sites
CH_OFFSET = 0.4 # Vertical distance between chromosomes

SPB_COLOR = QtGui.QColor(255,0,0)
CH_COLOR = QtGui.QColor(0,100,100, alpha = 200)
ACTIVE_SAC_COLOR = QtGui.QColor(100,10,100, alpha = 200)

GOOD_PLUGSITE_COLOR = QtGui.QColor(0,250,20, alpha = 200)
BAD_PLUGSITE_COLOR = QtGui.QColor(255,0,0, alpha = 255)
UNPLUGED_COLOR = QtGui.QColor(0,20,250, alpha = 200)
        


class GraphCell(QtGui.QGraphicsItem):
    '''
    This is the parent item containing all the objects within the cell
    '''
    
    def __init__(self, mt, parent = None):

        QtGui.QGraphicsItem.__init__(self, parent = parent)

        N = mt.KD.params['N']
        Mk = mt.KD.params['Mk']
        self.mt = mt # A SimMetaphase instance
        self.items = []
        self.spbR = GraphSPB(0, parent = self)
        self.items.append(self.spbR)
        self.spbL = GraphSPB(-1, parent = self)
        self.items.append(self.spbL)
       
        for n in range(N):
            left_kineto = GraphKinetochore(n, -1, parent = self)
            self.items.append(left_kineto)
            right_kineto = GraphKinetochore(n, 0, parent = self)
            self.items.append(right_kineto)            
        
        self.time_point = 0
        for item in self.items:
            item.setPos(item.getSimPos(0))

    def advance(self):

        self.time_point += 1
        return self.gotoTime(self.time_point)

    def gotoTime(self, time_point):
        
        self.time_point = time_point
        for item in self.items:
            newPos = item.getSimPos(self.time_point)
            if newPos == item.pos():
                return False
            item.setPos(newPos)
            if isinstance(item, GraphKinetochore):
                for plugsite in item.plugsites:
                    newPos = plugsite.getSimPos(self.time_point)
                    plugsite.setPos(newPos)
                    if plugsite.newPlug != plugsite.sim.state_hist[self.time_point]:
                        plugsite.newPlug = plugsite.sim.state_hist[self.time_point]
                        if plugsite.newPlug == 1:
                            brush = QtGui.QBrush(GOOD_PLUGSITE_COLOR)
                        elif plugsite.newPlug == -1:
                            brush = QtGui.QBrush(BAD_PLUGSITE_COLOR)
                        else: 
                            brush = QtGui.QBrush(UNPLUGED_COLOR)
                        plugsite.color = brush
        return True


    def boundingRect(self):
        '''
        As this is not easily changed (or not supposed to, we give a
        large box
        '''
        N = self.mt.KD.params['N']
        Mk = self.mt.KD.params['Mk']
        
        height = 5.#N * ( Mk * SITE_OFFSET + CH_OFFSET)
        width = 14. # microns, shall be enough
        
        return QtCore.QRectF(-width/2, -height/2, width, height)

    def paint(self, painter, option, widget):

        painter.setBrush(QtCore.Qt.white)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.1))
        painter.drawRoundedRect(self.boundingRect(), 30, 100, QtCore.Qt.RelativeSize)


class GraphKinetochore(QtGui.QGraphicsItem):

    def __init__(self, n, side, parent = None):

        QtGui.QGraphicsItem.__init__(self, parent = parent)
        self.graphcell = parent

        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)

        N = self.graphcell.mt.KD.params['N']
        Mk = self.graphcell.mt.KD.params['Mk']

        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
        self.n = n
        self.side = side
        self.ch = self.graphcell.mt.KD.chromosomes[n]

        if side == 0:
            x = self.ch.rightpos
            self.traj = self.ch.righttraj
        else:
            x = self.ch.leftpos
            self.traj = self.ch.lefttraj
        self.y = ( n - (N - 1) / 2. ) * 0.4 #* ( Mk * SITE_OFFSET + CH_OFFSET)

        self.width = 0.2 
        self.height = 0.1* Mk #Mk * SITE_OFFSET + CH_OFFSET
        self.newPos = QtCore.QPointF(x, self.y)

        self.plugsites = []
        self.setZValue(n)
        for m in range(Mk):
            self.plugsites.append(GraphPlugSite(m, self, parent))
        
        
    def getSimPos(self, time_point):

        try:
            x = self.traj[time_point]
        except IndexError:
            x = self.traj[-1]
        return  QtCore.QPointF(x, self.y)

    def advance(self):

        self.time_point += 1
        for item in self.items:
            newPos = item.getSimPos(self.time_point)
            if newPos == item.pos():
                return False
            item.setPos(newPos)

        return True

    def shape(self):
        
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(self.width, self.height,
                             2*self.width, 2*self.height)
        return self.path
    
    def paint(self, painter, option, widget):
        if self.ch.active_sac:
            brush = QtGui.QBrush(ACTIVE_SAC_COLOR)
        else:
            brush = QtGui.QBrush(CH_COLOR)
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.newPos
        painter.drawEllipse(center, self.width, self.height)

    def boundingRect(self):
        adjust = 1.
        return QtCore.QRectF(-self.width - adjust, -self.height - adjust,
                             2*self.width  + adjust, 2* self.height + adjust)

class GraphPlugSite(QtGui.QGraphicsItem):

    def __init__(self, m, kineto, parent = None):

        QtGui.QGraphicsItem.__init__(self, parent = parent)
        self.kineto = kineto
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)

        Mk = self.kineto.graphcell.mt.KD.params['Mk']

        self.m = m
        side = self.kineto.side
        if side == 0:
            self.sim = self.kineto.ch.rplugs[m]
        else:    
            self.sim = self.kineto.ch.lplugs[m]
        
        self.width = 0.1
        self.height = self.kineto.height/Mk 

        self.y = self.kineto.y + ( m - (Mk -1 )/2.) * 0.1
        x = self.sim.pos 
        self.newPos = QtCore.QPointF(x, self.y)
        self.newPlug = self.sim.plug
        self.setZValue(self.kineto.n * (1 + m))

        if self.sim.plug == 1:
            brush = QtGui.QBrush(GOOD_PLUGSITE_COLOR)
        elif self.sim.plug == -1:
            brush = QtGui.QBrush(BAD_PLUGSITE_COLOR)
        else: 
            brush = QtGui.QBrush(UNPLUGED_COLOR)
        self.color = brush


    def getSimPos(self, time_point):

        try:
            x = self.sim.traj[time_point]
        except IndexError:
            x = self.sim.traj[-1]

            
        return  QtCore.QPointF(x, self.y)

    def advance(self):

        self.time_point += 1
        for item in self.items:
            newPos = item.getSimPos(self.time_point)
            if newPos == item.pos():
                return False
            item.setPos(newPos)
        return True

    def shape(self):
        
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(-self.width/2., -self.height/2,
                             self.width, self.height)
        return self.path
        
    def paint(self, painter, option, widget):

        brush = self.color
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.newPos
        painter.drawEllipse(center, self.width, self.height)

    def boundingRect(self):
        adjust = 1.
        return QtCore.QRectF(-self.width - adjust, -self.height - adjust,
                             2*self.width  + adjust, 2* self.height + adjust)


        
class GraphSPB(QtGui.QGraphicsItem):

    def __init__(self, side, parent = None):

        QtGui.QGraphicsItem.__init__(self, parent = parent)
        self.graphcell = parent
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)

        if side == 0:
            self.sim = self.graphcell.mt.KD.spbR
        else:
            self.sim = self.graphcell.mt.KD.spbL
        self.newPos = QtCore.QPointF(self.sim.pos, 0.)

    
    def getSimPos(self, time_point):
        try:
            x = self.sim.traj[time_point]
        except IndexError:
            x = self.sim.traj[-1]

        y = 0.
        return  QtCore.QPointF(x, y)

    def advance(self):
        if self.newPos == self.pos():
            return False

        self.setPos(self.newPos)
        return True

    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(-0.1, -0.3, 0.2, 0.6)

        return self.path
        
    def paint(self, painter, option, widget):
        brush = QtGui.QBrush(SPB_COLOR)
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.newPos
        painter.drawEllipse(center, 0.2, 0.6)

    def boundingRect(self):
        adjust = 5.
        return QtCore.QRectF(-0.2 - adjust, -0.6  - adjust,
                             0.4  + adjust, 1.2 + adjust)

    
        
class NakedWidget(QtGui.QGraphicsView):

    def __init__(self, mt):

        #QtGui.QGraphicsView.__init__()
        super(NakedWidget, self).__init__()
        self.timerId = 0

        self.cell = GraphCell(mt)

        scene = QtGui.QGraphicsScene(self)

        #self.setCacheMode(QtGui.QGraphicsView.CacheBackground)
        self.setViewportUpdateMode(QtGui.QGraphicsView.BoundingRectViewportUpdate)
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

        # # Text.
        # textRect = QtCore.QRectF(sceneRect.left() + 4, sceneRect.top() + 4,
        #         sceneRect.width() - 4, sceneRect.height() - 4)
        # message = "Click and drag the nodes around, and zoom with the " \
        #         "mouse wheel or the '+' and '-' keys"

        

    def scaleView(self, scaleFactor):
        factor = self.matrix().scale(scaleFactor, scaleFactor).mapRect(
            QtCore.QRectF(0, 0, 1, 1)).width()

        if factor < 0.1 or factor > 200:
            return

        self.scale(scaleFactor, scaleFactor)


    def startAnim(self):
        if not self.timerId:
            self.timerId = self.startTimer(1000 / FRAME_RATE)
    

    def timerEvent(self, event):

        itemsMoved = False
        
        if self.cell.advance():
            itemsMoved = True

        if not itemsMoved:
            self.killTimer(self.timerId)
            self.timerId = 0



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
    widget = NakedWidget(mt)
    widget.setRenderHint(QtGui.QPainter.Antialiasing)
    widget.show()
    
    
    
    sys.exit(app.exec_())
