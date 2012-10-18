#!/usr/bin/env python -*- coding: utf-8 -*-

## Title:
## Description:
## Author:uillaume Gay<elagachado AT  gmail DOT com>
## Commentary:

from PySide import QtCore, QtGui

import math
import sys
import os

FRAME_RATE = 25  # image/seconds

SITE_OFFSET = 0.2  # Vertical distance between attachment sites
CH_OFFSET = 0.4  # Vertical distance between chromosomes

SPB_COLOR = QtGui.QColor(255, 0, 0)
CH_COLOR = QtGui.QColor(0, 100, 100, alpha=200)
ACTIVE_SAC_COLOR = QtGui.QColor(100, 10, 100, alpha=200)

GOOD_PLUGSITE_COLOR = QtGui.QColor(0, 250, 20, alpha=200)
BAD_PLUGSITE_COLOR = QtGui.QColor(255, 0, 0, alpha=255)
UNPLUGED_COLOR = QtGui.QColor(0, 20, 250, alpha=200)


class GraphCell(QtGui.QGraphicsItem):
    """
    This is the parent item containing all the objects within the cell
    """
    def __init__(self, mt, parent=None):
        QtGui.QGraphicsItem.__init__(self, parent=parent)

        self.N = int(mt.KD.params['N'])
        self.Mk = int(mt.KD.params['Mk'])
        self.mt = mt  # Metaphase instance
        self.items = []
        self.spbR = GraphSPB(0, parent=self)
        self.items.append(self.spbR)
        self.spbL = GraphSPB(-1, parent=self)
        self.items.append(self.spbL)

        for n in range(self.N):
            left_kineto = GraphCentromere(n, -1, parent=self)
            self.items.append(left_kineto)
            right_kineto = GraphCentromere(n, 0, parent=self)
            self.items.append(right_kineto)

        self.time_point = 0
        for item in self.items:
            item.setPos(item.get_pos(0))

        print self.mt.KD.spbR.traj[0]
        print self.spbR.get_pos(0)

        print self.mt.KD.spbL.traj[0]
        print self.spbL.get_pos(0)

    def advance(self):
        self.time_point += 1
        return self.gotoTime(self.time_point)

    def gotoTime(self, time_point):
        if time_point > self.mt.KD.num_steps:
            return False
        self.time_point = time_point
        for item in self.items:
            newPos = item.get_pos(self.time_point)
            if newPos == item.pos():
                return False
            item.setPos(newPos)
            if isinstance(item, GraphCentromere):
                for plugsite in item.plugsites:
                    newPos = plugsite.get_pos(self.time_point)
                    plugsite.setPos(newPos)
                    if plugsite.newPlug != plugsite.sim.state_hist[self.time_point]:
                        plugsite.newPlug = plugsite.sim.state_hist[self.time_point]
                        if plugsite.newPlug == 0:
                            brush = QtGui.QBrush(UNPLUGED_COLOR)
                        elif plugsite.sim.is_correct(self.time_point):
                            brush = QtGui.QBrush(GOOD_PLUGSITE_COLOR)
                        else:
                            brush = QtGui.QBrush(BAD_PLUGSITE_COLOR)
                        plugsite.color = brush
        return True

    def boundingRect(self):
        """
        As this is not easily changed (or not supposed to, we give a
        large box
        """
        height = 5.  # self.N * ( self.Mk * SITE_OFFSET + CH_OFFSET)
        width = 14.  # microns, shall be enough
        return QtCore.QRectF(-width / 2, - height / 2, width, height)

    def paint(self, painter, option, widget):
        painter.setBrush(QtCore.Qt.white)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.1))
        painter.drawRoundedRect(self.boundingRect(), 30, 100, QtCore.Qt.RelativeSize)


class GraphCentromere(QtGui.QGraphicsItem):

    def __init__(self, n, side, parent=None):

        QtGui.QGraphicsItem.__init__(self, parent=parent)
        self.graphcell = parent
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)

        N = int(self.graphcell.mt.KD.params['N'])
        Mk = int(self.graphcell.mt.KD.params['Mk'])
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
        self.n = n
        self.side = side
        self.ch = self.graphcell.mt.KD.chromosomes[n]
        if side == 0:
            x = self.ch.cen_A.pos
            self.traj = self.ch.cen_A.traj
        else:
            x = self.ch.cen_B.pos
            self.traj = self.ch.cen_B.traj
        self.y = (n - (N - 1) / 2.) * 0.4  # * ( Mk * SITE_OFFSET + CH_OFFSET)
        self.width = 0.2
        self.height = 0.1 * Mk  # Mk * SITE_OFFSET + CH_OFFSET
        self.newPos = QtCore.QPointF(x, self.y)
        self.plugsites = []
        self.setZValue(n)
        for m in range(Mk):
            self.plugsites.append(GraphPlugSite(m, self, parent))

    def get_pos(self, time_point):
        try:
            x = self.traj[time_point]
        except IndexError:
            x = self.traj[-1]
        return QtCore.QPointF(x, self.y)

    def advance(self):
        self.time_point += 1
        for item in self.items:
            newPos = item.get_pos(self.time_point)
            if newPos == item.pos():
                return False
            item.setPos(newPos)
        return True

    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(self.width, self.height,
                             2 * self.width, 2 * self.height)
        return self.path

    def paint(self, painter, option, widget):
        if self.ch.cen_A.is_attached() and self.ch.cen_B.is_attached():
            brush = QtGui.QBrush(CH_COLOR)
        else:
            brush = QtGui.QBrush(ACTIVE_SAC_COLOR)
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.newPos
        painter.drawEllipse(center, self.width, self.height)

    def boundingRect(self):
        adjust = 1.
        return QtCore.QRectF(-self.width - adjust, - self.height - adjust,
                             2 * self.width + adjust, 2 * self.height + adjust)


class GraphPlugSite(QtGui.QGraphicsItem):

    def __init__(self, m, kineto, parent=None):

        QtGui.QGraphicsItem.__init__(self, parent=parent)
        self.kineto = kineto
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
        Mk = int(self.kineto.graphcell.mt.KD.params['Mk'])
        self.m = m
        side = self.kineto.side
        if side == 0:
            self.sim = self.kineto.ch.cen_A.plugsites[m]
        else:
            self.sim = self.kineto.ch.cen_B.plugsites[m]

        self.width = 0.1
        self.height = self.kineto.height/Mk

        self.y = self.kineto.y + ( m - (Mk -1 )/2.) * 0.1
        x = self.sim.pos
        self.newPos = QtCore.QPointF(x, self.y)
        self.newPlug = self.sim.plug_state
        self.setZValue(self.kineto.n * (1 + m))

        if self.sim.plug_state == 0:
            brush = QtGui.QBrush(UNPLUGED_COLOR)
        elif self.sim.is_correct(0):
            brush = QtGui.QBrush(GOOD_PLUGSITE_COLOR)
        else:
            brush = QtGui.QBrush(BAD_PLUGSITE_COLOR)
        self.color = brush

    def get_pos(self, time_point):
        try:
            x = self.sim.traj[time_point]
        except IndexError:
            x = self.sim.traj[-1]
        return  QtCore.QPointF(x, self.y)

    def advance(self):
        self.time_point += 1
        for item in self.items:
            newPos = item.get_pos(self.time_point)
            if newPos == item.pos():
                return False
            item.setPos(newPos)
        return True

    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(-self.width / 2., - self.height / 2,
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
        return QtCore.QRectF(-self.width - adjust, - self.height - adjust,
                             2 * self.width + adjust, 2 * self.height + adjust)


class GraphSPB(QtGui.QGraphicsItem):

    def __init__(self, side, parent=None):
        QtGui.QGraphicsItem.__init__(self, parent=parent)
        self.graphcell = parent
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
        if side == 0:
            self.sim = self.graphcell.mt.KD.spbR
        else:
            self.sim = self.graphcell.mt.KD.spbL
        self.newPos = QtCore.QPointF(self.sim.pos, 0.)

    def get_pos(self, time_point):
        try:
            x = self.sim.traj[time_point]
        except IndexError:
            x = self.sim.traj[-1]
        y = 0.
        return QtCore.QPointF(x, y)

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
        super(NakedWidget, self).__init__()
        self.timerId = 0
        self.cell = GraphCell(mt)
        scene = QtGui.QGraphicsScene(self)
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


class InteractiveCellWidget(QtGui.QWidget):
    """
    """

    def __init__(self, mt):
        """
        """
        super(InteractiveCellWidget, self).__init__()

        self.timerId = 0
        self.view = ViewCellWidget(mt)
        self.ccw = ControlCellWidget(mt)
        self.connect(self.ccw.playButton, QtCore.SIGNAL('clicked()'),
                     self.play)
        self.connect(self.ccw.pauseButton, QtCore.SIGNAL('clicked()'),
                     self.pause)
        self.ccw.slider.valueChanged.connect(self.gotoTime)

        vbox = QtGui.QVBoxLayout()
        vbox.setSpacing(5)
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
            #self.pause()
            rects = [self.view.cell.boundingRect()]
            self.view.updateScene(rects)
            self.view.cell.gotoTime(time)
            #self.timerId = 0

    def play(self):
        # I'll re-implement this when user can move things in graphcell
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
        super(ControlCellWidget, self).__init__(parent)
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
        hbox.setSpacing(5)
        hbox.addWidget(self.playButton)
        hbox.addWidget(self.pauseButton)
        hbox.addWidget(self.slider)
        self.setLayout(hbox)


class ViewCellWidget(QtGui.QGraphicsView):

    def __init__(self, mt):
        super(ViewCellWidget, self).__init__()

        self.setRenderHint(QtGui.QPainter.Antialiasing)

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

        # Fill
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
