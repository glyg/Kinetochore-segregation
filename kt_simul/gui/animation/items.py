#-*- coding: utf-8 -*-

from PySide import QtCore, QtGui

SITE_OFFSET = 0.2  # Vertical distance between attachment sites
CH_OFFSET = 0.4  # Vertical distance between chromosomes

SPB_COLOR = QtGui.QColor(255, 0, 0)  # Red

CH_COLOR = QtGui.QColor(0, 100, 100, alpha=200)  # Green
ACTIVE_SAC_COLOR = QtGui.QColor(100, 10, 100, alpha=200)  # Purple

GOOD_PLUGSITE_COLOR = QtGui.QColor(0, 250, 20, alpha=200)  # Green
BAD_PLUGSITE_COLOR = QtGui.QColor(255, 0, 0, alpha=255)  # Red
UNPLUGED_COLOR = QtGui.QColor(0, 20, 250, alpha=200)  # Blue


class CellItem(QtGui.QGraphicsItem):
    """
    This is the parent item containing all the objects within the cell
    """

    def __init__(self, metaphase, parent=None):
        QtGui.QGraphicsItem.__init__(self, parent=parent)

        self.N = int(metaphase.KD.params['N'])
        self.Mk = int(metaphase.KD.params['Mk'])
        self.mt = metaphase  # Metaphase instance
        self.items = []
        self.spbR = SPBItem(0, parent=self)
        self.items.append(self.spbR)
        self.spbL = SPBItem(-1, parent=self)
        self.items.append(self.spbL)

        for n in range(self.N):
            left_kineto = CentromereItem(n, -1, parent=self)
            self.items.append(left_kineto)
            right_kineto = CentromereItem(n, 0, parent=self)
            self.items.append(right_kineto)

        self.time_point = -1
        self.gotoTime(0)

    def advance(self):
        self.time_point += 1
        return self.gotoTime(self.time_point)

    def gotoTime(self, time_point):
        if time_point >= self.mt.KD.num_steps - 1:
            return False

        self.time_point = time_point
        for item in self.items:
            newPos = item.get_pos(self.time_point)
            if newPos == item.pos():
                return False
            item.setPos(newPos)
            if isinstance(item, CentromereItem):
                self.update_plugsite_state(item)

        return True

    def update_plugsite_state(self, centromere):
        """
        Update plugsite color according to their state
        """
        for plugsite in centromere.plugsites:
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


class CentromereItem(QtGui.QGraphicsItem):

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
            self.plugsites.append(PlugSiteItem(m, self, parent))

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


class PlugSiteItem(QtGui.QGraphicsItem):

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
        self.height = self.kineto.height / Mk

        self.y = self.kineto.y + (m - (Mk - 1) / 2.) * 0.1
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


class SPBItem(QtGui.QGraphicsItem):

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
        return QtCore.QRectF(-0.2 - adjust, -0.6 - adjust,
                             0.4 + adjust, 1.2 + adjust)
