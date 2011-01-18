#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:uillaume Gay<elagachado AT  gmail DOT com>
## Commentary:

from PyQt4 import QtCore, QtGui


sys.path.append('/home/guillaume/Python')

from kineto_simulation import SigMetaphase 

SITE_OFFSET = 0.2 #Vertical distance between attachment sites
CH_OFFSET = 0.4 # Vertical distance between chromosomes

SPB_COLOR = QtGui.QColor(1,0,1)
CH_COLOR = QtGui.QColor(0.2,0.2,0.2, alpha = 200)
GOOD_PLUGSITE_COLOR = QtGui.QColor(0,1,0, alpha = 200)
BAD_PLUGSITE_COLOR = QtGui.QColor(1,0,0, alpha = 255)

class PlayGround(QtGui.QGraphiscScene):

    def __init__(self, spindle, parent = None):
        QtGui.QGraphiscScene.__init__(0, 0, 600, 200, parent = parent)
        self.addItem(spindle)
        
        # self.spindle = spindle
        # self.draw_spindle()
        

        
    def draw_spindle(self):
        pass 
        


class GraphSpindle(QtGui.QGraphicsItem):
    '''
    This is the parent item containing all the objects within the spindle
    '''
    
    def __init__(self, mt, parent = None):

        QtGui.QGraphicsItem.__init__()
        self.mt = mt # A SimMetaphase instance
        

    def boundingRect(self):
        '''
        As this is not easily changed (or not supposed to, we give a
        large box
        '''
        N = self.mt.KD.params['N']
        Mk = self.mt.KD.params['Mk']
        
        height = N * ( Mk * SITE_OFFSET + CH_OFFSET)
        width = 12 # microns, shall be enough
        
        return QtCore.QRectF(-width/2, -height/2, width, height)

class GraphSPB(QtGui.QGraphicsItem):

    def __init__(self, parent = graphspindle):

        QtGui.QGraphicsItem.__init__()


    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(-0.1, -0.3, 0.2, 0.6)
            
    def paint(self):
        

    def boundingRect(self):
        
class GraphicsView(QtGui.QGraphicsView):

    def __init__(self):
        

