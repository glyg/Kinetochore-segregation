#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:Guillaume Gay<elagachado AT  gmail DOT com>
## Commentary:


'''
QtWidgets to set the values of the parameters and the measures
used in the simulation
'''

import sys
from PyQt4 import QtGui, QtCore

sys.path.append('/home/guillaume/Python')
from kt_simul.xml_handler import *
from numpy import log10, floor

__all__ = ['SetParameters', 'SetMeasures']

class ParamBox(QtGui.QSpinBox):
    def __init__(self, param, parent = None):
        QtGui.QSpinBox.__init__(self, parent)
        self.param = param
    #I'm using PyQt specific short circuit signals, see for exemple:
    #http://techbase.kde.org/Development/Tutorials/Python_introduction_to_signals_and_slots
    def my_emiter(self, val):
        self.emit(QtCore.SIGNAL('paramValueChanged'), self.param, val)

class DoubleParamBox(QtGui.QDoubleSpinBox):
    def __init__(self, param, parent = None):
        QtGui.QDoubleSpinBox.__init__(self, parent)
        self.param = param
        
    def my_emiter(self, val):
        self.emit(QtCore.SIGNAL('doubleparamValueChanged'), self.param, val)



class SetParameters(QtGui.QWidget):

    def __init__(self, paramtree, parent=None, adimentionalized = True):
        QtGui.QWidget.__init__(self, parent)

        self.setWindowTitle('Parameters')
        self.paramtree = paramtree
        self.parameters = self.paramtree.create_dic(adimentionalized)
        self.isModified = False

        hbox = QtGui.QHBoxLayout()
        
        vbox1 = QtGui.QVBoxLayout()
        vbox2 = QtGui.QVBoxLayout()
        k = 0

        params = self.paramtree.root.findall("param")
        nb_params = len(params)

        for param in params:
            description = param.find("description").text
            p_step = param.get("step")
            # If step is 0, the parameter is not displayed 
            # (it is set from the measures, via reduce_params)
            if p_step == '0':
                continue 
            unit = param.find('unit').text

            value = param.get("value")
            name = param.get("name")

            p_min = param.get("min")
            p_max = param.get("max")
            label = QtGui.QLabel(description, self)
            if not '.' in value:
                value = int(value)
                p_min = int(p_min)
                p_max = int(p_max)
                p_step = int(p_step)
                
                spinbox = ParamBox(param, self)
                self.connect(spinbox,  QtCore.SIGNAL('valueChanged(int)'), spinbox.my_emiter)
                self.connect(spinbox, QtCore.SIGNAL('paramValueChanged'), self.print_new_val)
                self.connect(spinbox, QtCore.SIGNAL('valueChanged(int)'), self.setModified)
                
            else:
                value = float(value)
                p_min = float(p_min)
                p_max = float(p_max)
                p_step = float(p_step)

                if value != 0 :
                    dec =  - int(floor(log10(abs(value)))) + 2
                else : 
                    dec = 2
                spinbox = DoubleParamBox(param, self)
                spinbox.setDecimals(dec)
                self.connect(spinbox,  QtCore.SIGNAL('valueChanged(double)'), spinbox.my_emiter)
                self.connect(spinbox, QtCore.SIGNAL('doubleparamValueChanged'), self.print_new_val)
                self.connect(spinbox, QtCore.SIGNAL('valueChanged(float)'), self.setModified)

            spinbox.setSingleStep(p_step)
            spinbox.setRange(p_min,p_max)
            spinbox.setSuffix(' '+unit)
            spinbox.setValue(value)

            vbox1.addWidget(label)
            vbox2.addWidget(spinbox)

        hbox.addLayout(vbox1)
        hbox.addLayout(vbox2)

        vbox1.setSizeConstraint(QtGui.QLayout.SetNoConstraint)
        vbox2.setSizeConstraint(QtGui.QLayout.SetNoConstraint)
        self.setLayout(hbox)
        hbox.setSizeConstraint(QtGui.QLayout.SetNoConstraint)

    def setModified(self, val):
        if val:
            self.isModified = True
        else :
            self.isModified = False

    def print_new_val(self, param, val):

        name = param.get("name")
        old_value = param.get('value')
        param.set("value", str(val))
        self.paramtree.change_dic(name, val, write = False, back_up = False, verbose = False)


class SetMeasures(SetParameters):

    def __init__(self, measuretree, parent = None):
        
        SetParameters.__init__(self, measuretree, parent, adimentionalized = False)
        self.setWindowTitle('Measures')
        
        params = self.paramtree.root.findall("param")
        nb_params = len(params)



if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    if len(sys.argv) > 1:
        paramfile = sys.argv[1]
        paramtree = ParamTree(paramfile)
    else :
        print "usage python param_seter.py parameter_file"
        sys.exit(0)
    setparameters = SetParameters(paramtree)
    scrollArea = QtGui.QScrollArea()
    scrollArea.setWidget(setparameters)
    scrollArea.show()

    sys.exit(app.exec_())
