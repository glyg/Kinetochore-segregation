#!/usr/bin/python
# -*- coding: utf-8 -*-

# Handler for the parameters xml files -- to allow phase space exploration

import sys, os

from xml.etree.ElementTree import Element
from xml.etree.ElementTree import parse, tostring
from numpy import loadtxt, load

#Those strings should be respected in the xml file
SPRING_UNIT=u'pN/µm' 
DRAG_UNIT = u'pN.s/µm' 
LENGTH_UNIT = u'µm' 
FREQ_UNIT = u'Hz' 
FORCE_UNIT = u'pN' 
SPEED_UNIT = u'µm/s'


__all__ = ["ParamTree", "indent", "ResultTree"]

def indent(elem, level=0):
    '''utility to have a nice printing of the xml tree
    '''
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i



class ParamTree(object):

    '''
    This class defines the container for the simulation parameters.
    It wraps an ElementTree instance whose elements contains the
    name, value (as a string), description and unit of each parameter and a
    dictionnary that is used during the simulation.


    The 'value' attribute of the tree is not modified by adimentionalization,
    whereas the value in the dictionnary is changed.
    '''


    def __init__(self, filename = None, root = None, adimentionalized = True):

        if filename is not None:
            self.filename = filename
            source = file(filename, "r")
            self.tree = parse(source)
            self.root = self.tree.getroot()
            source.close()

        elif root is not None:
            self.root = root
        else:
            print 'A etree root or a filename should be provided'

        list=[]
        a = self.root.findall("param")
        for i in a:
            n = i.get("name")
            v = i.get("value")
            if '.' in v or 'e' in v:
                v = float(v)
            else:
                v = int(v)
            list.append((n, v))

        self.absolute_dic = dict(list)
        self.relative_dic = dict(list)
        if adimentionalized :
            self.adimentionalize()
            
    def has_unit(self, param, UNIT):

        unit_str = param.find("unit").text
        if unit_str == UNIT:
            return True
        else :
            return False

    def adimentionalize(self):

        '''
        This function scales everything taking dt as unit time, Vk as
        unit speed, and Fk as unit force. It relies on a correct
        definition of the units of the elements of the param tree,
        thus a correct spelling in the xml file, so please beware

        Note that the "value" attribute of the ElementTree instance are
        NOT modified. Also, applying this function on an already
        adimentionalized  dictionnary won"t change any thing
        '''

        Vk = self.absolute_dic["Vk"]
        Fk = self.absolute_dic["Fk"]
        dt = self.absolute_dic["dt"]
        d0 = self.absolute_dic["d0"]
        params = self.root.findall("param")

        for param in params:
            key = param.get("name")
            val = self.absolute_dic[key]
            if self.has_unit(param, SPRING_UNIT):
                val /= Fk
            elif self.has_unit(param, DRAG_UNIT):
                val *= Vk/Fk  
            elif self.has_unit(param, SPEED_UNIT):
                val /=  Vk
            elif self.has_unit(param, FREQ_UNIT):
                val *= dt
            elif self.has_unit(param, FORCE_UNIT):
                val /= Fk
            self.relative_dic[key] = val


    def change_dic(self, key, new_value, write = True, back_up = False,
                   verbose = True):
        '''
        Changes the Element tree and re-creates the associated dictionnary.
        If write is True, re_writes the parameters files
        (older version is backed up if back_up is True)
        
        new_value is absolute - it has units
        '''
        
        
        if self.absolute_dic is None:
            self.create_dic()
        a = self.root.findall('param')
        for item in a:
            if item.attrib['name'] == key:
                item.attrib['value'] = str(new_value)
                try:
                    self.absolute_dic[key] = new_value
                except KeyError:
                    print "Couldn't find the parameter %s" %key
                    return 0
                break
        try:
            self.adimentionalize()
        except KeyError:
            if 'metaph_rate' in self.absolute_dic.keys():
                self.relative_dic = self.absolute_dic
            else:
                print 'Oooups'
                raise()
            
        if write:
            indent(self.root)
            xf = open(self.filename, 'w+')
            if back_up:
                bck = self.filename+'.bck'
                xfb = open(bck,'w+')
                for line in xf:
                    xfb.write
                xfb.close()
                xf.seek(0)
                print "Backed up old parameter file as %s" %bck
            else :
                print "Warning : %s changed without back up" %self.filename
            xf.write(tostring(self))
            xf.close
            print "Changed  %s to %03f in file %s" %(key,
                                                     new_value,
                                                     self.filename)
        elif verbose:
            print "Warning: parameter %s changed but not written!" %key


        
class ResultTree(ParamTree):    

    def __init__(self, xmlfname = "resuts.xml"):

        xmlfname = os.path.abspath(xmlfname)
        ParamTree.__init__(self, xmlfname, adimentionalized = False)
        
        datafname = self.root.get("datafile")
        if not os.path.isabs(datafname):
            if '/' not in datafname:
                self.datafname = os.path.join(os.path.dirname(xmlfname), 
                                              datafname)
            else:
                self.datafname = os.path.join(os.path.dirname(xmlfname), 
                                              datafname.split('/')[-1])
                
        else:
            self.datafname = datafname
        if self.datafname is None or not os.path.isfile(self.datafname):
            raise ValueError, "Corresponding data file not specified"

        if self.datafname.endswith('.npy'):
            self.data = load(self.datafname)
        else:
            self.data = loadtxt(self.datafname, delimiter=' ',
                                comments = '#')        
    
    def get_spb_trajs(self):

        Rcol = None 
        Lcol = None
        for traj in self.tree.getiterator("trajectory"):
            if traj.get("name") == "rightspb":
                Rcol = traj.get("column")
            if traj.get("name") == "leftspb":
                Lcol = traj.get("column")
            if Rcol is not None and Lcol is not None:
                break
        spb_trajs = self.data.take((Rcol, Lcol), axis = 1)
        return spb_trajs

    def get_kineto_trajs(self):

        N = int(self.relative_dic['N'])
        cols = []
        for traj in self.tree.getiterator("trajectory"):
            if traj.get("name") == "rightkineto":
                rcol = traj.get("column")
                cols.append(rcol)
            if traj.get("name") == "leftkineto":
                lcol = traj.get("column")
                cols.append(lcol)
            if len(cols) == N * 2:
                break
        cols = tuple(cols)
        
        kineto_trajs = self.data.take(cols, axis = 1)
        return kineto_trajs
        
    def get_array(self, name, index = None):

        ''' returns the array of
        the corresponding trajectory

        name can be : {right | left}{spb | kineto | kMT | iMT}
        for spbs index is None
        for kinetos and iMTs, index is an int
        for kMTs index is a pair of integers.
        i.e.: get_array(rightiMT, 3), get_array(leftkMT, (0, 2))
        '''
        f = file(self.datafname)
        header = f.readline()
        for traj in self.tree.getiterator("trajectory"):
            if traj.get("name") == name :
                if index is None or traj.get("index") == str(index) :
                    col = int(traj.get("column"))
                    print col
                    return loadtxt(f, delimiter = ' ', usecols= (col, )) 
        
        print "trajectory not Found !"
        
    def get_all_trajs(self):

        f = file(self.datafname)
        header = f.readline()
        cols = []
        for traj in self.tree.getiterator("trajectory"):
            col = int(traj.get("column"))
            cols.append(col)

        cols = tuple(cols)
        
        return self.data.take(cols, axis = 1) 

    def get_all_pluged(self):

        cols = []
        for traj in self.tree.getiterator("numberpluged"):
            col = int(traj.get("column"))
            cols.append(col)

        cols = tuple(cols)
        return self.data.take(cols, axis = 1) 

        
    def get_all_mero(self):

        cols = []
        for traj in self.tree.getiterator("numbermero"):
            col = int(traj.get("column"))
            cols.append(col)

        cols = tuple(cols)
        return self.data.take(cols, axis = 1) 



    def get_all_plug_state(self):


        cols = []
        for state_hist in self.tree.getiterator("state"):
            col = int(state_hist.get("column"))
            cols.append(col)

        cols = tuple(cols)
        return self.data.take(cols, axis = 1) 
