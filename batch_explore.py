#!/usr/bin/python
# -*- coding: utf-8 -*-
#from pylab import *
from numpy import *
from scipy import *
#import matplotlib.mlab as mlab
import sys, os
import time

import pyximport
pyximport.install()

from kt_simul.simul_spindle import *
from kt_simul.xml_handler import *
from kt_simul.attachment_state import *

try:
    paramfile = os.path.join(os.path.dirname(__file__), 'params.xml')
    measurefile = os.path.join(os.path.dirname(__file__), 'measures.xml')
except NameError: #Not importing from a module
    paramfile = 'params.xml'
    measurefile = 'measures.xml'


def full_simul(new_params, num_steps,  plug = 'monotelic'):

    m = Metaphase()
    for key, new_value in new_params.items(): 
        m.paramtree.change_dic(key, new_value, write = False,
                               back_up = False, verbose = False)
    #we don't want to execute anaphase during this kind of simulation
    m.paramtree.change_dic('trans', num_steps, write = False,
                           back_up = False, verbose = False)

    m.__init__(num_steps, m.paramtree,  plug = plug)
    m.simul()
    m.KD.num_steps = len(m.KD.spbR.traj)
    
    return m
    

def explore_2D(num_steps, num_ech, plug):

    new_params = {}

    pcs1s = linspace(0., 1., 41)[::-1]
    auroras = logspace(-1.5, .2, 41)
    
    all_defects = {}
    all_balance = None


    for pcs1 in pcs1s:
        for aurora in auroras:
            new_params['aurora'] = aurora
            new_params['orientation'] = pcs1
            
            defects, balance = explore_one(new_params, num_steps,
                                           num_ech, plug = plug)
            for key, value in defects.items():
                try:
                    all_defects[key] = vstack((all_defects[key],value))
                except KeyError:
                    all_defects[key] = value
            if all_balance is None:
                all_balance = balance
            else:
                all_balance = vstack((all_balance,balance))
                

    return all_defects, all_balance

        

def explore_one(new_params, num_steps, num_ech, plug):
    '''
    if full is True, will perform a full simulation rather than
    only the plug/unplug sequence
    '''

    all_defects = {}
    for i in range(num_ech):

        mp = full_simul(new_params, num_steps, plug = plug)
        defects, were_defects = defect_histories(mp.KD)
        balance = balance_histories(mp.KD)
        if i == 0:
            all_defects = defects
            all_balance = balance
        else:
            for defect in defects.keys():
                all_defects[defect] = hstack((all_defects[defect],
                                              defects[defect]))
            all_balance = hstack((all_balance, balance))

    return all_defects, all_balance

 
if __name__ == "__main__":


    base_name = sys.argv[1]

    plugs = ['monotelic', 'null', 'merotelic', 'amphitelic','syntelic']
    
    for plug in plugs:
        defects, balance = explore_2D(800, 15, plug = plug)
        for key, value in defects.items():
            mname = '%sdefect_%s_%s.npy' %(base_name, key, plug)
            save(mname, value)
        mname = '%sbalance_%s.npy' %(base_name, plug)
        save(mname, balance)

    