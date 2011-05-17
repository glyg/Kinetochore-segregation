#!/usr/bin/python
# -*- coding: utf-8 -*-
#from pylab import *
from numpy import *
#from scipy import *
#import matplotlib.mlab as mlab
import sys, os
import time

import pyximport
pyximport.install()

sys.path.append('/home/gay/python')
sys.path.append('/home/gay/python/lib')

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
    

def explore_2D(pcs1s, auroras, num_steps, num_ech, plug, logfile = None):

    new_params = {}

    
    all_defects = {}
    all_balance = None
    all_trans_mat = None

    for pcs1 in pcs1s:
        for aurora in auroras:
            new_params['aurora'] = aurora
            new_params['orientation'] = pcs1
            new_params['dt'] = 2.
            defects, balance, trans_mat = explore_one(new_params, num_steps,
                                                      num_ech, plug = plug)
            
            logfile.write('ran plug = %s, pcs1 = %03f, aurora = %03f\n' %(plug, pcs1, aurora))
                
            for key, value in defects.items():
                try:
                    all_defects[key] = vstack((all_defects[key],value))
                except KeyError:
                    all_defects[key] = value
            if all_balance is None:
                all_balance = balance
                all_trans_mat = trans_mat
            else:
                all_balance = vstack((all_balance,balance))
                all_trans_mat = vstack((all_trans_mat, trans_mat))

    return all_defects, all_balance, all_trans_mat

        

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
        trans_mat = transition_matrix(mp.KD)
        if i == 0:
            all_defects = defects
            all_balance = balance
            all_trans_mat = trans_mat
        else:
            for defect in defects.keys():
                all_defects[defect] = hstack((all_defects[defect],
                                              defects[defect]))
            all_balance = hstack((all_balance, balance))
            all_trans_mat = hstack((all_trans_mat, trans_mat))
            del mp, defects, were_defects, balance, trans_mat
            
    return all_defects, all_balance, all_trans_mat

 
if __name__ == "__main__":


    base_name = sys.argv[1]

<<<<<<< HEAD
    pcs1s = linspace(0., 1., 21)[::-1]
    auroras = logspace(-2, .2, 41)
=======
    pcs1s = linspace(0.95, 1., 5)[::-1]
    auroras = logspace(-2.5, .2, 41)
>>>>>>> 147a9742874e4606cdb343b90ebd239ab9f93e1e
    pcs1name = '%s_pcs1s.npy' %base_name
    auroraname = '%s_auroras.npy' %base_name
    save(pcs1name, pcs1s)
    save(auroraname, auroras)
    logfile_name = '%s_log.txt'  %base_name
    logfile = open(logfile_name, 'w+')

<<<<<<< HEAD
    plugs = ['null', 'merotelic', 'random', 'amphitelic', 'monotelic', 'syntelic']
    for plug in plugs:
        defects, balance, trans_mat = explore_2D(pcs1s, auroras, 800, 25, plug, logfile)
=======
    print 'saved auras'
    plugs = ['monotelic', 'null', 'random', 'merotelic']
    for plug in plugs:
        defects, balance, trans_mat = explore_2D(pcs1s, auroras, 800, 100, plug = plug)
>>>>>>> 147a9742874e4606cdb343b90ebd239ab9f93e1e
        for key, value in defects.items():
            mname = '%sdefect_%s_%s.npy' %(base_name, key, plug)
            save(mname, value)
        mname = '%sbalance_%s.npy' %(base_name, plug)
        save(mname, balance)
        mname = '%stransition_%s.npy' %(base_name, plug)
        save(mname, trans_mat)
        del defects, balance, trans_mat

    
