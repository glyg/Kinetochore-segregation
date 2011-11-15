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

# sys.path.append('/home/guillaume/python')
# sys.path.append('/home/guil/python/lib')

from kt_simul.simul_spindle import *
from kt_simul.xml_handler import *
from kt_simul.attachment_state import *

try:
    paramfile = os.path.join(os.path.dirname(__file__), 'params.xml')
    measurefile = os.path.join(os.path.dirname(__file__), 'measures.xml')
except NameError: #Not importing from a module
    paramfile = 'params.xml'
    measurefile = 'measures.xml'


def full_simul(new_params, plug = 'monotelic'):

    m = Metaphase()
    for key, new_value in new_params.items(): 
        m.paramtree.change_dic(key, new_value, write = False,
                               back_up = False, verbose = False)
    #we don't want to execute anaphase during this kind of simulation
    m.paramtree.change_dic('trans', num_steps, write = False,
                           back_up = False, verbose = False)

    m.__init__(paramtree = m.paramtree,  plug = plug)
    m.simul()
    #m.KD.num_steps = len(m.KD.spbR.traj)
    
    return m
    

def explore_2D(pcs1s, auroras, num_steps, num_ech, plug, dt = 1):

    new_params = {}
    # all_defects = {}
    # all_balance = None
    # all_trans_mat = None

    logfile = file('%s.log' %base_name, 'w+')

    for n, pcs1 in enumerate(pcs1s[20:]):
        for m, aurora in enumerate(auroras):
            new_params['aurora'] = aurora
            new_params['orientation'] = pcs1
            new_params['span'] = num_steps * dt
            new_params['dt'] = dt

            for i in range(num_ech):
                mp = full_simul(new_params, plug = plug)
                xmlfname = '%s_res_pcs1-%03i_auroras-%03i_%03i.xml' %(base_name, n+20, m, i)
                datafname = '%s_data_pcs1-%03i_auroras-%03i_%03i.npy' %(base_name, n+20, m, i)
                mp.write_results(xmlfname, datafname)
                del mp
            
            logfile.write('ran plug = %s, pcs1 = %03f, aurora = %03f\n' %(plug, pcs1, aurora))

    return 0
            
    #         for key, value in defects.items():
    #             try:
    #                 all_defects[key] = vstack((all_defects[key],value))
    #             except KeyError:
    #                 all_defects[key] = value
    #         if all_balance is None:
    #             all_balance = balance
    #             all_trans_mat = trans_mat
    #         else:
    #             all_balance = vstack((all_balance,balance))
    #             all_trans_mat = vstack((all_trans_mat, trans_mat))

    # return all_defects, all_balance, all_trans_mat


def explore_aurora2D(fd0s, auroras, num_steps, num_ech, plug):

    new_params = {}
    for i in range(num_ech):
        for n, fd0 in enumerate(fd0s):
            for m, aurora in enumerate(auroras):
                new_params['aurora'] = aurora
                new_params['fd0'] = fd0
                new_params['dt'] = 2.

                mp = full_simul(new_params, num_steps, plug = plug)
                xmlfname = '%s_res_fd0-%03i_auroras-%03i_%03i.xml' %(base_name, n, m, i)
                datafname = '%s_data_fd0-%03i_auroras-%03i_%03i.npy' %(base_name, n, m, i)
                mp.write_results(xmlfname, datafname)
                del mp

    return 0

                



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
    t0 = time.time()
    #fd0s = linspace(0.01, 1., 21)
    pcs1s = linspace(0., 1., 31)
    auroras = logspace(-3., .2, 51)

    #fd0sname = '%s_fd0s.npy' %base_name
    #save(fd0sname, fd0s)
    pcs1sname = '%s_pcs1s.npy' %base_name
    save(pcs1sname, pcs1s)
    auroraname = '%s_auroras.npy' %base_name
    save(auroraname, auroras)

    plug = 'random'
    num_steps = 800
    dt = 1.
    num_ech = 100
    explore_2D(pcs1s, auroras, num_steps, num_ech, plug, dt)
    print 'time: %.3f' %(time.time() - t0)

    ##### fin prévue dimanche 14h ####
