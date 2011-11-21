#!/usr/bin/python
# -*- coding: utf-8 -*-
#from pylab import *
from numpy import *
from scipy import *
#import matplotlib.mlab as mlab
import sys, os


import pyximport
pyximport.install()

from spindle_dynamics import *
from simul_spindle import *
from xml_handler import *
from explore_params import *

from numpy import exp, logspace

try:
    paramfile = os.path.join(os.path.dirname(__file__), 'params.xml')
    measurefile = os.path.join(os.path.dirname(__file__), 'measures.xml')
except NameError: #Not importing from a module
    paramfile = 'params.xml'
    measurefile = 'measures.xml'

def load_paramtree():

    paramtree = ParamTree(paramfile)
    paramtree.create_dic()

    base = {}#paramtree.dic
    base.update({'dt':2.})#, 'fa':0.01, 'fd':0.005})

    pombe = base.copy()
    pombe.update({'N': 3, 'Mk':4, 'name':'Fission Yeast (Mk = 4)'})

    drosoS2 = base.copy()
    drosoS2.update({'N':4, 'Mk':15, 'name':'Drosophila S2 cell (Mk = 15)',
                    'l0':2., 'd0':1.})

    haemanthus = base.copy()
    haemanthus.update({'N':14, 'Mk':80, 'name':'Haemanthus endosperm (Mk = 80)',
                       'l0':10., 'd0':1.})

    hela = base.copy()
    hela.update({'N':23, 'Mk':30, 'name':'Homo sapiens (Mk = 30)'})

    ptk1 = base.copy()
    ptk1.update({'N':6, 'Mk':40, 'name':'Ptk1 (Mk = 40)', 'l0':2., 'd0':1.})

    newt = base.copy()
    newt.update({'N':5, 'Mk':24, 'name':'Newt lung cell (Mk = 24)',
                 'l0':2., 'd0':1.})

    organisms = [pombe, drosoS2,  ptk1,  haemanthus]
    return paramtree, organisms

paramtree, organisms = load_paramtree()
# #hela, drosoS2, pombe, haemanthus, ptk1, newt = organisms
pombe, drosoS2,  ptk1,  haemanthus = organisms
# #organisms = [drosoS2, pombe, haemanthus]

MEASURES = {}
pombe['measures'] = {'metaph_rate': 0.001,
                     'poleward_speed':  0.018,
                     'anaph_rate':0.022,
                     'metaph_k_dist':0.8,
                     'oi_dist':0.05,
                     'tau_o':10.,
                     'tau_i':10.,
                     'obs_d0':0.2}
ptk1['measures'] = {'metaph_rate': 0.001,
                    'poleward_speed':  0.018,
                    'anaph_rate':0.022,
                    'metaph_k_dist':2.1,
                    'oi_dist':0.05,
                    'tau_o':10.,
                    'tau_i':10.,
                    'obs_d0':1.}

#hela['measures'] = ptk1['measures']
drosoS2['measures'] = ptk1['measures']
haemanthus['measures'] = ptk1['measures']
#newt['measures'] = ptk1['measures']



def markov_plug(organism, num_steps):

    m = Metaphase(num_steps, paramfile = paramfile)
    for key, new_value in organism.items(): 
        m.paramtree.change_dic(key, new_value, write = False,
                               back_up = False, verbose = False)
    #we don't want to execute anaphase during this kind of simulation
    m.paramtree.change_dic('trans', num_steps, write = False,
                           back_up = False, verbose = False)
    
    m.__init__(num_steps, m.paramtree)
    
    kd = m.KD
    dt = kd.params['dt']
    
    for t in range(num_steps // dt):
        kd.plug_unplug()
        for ch in kd.chromosomes.values():
            ch.pluged_history.append(ch.pluged)
            ch.mero_history.append(ch.mero)

            
    kd.num_steps = num_steps //dt + 1

    return kd

def full_simul(organism, num_steps, reduced = True):

    m = Metaphase(num_steps, paramfile = paramfile,
                  measurefile = measurefile, reduce_p = reduced)
    for key, new_value in organism.items(): 
        m.paramtree.change_dic(key, new_value, write = False,
                               back_up = False, verbose = False)
    #we don't want to execute anaphase during this kind of simulation
    m.paramtree.change_dic('trans', num_steps, write = False,
                           back_up = False, verbose = False)


    m.__init__(num_steps, m.paramtree)
    m.simul()
    m.KD.num_steps = len(m.KD.spbR.traj)
    
    return m.KD
    
def get_history(kd, fmt = '-'):

    '''
    returns the history of the various attachment states.
    Takes a KinetochoreDynamics instance as unique argument

    output:
    plugs : a (2*N, num_steps) shaped array of plug_history
    meros : a (2*N, num_steps) shaped array of mero_history
    The following outputs are established on a CHROMOSOMES
    (i.e. kt pairs) basis
    num_plugs: the number of correctely plugged chromosomes
    num_synt: the number of syntelic chromosomes
    num_mono: the number of monotelic chromosomes
    num_mero: the number of merotelic chromosomes
    num_unplugs: the number of unpluged chromosomes
    Those outputs are returned in a dictionnary called defects:
    defects = {"amphitelic":num_plugs,
               "merotelic":num_meros,
               "monotelic":num_mono,
               "syntelic":num_synt
               "unattached":num_unplugs}

    Normaly num_plugs+num_synt+num_mono+num_mero+num_unplugs = N
    '''

    Mk = int(kd.params['Mk'])
    N = int(kd.params['N'])
    
    meros = []
    plugs = []
    time_lapse = arange(kd.num_steps) * kd.params['dt']

    for ch in kd.chromosomes.values():
        mh = array(ch.mero_history)
        meros.append(mh)
        ph = array(ch.pluged_history)
        plugs.append(ph)
        
    meros = hstack(meros)
    plugs = hstack(plugs)
    mh_b = meros.astype(bool)
    ph_b = plugs.astype(bool)

    #The number of amphitelic chromosomes
    #             Correct  Merotelic
    # -------------------------------
    #      Left  |  True     False    
    # AND  Right |  True     False  
    #OR
    #      Left  |  False    True
    # AND  Right |  False    True  


    correct_r = ph_b[:,::2] & logical_not(mh_b[:,::2])
    correct_l = ph_b[:,1::2] & logical_not(mh_b[:,1::2])
    and_x = correct_r & correct_l    

    mero_r = mh_b[:,::2] & logical_not(ph_b[:,::2])
    mero_l = mh_b[:,1::2] & logical_not(ph_b[:,1::2])
    and_x_mero = mero_r & mero_l    

    num_plug = (and_x | and_x_mero).sum(axis = 1)

    #The number of syntelic chromosomes
    #        Correct  Merotelic OR  Correct  Merotelic
    #------------------------------------------------
    # Left  |  True     False        False    True
    # Right | False     True         True     False

    xor_l = mh_b[:,::2] ^ ph_b[:,::2]   #left
    xor_r = mh_b[:,1::2] ^ ph_b[:,1::2] #right 
    xor_x =  ph_b[:,::2] ^ ph_b[:,1::2] #crossed
    xor_t = xor_l & xor_r & xor_x
    num_synt = xor_t.sum(axis = 1)

    #The number of merotelic chromosomes
    and_r = mh_b[:,::2] & ph_b[:,::2]   #          Correct  Merotelic
    and_l = mh_b[:,1::2] & ph_b[:,1::2] #     Left   True     True    
    or_x = and_r | and_l                # OR  Right  True     True  
    num_mero = or_x.sum(axis = 1)
    
    #The number of monotelic chromosomes
    nand_r = logical_not(mh_b[:,::2] | ph_b[:,::2]) & xor_r   #            Correct  Merotelic
    nand_l = logical_not(mh_b[:,1::2] | ph_b[:,1::2]) & xor_l #      Left   False     False    
    xor_x = nand_r ^ nand_l                                   # XOR  Right  False     False  
    num_mono = xor_x.sum(axis = 1)

    #The number of unattached chromosomes
    #             Correct  Merotelic
    # -------------------------------
    #      Left  |  False     False    
    # AND  Right |  False     False  
    and_x = nand_r & nand_l    
    num_unat = and_x.sum(axis = 1)

    #Put all this in a dictionary
    defects = {'amphitelic':num_plug,
               'merotelic':num_mero,
               'monotelic':num_mono,
               'syntelic':num_synt,
               'unattached':num_unat }

    return meros , plugs, defects

    

def explore_aurora(organism, num_steps, num_ech):

    #auroras = logspace(-1.5, .5, 20)
    auroras = logspace(-1.5, .2, 40)
    all_defects = {}
    for aurora in auroras:
        organism['aurora'] = aurora
        ms, ps, defects = explore_one(organism, num_steps,
                                      num_ech, display = False, full = True)
        for key, value in defects.items():
            #value_mean = value.mean(axis = 0)
            try:
                all_defects[key] = vstack((all_defects[key], value))
            except KeyError:
                all_defects[key] = value
    return all_defects

def explore_orientation(organism, num_steps, num_ech):

    organism['aurora'] = 0
    fds = logspace(-2, 0., 40)
    all_defects = {}

    for fd in fds:
        organism['fd'] = 0.05 * fd
        ms, ps, defects = explore_one(organism, num_steps,
                                      num_ech, display = False, full = True)
        for key, value in defects.items():
            #value_mean = value.mean(axis = 0)
            try:
                all_defects[key] = vstack((all_defects[key],value))# value_mean))
            except KeyError:
                all_defects[key] = value

    return all_defects
    

def explore_one(organism, num_steps, num_ech, display = False, full = False):
    '''
    if full is True, will perform a full simulation rather than
    only the plug/unplug sequence
    '''

    N = int(organism['N'])
    dt = organism['dt']
    all_defects = {}
    for i in range(num_ech):
        if full:
            mp = full_simul(organism, num_steps, reduced = True)

        else:
            mp =  markov_plug(organism, num_steps )

        # meros, plugs, unplugs, num_meros, num_synt = get_history(mp, fmt = 'k,',
        #                                                           display = False)
        meros, plugs, defects = get_history(mp)
        if i == 0:
            all_meros = meros
            all_plugs = plugs
            all_defects = defects
        else:
            all_meros = hstack((all_meros, meros))
            all_plugs = hstack((all_plugs, plugs))
            for defect in defects.keys():
                all_defects[defect] = hstack((all_defects[defect],
                                              defects[defect]))
            
    if display: #FIXME
        print "Broken now"
        #elapsed = arange(num_steps + dt, step = dt)
        #figure()
        # figure(101)
        # errorbar(elapsed, all_plugs.mean(axis = 0),
        #          yerr = all_plugs.std(axis = 0)/sqrt(num_ech),
        #          fmt = '-', label = organism['name'])
        # xlabel('Time (s)')
        # ylabel('pluged kt frequency')
        # figure(102)
        # errorbar(elapsed, all_unplugs.mean(axis = 0),
        #          yerr = all_unplugs.std(axis = 0)/sqrt(num_ech),
        #          fmt = '-', label = organism['name'])
        # xlabel('Time (s)')
        # ylabel('Unpluged kt frequency')
    return all_meros, all_plugs, all_defects

#BROKEN
def explore_all(num_steps, num_ech, display = False, base_name = None,
                full = False):

    for organism in organisms:
        print 'Starting organism %s' %organism['name']

        auroras = logspace(-1, 1, 10)
        n = 0
        for aurora in auroras:
            organism['aurora'] = aurora
        
            num_ech_eff = int(  num_ech /  float( organism['N'] ))
            if num_ech_eff < 2:
                print 'Too few samples'
                num_ech_eff = 2
            organism['meros'], organism['plugs'], organism['unplugs'],\
            organism['all_num_meros'] = explore_one(organism, num_steps,
                                                    num_ech_eff,
                                                    display = display,
                                                    full = full)


            if base_name is not None:
                mname = '%s_meros_%s_auro%i.txt.gz' %(base_name,
                                                      organism['name'].split()[0], n)
                savetxt(mname, organism['meros'])
                mname = '%s_plugs_%s_auro%i.txt.gz' %( base_name,
                                                       organism['name'].split()[0], n)
                savetxt(mname, organism['plugs'])
                mname = '%s_unplugs_%s_auro%i.txt.gz' %(base_name,
                                                        organism['name'].split()[0], n)
                savetxt(mname, organism['unplugs'])
                mname = '%s_nmeros_%s_auro%i.txt.gz' %(base_name,
                                                       organism['name'].split()[0], n)
                savetxt(mname, organism['all_num_meros'])
                n += 1

        print 'Finished organism %s' %organism['name']        
    return organisms 




def power_err((N0, Ninf, tau), t, plugs):

    ffunc = N0 + Ninf * (1 - exp(-t/tau))
    return plugs  - ffunc

def exp_err((N0, tau), t, meros):
    ffunc = N0 * exp( - t/ tau)
    return meros - ffunc
    

def fit_plugs(organism, plugs, num_steps):

    elapsed = arange(num_steps + 1)*organism['dt']
    Mk = organism['Mk']
    init = Mk/2, Mk/2, 100
    fit_param = leastsq(power_err, init, (elapsed, plugs))

    print fit_param

def fit_meros(organism, meros, num_steps):

    elapsed = arange(num_steps + 1)*organism['dt']
    Mk = int(organism['Mk'])
    init = Mk/2, 100
    fit_param = leastsq(exp_err, init, (elapsed, meros))

    print fit_param


def reload_sim(base_name, attch = 'meros', term = '.txt'):
    
    histories = {}
    for organism in organisms:
        fname = '%s_%s_%s%s' %((base_name, attch, organism['name'].split()[0],
                                term))
        histories[organism['name']] = loadtxt(fname)

    return histories

def get_biorientation_times():

    for organism in organisms:
        organism['bio_times'] = []
    for aurora in auroras:
        for organism in organisms:
            m_p = vstack(organism[aurora]).mean(axis = 0)
            organism['bio_times'].append(more_than(m_p, 0.01))
            

            

def reload_aurora(attch):

    auroras = logspace(-1, 1, 10)
    for n, aurora in enumerate(auroras):
        base_name = 'attach_only/aurora'
        term = '_auro%i.txt.gz' %n
        reload_allsims(base_name, attch, term = term)
        for organism in organisms:
            organism[aurora] = organism[attch]

def more_than(a, lim, axis = None):
    ''' 
    returns the index of the last value of a higher than lim
    '''
    a = asarray(a)
    
    if axis is None:
        if a.min() >= lim:
            return 0
        elif a.max() < lim:
            return -1
        else:
            return find(a > lim).max()
    else:
        sw_a = a.swapaxes(0, axis)
        lps = zeros(sw_a.shape[0]) 
        for n, col in enumerate(sw_a):
            lps[n] = more_than(col, lim)

        return lps


    
def reload_allsims(base_name, attch = 'meros', term = '.txt.gz'):

    for organism in organisms:
        organism[attch] = []
        
 
    for i in range(1, 5):

        base_name_i = '%s%i' %(base_name, i)
        histories = reload_sim(base_name_i, attch, term = term)
        for organism in organisms :
            organism[attch].append(histories[organism['name']])
    
    for organism in organisms:
        all_hist = vstack(organism[attch])
        dt  = organism['dt']
        elapsed = arange(all_hist.shape[1]) * dt
        avrg = all_hist.mean(axis = 0)
        error = all_hist.std(axis = 0)/sqrt(all_hist.shape[0])
        upper = avrg + error
        lower = avrg - error
#        xs, ys = mlab.poly_between(elapsed[lower >0], lower[lower >0],
#                                   upper[lower >0])
#        figure(102) 
#        plot(elapsed, avrg, label = organism['name'])
#        fill(xs, ys, alpha = 0.1, lw = 0)
#        figure(101)
#        semilogy(elapsed, avrg, label = organism['name'])
 
# if __name__ == "__main__":


#     base_name = sys.argv[1]
#     defects = explore_aurora(pombe, 2000, 250)
#     for key, value in defects:
#         mname = '%spombe_%s_auro.npy' %(base_name, key)
#         save(mname, value)

#     defects = explore_orientation(pombe, 2000, 250)
#     for key, value in defects:
#         mname = '%spombe_%s_orient.npy' %(base_name, key)
#         save(mname, value)
