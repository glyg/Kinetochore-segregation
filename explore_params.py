#!/usr/bin/python
# -*- coding: utf-8 -*-


import sys


from simul_spindle import *
from eval_simul import *

from numpy import logspace, linspace, log10, ndindex
from pylab import Figure, figure, plot, semilogx, show, pcolor
from time import time

def check_all(param_name, span, num_ech = 50, repeats = 10,
              logscale = False, reduced = True, *args):
    """
    This will just run the simulation repeats times for each of the
    num_ech values of the parameter param_name. It will run the
    observations given by Metaphsase.evalute and record them.
    It will also record the report attribute  for inspection

    span is a tuple (inf, sup), relative to the parameter default
    value (i.e, if (inf, sup) = (0.1, 10.), and the parameter is 8,
    we'll look from 0.8 to 80

    """

    observations = {'anaphase_rate':[],
                    'metaph_rate': [],
                    'metaph_k_dist': [],
                    'pitch': [],
                    'poleward_speed': [],
                    'kt_rms_speed':[],
                    'times_of_arrival':[],
                    'pluged_stats': []}

    reports = []
    param_vals = spaced(span, logscale, num_ech)    
    t0 = time()
    T = Metaphase(800)
    back = T.KD.params[param_name]
    for val in param_vals :
        
        for rep in range(repeats):
            T.paramtree.change_dic(param_name, val * back,
                                   write = False, verbose = False)
            if reduced:
                reduce_params(T.paramtree, MEASURES )
            T.__init__(800, T.paramtree)
            T.simul()
            for key, result in T.observations.items():
                observations[key].append(result)
            reports.append(T.report)
            print 'done %s = %.2f test number %d' %(param_name, val,
                                                    rep + 1)
    for obs, results in observations.items():
        observations.update({obs:array(results)})

    
    elt = (time() - t0)/float(num_ech * repeats)
    print "mean elapsed time by run : %2.5f seconds " %elt
    param_vals *= back
    param_vals = array(param_vals)
    return (param_vals, observations, reports)
    

def check_allND(param_names, spans, num_echs, repeats, 
                logscale = False, reduced = True, *args):
    '''
    Same as check_all except param_names, spans, and num_echs are
    lists whose elements are the arguments of check_all.

    The simulation will be ran repeats times on each cell of an array
    with num_echs shape. Thus if len(param_names) = 3 and num_echs =
    [4, 5, 4], we will look at every 4x5x4=80 possible cases (i.e:
    beware it might become resource intensive rapidly)

    '''

    if len(param_names) != len(spans) or len(spans) != len(num_echs):
        print 'Ill formed arguments'
        return 

    
    observations = {'anaphase_rate':[],
                    'metaph_rate': [],
                    'metaph_k_dist': [],
                    'pitch': [],
                    'poleward_speed': [],
                    'kt_rms_speed':[],
                    'times_of_arrival':[],
                    'pluged_stats': []}

    param_vals = []
    for span, num_ech in zip(spans, num_echs):
        param_vals.append(spaced(span, logscale, num_ech))
    param_vals = array(param_vals)

    ndi = ndindex(tuple(num_echs))
    t0 = time()
    T = Metaphase(800)
    
    backs = {}
    for n, param in enumerate(param_names):
        backs[n] = T.KD.params[param]
    

    for index in ndi:
        for n, param in enumerate(param_names):
            val = param_vals[n,index[n]]
            T.paramtree.change_dic(param, val * backs[n],
                                   write = False, verbose = False)
        if reduced:
            reduce_params(T.paramtree, MEASURES )
        
        for rep in range(repeats):
            T.__init__(800, T.paramtree)
            T.simul()
            for key, result in T.observations.items():
                observations[key].append(result)

    for obs, results in observations.items():
        observations.update({obs:array(results)})

    
    for n, p_vals in enumerate(param_vals):

        p_vals *= backs[n]
    

    return (param_names, param_vals, observations)
    



def check_carac(carac, param_name, span, num_ech = 50,
                logscale = False, *args):
    ''' This a generic function to check dependance of a caracteristic
    on a given parameter. returns a 2D array [[param_vals],[observations]]
    
    carac is a function returning a single value (e.g. an elongation
    rate). The only requirement for this function is to take a
    Kineto_Dynamics instance as it"s first argument. param_name is
    valid parameter name, span is a tuple (inf, sup), relative to the
    parameter default value (i.e, if (inf, sup) = (0.1, 10.), and the
    parameter is 8, we"ll look from 0.8 to 80) num_ech : number of
    samples.
    '''
    
    observations = []



    param_vals = spaced(span, logscale, num_ech)    
    t0 = time()
    T = Metaphase(800)
    back = T.KD.params[param_name]
    for val in param_vals :
        T.paramtree.change_dic(param_name, val * back,
                               write = False, verbose = False)
        reduce_params(T.paramtree, MEASURES )
        T2 = Metaphase(800, T.paramtree)
        T2.simul()
        obs = carac(T2.KD, *args)
        observations.append(obs)
    elt = (time() - t0)/float(num_ech)
    print "mean elapsed time by run :: %2.5f seconds " %elt
    param_vals *= back
    param_vals = array(param_vals)
    observations = array(observations)
    if len(observations.shape) > 1: 
        ndim = observations.shape[1]
        if logscale :
            for i in range(ndim):
                figure()
                semilogx(param_vals, observations[:, i], 'o')
        else:
            for i in range(ndim):
                figure()
                plot(param_vals, observations[:, i], 'o')
    else :
        if logscale :
            semilogx(param_vals, observations, 'o')
        else :
            plot(param_vals, observations, 'o')
            
    show()
    return (param_vals, observations)

def spaced(span, logscale, num_ech):
    
    if logscale :
        inf = log10(span[0])
        sup = log10(span[1])
        param_vals = logspace(inf, sup, num = num_ech)
    else :
        inf = span[0]
        sup = span[1]
        param_vals = linspace(inf, sup, num_ech)

    return param_vals

def check_2d(carac, param1_name, span1, param2_name, span2,
             num_ech1 = 50, logscale1 = False,
             num_ech2 = 50, logscale2  = False, *args):

    ''' Same as chek_carac, but displays a 2D colorplot of the choosen
    carac for the two parameters
    '''

    param1_vals = spaced(span1, logscale1, num_ech1)
    param2_vals = spaced(span2, logscale2, num_ech2)

    lx = len(param1_vals)
    ly = len(param2_vals)


    T = Metaphase(800)
    back1 = T.KD.params[param1_name]
    back2 = T.KD.params[param2_name]

    for val1 in param1_vals :
        for val2 in param2_vals :
            T = Metaphase(800)
            T.KD.params[param1_name] = val1 * back1
            T.KD.params[param2_name] = val2 * back2
            T.simul()
            obs = carac(T.KD, *args)
           # if obs > 5:
            observations.append(obs)

    T.KD.params[param1_name] = back1
    T.KD.params[param2_name] = back2

    param1_vals *= back1
    param2_vals *= back2
    observations = array(observations)
    # How to scale axes properly -considering logscales and all???
    # To do: look at imshow
    if len(observations.shape) > 1 : 
        ndim = observations.shape[1]
        for i in range(ndim):
            obs_i = reshape(observations[:,i], (lx,ly))
            figure(i)
            pcolor(param1_vals, param2_vals, obs_i)
    else :
        pcolor(param1_vals, param2_vals, observations)
        
    show()
    return (param1_vals, param2_vals, observations)
    
def check_fO( span, num_ech = 50, num_rep = 1, logscale = False):
    
    observations = []

    param_vals = spaced(span, logscale, num_ech)    
    t0 = time()
    for val in param_vals :
        samples = []
        for k in range(num_rep):
            T = Metaphase(800)
            back_fa = T.KD.params['fa']
            back_fd = T.KD.params['fd']
            T.KD.params['fa'] = val * back_fa
            #T.KD.params['fd'] = val * back_fd
            T.simul()
            obs = auto_corel(T.KD)
            samples.append(obs)
            T.KD.params['fa'] = back_fa
            #T.KD.params['fd'] = back_fd
        observations.append(array(samples))

    elt = (time() - t0)/float(num_ech)
    print "mean elapsed time by run :: %2.5f seconds " %elt
    param_vals *= back_fa
    param_vals = array(param_vals)
    observations = array(observations)
    if len(observations.shape) > 1: 
        ndim = observations.shape[-1]
        if logscale :
            for i in range(ndim):
                figure()
                semilogx(param_vals, observations[..., i], 'o')
        else:
            for i in range(ndim):
                figure()
                plot(param_vals, observations[..., i], 'o')
    else :
        if logscale :
            semilogx(param_vals, observations, 'o', alpha = 1./num_rep)
        else :
            plot(param_vals, observations, 'o', alpha = 1./num_rep)
            
    show()
    return (param_vals, observations)




    


