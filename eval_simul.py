#!/usr/bin/python
# -*- coding: utf-8 -*-
## imports

import os
import sys
from numpy import correlate, sqrt, diff, r_, hanning, convolve
from scipy import array, arange, polyfit, argmax
#from scipy.stats import mean, std
from numpy import mean, std
from scipy.interpolate import splrep, splev
from numpy.fft import rfft, fftfreq
from numpy.random import normal
import pyximport
pyximport.install()

#from pylab import *
#

from glygs_utils.array_utils import first_min, fwhm


## local imports
from kt_simul.xml_handler import *
from kt_simul.spindle_dynamics import *


__all__ = ["anaphase_rate", "metaph_rate",
           "metaph_kineto_dist", "auto_corel", "poleward_speed",
           "kt_rms_speed", "time_of_arrival", "pluged_stats"]
           


def get_kt_speeds(KD, step = 2):
    
    speeds = []
    for ch in KD.chromosomes.values():

        speeds.append(diff(array(ch.righttraj[::step]))/step)
        speeds.append(diff(array(ch.lefttraj[::step]))/step)

    return speeds
    
def get_spb_speed(KD, step = 2):
    
    speed = diff(array(KD.spbR.traj[::step]))/step
    return speed


def anaph_transAB(KD):
    '''
    As anaphase A -> B transition is defined (at least in the lab)
    by the date at which the first kinetochore reaches and stay at the spb
    we have to retrieve this moment. Fortunately, we already know times of arrival
    '''
    toas = []
    for ch in KD.chromosomes.values():
        toas.append(ch.right_toa)
        toas.append(ch.left_toa)
    trans_AB = min(toas)
    return trans_AB
    
def anaphase_rate(KD):
    '''Calculates anaphase B elongation rate.
    '''
    t_A = KD.params['t_A']
    dt = KD.params['dt']
    trans_AB = anaph_transAB(KD)
    stop = len(KD.spbR.traj)*dt
    if not trans_AB:
        return 0
    else:
        spindle_length = array(KD.spbR.traj) - array(KD.spbL.traj)
        elapsed = arange(trans_AB, stop, dt)
        (a, b)=polyfit(elapsed, spindle_length[int(trans_AB/dt):], 1)

    return a

def metaph_rate(KD):
    '''Calculates metaphase elongation rate.
    '''
    N = int(KD.params['N'])
    t_A = KD.params['t_A']
    dt = KD.params['dt']
    spindle_length = array(KD.spbR.traj) - array(KD.spbL.traj)
    if spindle_length.size >= int(t_A/dt):
        elapsed = arange(t_A, step = dt)
        (a, b)=polyfit(elapsed, spindle_length[:elapsed.shape[0]], 1)
        print dt

    else:
        elapsed = arange(spindle_length.size*dt, step = dt)
        (a, b)=polyfit(elapsed, spindle_length, 1)
        
    return a

    
def metaph_kineto_dist(KD):
    '''calculates the mean and sdev of the distance between the kinetochores
    of each chromosome during metaphase. Returns a tuple (mean, sdev)
    '''
    N = int(KD.params['N'])
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    trans_MA = int(trans_MA/dt)
    average = 0
    stdev = 0
    for n in range(N):
        kinetoR = array(KD.chromosomes[n].righttraj)
        kinetoL = array(KD.chromosomes[n].lefttraj)

        dist = kinetoR - kinetoL
        dist = dist[:trans_MA]

        average += mean(dist)
        stdev += std(dist)


    average /= N
    stdev /= N

    return average, stdev

def add_noise(KD, detect_noise = 65e-3, smth = 5):

    wd = hanning(smth)
    KD.spbR.traj = array(KD.spbR.traj)
    
    KD.spbR.traj += normal(scale = detect_noise, size = KD.spbR.traj.shape)
    KD.spbR.traj = convolve(KD.spbR.traj, wd/sum(wd), mode = 'same')
    for ch in KD.chromosomes.values():
        ch.righttraj = array(ch.righttraj)
        ch.righttraj += normal(scale = detect_noise, size = ch.righttraj.shape)
        ch.righttraj = convolve(ch.righttraj, wd/sum(wd), mode = 'same')
        ch.lefttraj = array(ch.lefttraj)
        ch.lefttraj += normal(scale = detect_noise, size = ch.lefttraj.shape)
        ch.lefttraj = convolve(ch.lefttraj, wd/sum(wd), mode = 'same')

def metaph_kineto_dist_list(KD):
    N = int(KD.params['N'])
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    trans_MA = int(trans_MA/dt)
    dist_list = []
    for n in range(N):
        kinetoR = array(KD.chromosomes[n].righttraj)
        kinetoL = array(KD.chromosomes[n].lefttraj)
        dist = kinetoR - kinetoL
        dist = dist[:trans_MA]
        dist_list.append(dist)

    return dist_list

def poleward_speed(KD):

    N = int(KD.params['N'])
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    start = int(trans_MA/dt)
    
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    trans_AB = anaph_transAB(KD)

    pole_speeds = []
    for ch in KD.chromosomes.values():
        r_stop = int(ch.right_toa/dt)
        if r_stop - start > 2:
            r_dist = - array(KD.spbR.traj)[start:r_stop] + array(ch.righttraj)[start:r_stop]
            elapsed = r_[trans_MA:ch.right_toa:dt]
            (ra,rb) = polyfit(elapsed, r_dist, 1)
            pole_speeds.append(ra)
        l_stop = int(ch.left_toa/dt)
        if l_stop - start > 2:
            l_dist = array(KD.spbL.traj)[start:l_stop] - array(ch.lefttraj)[start:l_stop]
            elapsed = r_[trans_MA:ch.left_toa:dt]
            (la,lb) = polyfit(elapsed, l_dist, 1)
            pole_speeds.append(la)
            
    pole_speeds = array(pole_speeds)

    return pole_speeds.mean(), pole_speeds.std()

def kt_rms_speed(KD):

    N = int(KD.params['N'])
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    stop = int(trans_MA/dt)
    rms_speeds = []

    for ch in KD.chromosomes.values():
        r_speed = diff(array(ch.righttraj)[:stop])
        r_rmss = sqrt((r_speed**2).mean())
        rms_speeds.append(r_rmss)

        l_speed = diff(array(ch.lefttraj)[:stop])
        l_rmss = sqrt((l_speed**2).mean())
        rms_speeds.append(l_rmss)

    rms_speeds = array(rms_speeds)

    return rms_speeds.mean(), rms_speeds.std()
    
def time_of_arrival(KD):

    toas = []
    for ch in KD.chromosomes.values():

        toas.append(ch.right_toa)
        toas.append(ch.left_toa)

    toas = array(toas)
    if any(toas):
        toas -= toas[toas>0].min()
    else :
        toas -= 1.
    return toas

def pluged_stats(KD):

    N = int(KD.params['N'])
    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    trans_MA = int(trans_MA/dt)
    tot_avg = 0
    delta_avg = 0
    for ch in KD.chromosomes.values():
        plugedR = array(ch.pluged_history)[:, 0]
        plugedL = array(ch.pluged_history)[:, 1]

        tot = plugedR + plugedL
        tot = tot[:trans_MA]
        delta = abs(plugedR - plugedL)
        delta = delta[:trans_MA]

        tot_avg += mean(tot)/2
        delta_avg += mean(delta)

    tot_avg /= N
    delta_avg /= N

    return (tot_avg, delta_avg)

def auto_corel(KD, smooth = 10.):
    '''
    Calculates the "pitch" of chromosomes"trajectory (kind of a
    Pitch detection algorithm
    '''

    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    pitches = []
    if len(KD.spbR.traj) <= trans_MA/dt: #no anaphase execution
        elapsed = arange(len(KD.spbR.traj)*dt, step = dt)
    else:
        elapsed = arange(0, trans_MA, dt)
    smth = int(smooth/dt)
    
    for ch in KD.chromosomes.values():

        ktR = array(ch.righttraj)[:elapsed.shape[0]]
        # In order to compare with the 'real world' tracked kinetochores
        # we smooth by a factor smooth/dt
        ktR_tck = splrep(elapsed, ktR, t = elapsed[smth:-smth:smth])
        ktR_s = splev(elapsed, ktR_tck, der = 1)
        m_speed, st_speed = ktR_s.mean(), ktR_s.std()
        ktR_sc = (ktR_s - m_speed)/st_speed
        co_ktR = correlate(ktR_sc, ktR_sc, 'full') / ktR_sc.size
        pitches.append(first_min(co_ktR[-co_ktR.size//2:]))

        ktL = array(ch.lefttraj)[:elapsed.shape[0]]
        # In order to compare with the 'real world' tracked kinetochores
        # we smooth by a factor smooth/dt
        ktL_tck = splrep(elapsed, ktL, t = elapsed[smth:-smth:smth])
        ktL_s = splev(elapsed, ktL_tck, der = 1)
        m_speed, st_speed = ktL_s.mean(), ktL_s.std()
        ktL_sc = (ktL_s - m_speed)/st_speed
        co_ktL = correlate(ktL_sc, ktL_sc, 'full') / ktL_sc.size
        pitches.append(first_min(co_ktL[-co_ktL.size//2:]))
    try:    
        pitches = 1/(array(pitches)*dt)
    except:
        return 0
    return pitches.mean(), pitches.std()
    

def max_freqs(KD, show_fig = True):

    trans_MA = KD.params['t_A']
    dt = KD.params['dt']
    smth = int(10./dt)
    elapsed = arange(0, trans_MA , dt)
    if show_fig : figure(100)
    max_fs = []
    for ch in KD.chromosomes.values():
        cen = (array(ch.righttraj) + array(ch.lefttraj))/2
        cen_m = cen[:elapsed.shape[0]]
        knots = elapsed[smth:-smth:smth]
        tck = splrep(elapsed, cen_m, t = knots)
        sp_speed = splev(elapsed, tck, der = 1)
        sp_fft = rfft(sp_speed)
        freqs = fftfreq(elapsed.shape[0], dt)
        m_freq = freqs[argmax(sp_fft)]
        max_fs.append(m_freq)
        if show_fig:
            plot(freqs[:len(sp_fft)], abs(sp_fft), 'o', alpha = 0.6)

    return max_fs


