#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:Guillaume Gay<elagachado AT  gmail DOT com>
## Commentary:

### MT dynamics characteristics from Sagolla, Uzawa and Cande J. Cell Sc. 2003
vg = 0.01
vs = 0.08
dt = 1.



from numpy import array, diff, r_
from pylab import find

def infer_traj(plugsite):

    sh = array(plugsite.state_hist)
    traj = array(plugsite.traj)
    
    num_points = traj.shape[0]
    unpluged = find(sh==0)
    change_points = unpluged[diff(unpluged)>1]

    if sh[0] == 0 : #unplugged start
        ti = 0
        if sh[-1] == 0: #unpluged_end
            for tf, next_ti in zip( change_points[:-1:2],
                               change_points[1::2] ):
                free_MT_traj(traj, ti, tf) 
                ti = next_ti
            traj[ti:] = traj[ti] - vs * ( r_[ti:num_points] - ti )

        else: #pluged_end
            for tf, next_ti in zip( change_points[:-2:2],
                                    change_points[1:-1:2] ):
                free_MT_traj(traj, ti, tf) 
                ti = next_ti

            tf = change_points[-1]
            free_MT_traj(traj, ti, tf) 

    else: #pluged start

        if sh[-1] == 0: #unpluged_end
            for ti, tf in zip( change_points[:-2:2],
                               change_points[1:-1:2] ):
                free_MT_traj(traj, ti, tf) 
            traj[ti:] = traj[ti] - vs * ( r_[ti:num_points] - ti )

        else : #pluged_end
            for ti, tf in zip( change_points[:-1:2],
                               change_points[1::2] ):
                free_MT_traj(traj, ti, tf) 


    return traj
            

def free_MT_traj(traj, ti, tf) :

    xi = traj[ti]
    xf = traj[tf]
    tp = ((xf - xi) - (tf -ti) * vg) / (vg + vs) + ti
    traj[ti:tp] = xi - vs * ( r_[ti:tp] - ti )
    traj[tp:tf] = xf + vg * ( r_[tp:tf] - tf )

            

            
        
        
        
    
