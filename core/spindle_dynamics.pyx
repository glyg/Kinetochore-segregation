#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module defines all the objects considered in the simulation,
It uses cython used for the computer intensive bits

* Chromosome, Spb, Spindle, PlugSite are the spindle components
* KinetoDynamics wraps all the simulation internals: forces, vectors,
  and the system of equations.
"""

import random
import numpy as np
from scipy import sparse

cimport numpy as np
cimport cython

__all__ = ["KinetoDynamics", "Spb", "Chromosome", "PlugSite", "Spindle"]


class KinetoDynamics(object) :
    """ This class wraps all the simulation internals.

    Public methods:
    ---------------

    Public attributes:
    ------------------

    
    """
    
    def __init__(self, parameters, plug = None):
        """KinetoDynamics instenciation method

        Paramters:
        ----------
            parameters: dict
                A dictionnary of parameters as obtained from a
                xml_handler.ParamTree instance
            plug: string or None
                Defines globally the initial attachment states.
                This argument can have the following values: 
                *'null': all kinetochores are detached
                *'amphitelic': all chromosmes are amphitelic
                *'random': all attachement site can be bound to
                    either pole or deteched with equal prob.
                *'monotelic': right kinetochores are attached
                    to the same pole, left ones are detached
                *'syntelic' : all kinetochores are attached to the same pole
        """
        cdef int N, Mk
        self.params = parameters
        L0 = self.params['L0']
        N = int(self.params['N'])
        Mk = int(self.params['Mk'])
        d0 = self.params['d0']
        k_a = self.params['k_a']
        k_d0 = self.params['k_d0']

        self.spbR = Spb(1, L0) #right spb
        self.spbL = Spb(-1, L0) #left one
        self.spindle = Spindle(self.spbR, self.spbL)
        for n in range(N):
            ch = Chromosome(n, self.spindle, L0, Mk, d0, plug=plug)
            self.chromosomes[n] = ch
        self._idx = self._calc_idx()
        self.B_mat = self.write_B()
        k_a = self.params['k_a']
        self.Pat =  1 - np.exp(-k_a)
        
    def _calc_idx(self):
        """Returns the  index dictionnary
        """
        cdef int Mk, N
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        idxs = {}
        for n in range(N):
            idxs[(0,n)] =  2*n*(Mk+1) + 1
            idxs[(1,n)] = 2*n*(Mk+1) + 2 + Mk
            for m in range(Mk):
                idxs[(0,n,m)] = 2*n*(Mk + 1) + 2 + m
                idxs[(1,n,m)] =  2*n*(Mk + 1) + 3 + Mk + m
        return idxs

    def write_B(self):
        """Returns the matrix containing the linear terms
        of the equation set A\dot{X} + BX + C = 0
        """
        cdef N, Mk
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        kappa_c = self.params['kappa_c']
        kappa_k = self.params['kappa_k']
        dim = 1 + N * ( 1 + Mk ) * 2
        B = np.zeros((dim, dim), dtype = float)
        B[0,0] = 0
        for n in range(N) :
            B[self._idx[(0, n)],
              self._idx[(0, n)]] = -kappa_c - Mk*kappa_k
            B[self._idx[(1, n)],
              self._idx[(1, n)]] = -kappa_c - Mk*kappa_k
            B[self._idx[(1, n)],
              self._idx[(0, n)]] = kappa_c
            B[self._idx[(0, n)],
              self._idx[(1, n)]] = kappa_c
            for m in range(Mk):
                B[self._idx[(0, n, m)],
                  self._idx[(0, n, m)]] = - kappa_k
                B[self._idx[(1, n, m)],
                  self._idx[(1, n, m)]] = - kappa_k
                B[self._idx[(0, n)],
                  self._idx[(0, n, m)]] = kappa_k
                B[self._idx[(0, n, m)],
                  self._idx[(0, n)]] = kappa_k
                B[self._idx[(1, n)],
                  self._idx[(1, n, m)]] = kappa_k
                B[self._idx[(1, n, m)],
                  self._idx[(1, n)]] = kappa_k
        return B
    
    def get_state_vector(self):
        """
        return a vector of the positions of each components
        """
        N = int(self.params['N'])
        Mk = int(self.params['Mk'])
        X = np.zeros(1 + 2*N * ( Mk + 1 ))
        X[0] = self.spbR.pos
        for n in range(N):
            ch = self.chromosomes[n]
            X[self._idx[(0, n)]] = ch.rightpos
            X[self._idx[(1, n)]] = ch.leftpos        
            for m in range(Mk):
                X[self._idx[(0, n, m)]] = ch.rplugs[m].pos
                X[self._idx[(1, n, m)]] = ch.lplugs[m].pos
        return X
        
    #Wrighting the equations (see description.pdf for details)
    ##### --------------- METAPHASE -------------------------#######
    ### Spindle pole terms

    #@cython.profile(False)
    def calc_ldep(self, double plugsite_pos, double pole_pos):
        """Calculates the length dependency of the force applied
        at the plugsite.

        Parameters:
        -----------
        plugsite_pos: double
            The position of the plugsite
        pole_pose: double 
            The position of the attached spindle pole

        If the KD.param ld_slope is null, there is no length
        dependence, and 1. is returned


        Returns:
        --------
        ldep: double
            The prefactor to the force term
        """
        cdef double ld0, ld_slope
        ld_slope = self.params['ld_slope']
        if ld_slope == 0:
            return 1.
        ld0 = self.params['ld0']
        cdef double dist
        dist = abs(pole_pos - plugsite_pos)
        cdef double ldep
        ldep = ld_slope * dist + ld0
        # TODO: investigate: introduces a first order discontinuity
        #     suspected to trigger artifacts when the plugsite
        # is close to the pole
        if dist < 0.0001:
            ldep = dist  # No force when at pole
        return ldep
    
    def alpha(self, ldep, int p):
        """Returns 1 if well plugged
        length dependance implementation:
        alpha depends linearly on the spb-kt distance
        """
        
        
        return ldep * p * (1 + p) / 2

    #@cython.profile(False)
    def beta(self, double pspos, double spos, int p):
        """returns 1 if merotelic
        length dependance implementation:
        beta depends linearly on the spb-kt distance
        
        if ld_slope = 0 => No length dependance
        """
        cdef double ld0, ld_slope
        ld_slope = self.params['ld_slope']
        ld0 = self.params['ld0']
        cdef double dist
        dist = abs(spos - pspos)
        cdef double ldep
        ldep = ld_slope * dist + ld0
        if dist < 0.0001:
            ldep = dist  # No force when at pole
        return  ldep * p * (p - 1) / 2

    def calcA(self):
        
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef double muc = self.params['muc']
        cdef float muk = self.params['muk']
        cdef float mus = self.params['mus']
        cdef float Vmz = self.params['Vmz']
        cdef float Fmz = self.params['Fmz']

        #create a well- sized array
        cdef int dims
        dims = 1 + 2*N * ( Mk + 1 ) 
        #cdef np.ndarray[np.float64_t, ndim =2] A 
        A = np.zeros((dims, dims))
        #A = sparse.lil_matrix((dims, dims))
        
        cdef float sposL = self.spbR.pos
        cdef float sposL = self.spbL.pos
        cdef float pspos
        cdef int p, pR, pL

        A[0,0] = - 2 * mus - 4 * Fmz / Vmz

        # inner plate
        for n in range(N):
            ch = self.chromosomes[n]
            A[self._idx[(0,n)],
              self._idx[(0,n)]] = - Mk * muk - muc #self.am_nn() #right
            A[self._idx[(1,n)],
              self._idx[(1,n)]] = - Mk * muk - muc #self.am_nn() # left

            #outer plate
            for m in range(Mk):
                
                psposR = ch.rplugs[m].pos
                pR = ch.rplugs[m].plug
                alphaR = self.alpha(psposR, sposR, pR)
                betaR = self.beta(psposR, sposL, pR)

                psposL = ch.lplugs[m].pos
                pL = ch.lplugs[m].plug
                alphaL = self.alpha(psposL, sposL, pL)
                betaL = self.beta(psposL, sposR, pL)
                
                #spbs diag terms:
                A[0,0] += -alphaR - betaR
                A[0,0] += -alphaL - betaL

                #Right side
                pspos = ch.rplugs[m].pos
                p = ch.rplugs[m].plug

                A[self._idx[(0,n,m)],
                  self._idx[(0,n,m)]] = -muk - alphaR - betaR
                A[0, self._idx[(0,n,m)]] = alphaR - betaR
                A[self._idx[(0,n,m)], 0] = A[0, self._idx[(0,n,m)]]
                A[self._idx[(0,n,m)],
                  self._idx[(0,n)]] = muk
                A[self._idx[(0,n)],
                  self._idx[(0,n,m)]] = muk

                #Left side
                pspos = ch.lplugs[m].pos
                p = ch.lplugs[m].plug

                A[self._idx[(1,n,m)],
                  self._idx[(1,n,m)]] = -muk - alphaL - betaL
                A[0, self._idx[(1,n,m)]] = - alphaL + betaL
                A[self._idx[(1,n,m)], 0] = A[0, self._idx[(1,n,m)]]
                A[self._idx[(1,n,m)],
                  self._idx[(1,n)]] = muk
                A[self._idx[(1,n)],
                  self._idx[(1,n,m)]] = muk
            
        return A#.tocsr()

        # Constant vector
    #@cython.profile(False)
    def cm_S(self):
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef float Fmz = self.params['Fmz']
        cdef float sposR = self.spbR.pos
        cdef float sposL = self.spbL.pos
        cdef float psposR, psposL
        cdef int pR, pL

        cdef float b = 2 * Fmz
        for n in range(N):
            ch = self.chromosomes[n]
            for m in range(Mk):
                psposR = ch.rplugs[m].pos
                pR = ch.rplugs[m].plug
                psposL = ch.lplugs[m].pos
                pL = ch.lplugs[m].plug

                b += - self.alpha(psposR, sposR, pR) - self.beta(psposR, sposL, pR)
                b += - self.alpha(psposL, sposL, pL) - self.beta(psposL, sposR, pL)

        return b

    #@cython.profile(False)
    def cm_n(self, int side) :
        """ kineto """
        cdef float d0 = self.params['d0']
        cdef float kappa_c = self.params['kappa_c']
        return (1 - 2*side) * kappa_c * d0

    #@cython.profile(False)
    def cm_m(self, int side, int p, float pspos,
              float spos_good, float spos_bad):
        cdef float b
        b = (1 - side * 2) * ( self.alpha(pspos, spos_good, p)
                               - self.beta(pspos, spos_bad, p) )
        return b


    def calcc(self):

        cdef int Mk = int(self.params['Mk'])
        cdef int N = int(self.params['N'])
        
        C = np.zeros((1 + N * ( 1 + Mk ) * 2))

        C[0] = self.cm_S()
        cdef float sposR = self.spbR.pos
        cdef float sposL = self.spbL.pos
        cdef float pspos
        cdef int p

        for n in range(N) :
            C[self._idx[(0,n)]] = self.cm_n(0)
            C[self._idx[(1,n)]] = self.cm_n(1)
            ch = self.chromosomes[n]
            for m in range(Mk):
                pspos = ch.rplugs[m].pos
                p = ch.rplugs[m].plug
                C[self._idx[(0,n,m)]] = self.cm_m(0, p, pspos, sposR, sposL)

                pspos = ch.lplugs[m].pos
                p = ch.lplugs[m].plug
                C[self._idx[(1,n,m)]] = self.cm_m(1, p, pspos, sposL, sposR)

        return C
        

    def calcb(self):
        
        B = self.B_mat
        C = self.calcc()
        X = self.get_state_vector()

        ### --- THIS TRIGGERS BAD BEHAVIOUR - INVESTIGATE LATER 
        # ## Correction for the modeling  of chromatin: If
        # ## The distance between two kts is less than the equilibrium distance
        # ## The force should be 0 (it shall not oppose the congression
        # ## As seen experimentaly
        # N = self.params['N']
        # d0 = self.params['d0']
        # # for n in range(N):
        # #     ch = self.chromosomes[n]
        # #     if abs(ch.rightpos - ch.leftpos) < d0 :
        # #         B[self._idx[(0,n)], self._idx[(1,n)]] = 0
        # #         B[self._idx[(1,n)], self._idx[(0,n)]] = 0

        #product = np.dot(B.todense(), X)
        product = np.dot(B,X)
        pos_dep = product + C
        return pos_dep

    def test_anaphase_switch(self):
        N = int(self.params["N"])
        for n in range(N):
            ch = self.chromosomes[n]
            if ch.anaphase_switch[0] == 0:
                if ch.isatrightpole() and (ch.mero()[0] == 0):
                    #print 'switch!'
                    ch.anaphase_switch[0] = 1
            if  ch.anaphase_switch[1] == 0:                    
                if ch.isatleftpole() and (ch.mero()[1] == 0):
                    ch.anaphase_switch[1] = 1


    #@cython.profile(False)
    # def  Pat(self):#,  n, m, side):
    #     """Returns attachement probability for kineto n
    #     """

    #     cdef double k_a = self.params['k_a']
    #     # cdef double dt = self.params['dt']
    #     #k_at = fa
    #     return 1 - np.exp( - k_a ) # * dt <- unnecessary due to non-dimentionalization


    #@cython.profile(False)
    def Pdet(self, int n, double plugpos):
        """Calculates detachement frequency for kineto n
        side = 0 : right
        side = 1 : left
        """
        
        cdef int Mk = int(self.params['Mk'])
        cdef double k_d0 = self.params['k_d0']
        cdef double d_alpha = self.params['d_alpha']
        # cdef double dt = self.params['dt']

        ch = self.chromosomes[n]
        dist = ch.plug_dist(plugpos)
        
        ### Aurora ???
        cdef double k_dc
        if d_alpha != 0 :
            k_dc = k_d0  *  d_alpha / dist #np.exp( - dist / d_alpha)
            if k_dc > 1e4 : #* dt > 1e4:
                return 1.
        else:
            k_dc = k_d0
        return 1 - np.exp(-k_dc) 


    def Pmero(self, int n, int side):
        """
        
        """

        cdef int Mk
        cdef double p
        cdef double m
        cdef double orientation
        Mk = int(self.params['Mk'])
        ch = self.chromosomes[n]
        p = float(ch.pluged()[side]) #        p = np.sum(ch.pluged(), dtype=float) #
        m = float(ch.mero()[side]) #        m = np.sum(ch.mero, dtype=float)# 
        orientation = self.params['orientation']

        if orientation == 0:
            return 0.5

        if p + m  == 0:
            return 0.5
        # if m > p :
        #     return 0.5
        Pmero = 0.5 + orientation * (m - p) / (2 * (m + p))

        return Pmero
        
        
    def plug_unplug(self):
        """Let's play dices ...
        """
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])

        for n in range(N):
            ch = self.chromosomes[n]
            (right_pluged, left_pluged) = ch.pluged()
            (right_mero, left_mero) = ch.mero()
            for m in range(Mk):
                # Attached kMTs have a chance to unplug (until anaphaseB):
                if ch.rplugs[m].plug == 1:# and ch.anaphase_switch[0] == 0:
                    dice = random.random()
                    if  dice < self.Pdet(n, ch.rplugs[m].pos):
                        ch.rplugs[m].plug = 0
                elif ch.rplugs[m].plug == -1 :
                    dice = random.random()
                    if  dice < self.Pdet(n, ch.rplugs[m].pos):
                        ch.rplugs[m].plug = 0
                # Unattached kMTs have a chance to plug:                    
                else :
                    dice = random.random()
                    if  dice < self.Pat:#n, m, 0):
                        #Implementing the possibility to attach mero kts:
                        m_dice = random.random()
                        if m_dice < self.Pmero(n, 0):
                            ch.rplugs[m].plug = -1
                        else:
                            ch.rplugs[m].plug = 1

                if ch.lplugs[m].plug == 1 :# and ch.anaphase_switch[1] == 0:
                    dice = random.random()
                    if  dice < self.Pdet(n, ch.lplugs[m].pos):
                        ch.lplugs[m].plug = 0
                elif ch.lplugs[m].plug == -1:
                    dice = random.random()
                    if  dice < self.Pdet(n, ch.lplugs[m].pos):
                        ch.lplugs[m].plug = 0
                # Unattached kMTs have a chance to plug:                    
                else :
                    dice = random.random()
                    if  dice < self.Pat:
                        #Implementing the possibility to attach mero kts:
                        m_dice = random.random()
                        if m_dice < self.Pmero(n, 1):
                            ch.lplugs[m].plug = -1
                        else:
                            ch.lplugs[m].plug = 1
    
            #update
            # ch.pluged() = (right_pluged, left_pluged)
            # ch.mero = (right_mero, left_mero)

            #swap
            if np.sum(ch.pluged()) < np.sum(ch.mero()):
                ch.swap()


    ### So now, let's update the positions

    def position_update(self, speeds):
        """given the speeds obtained by solving Atot.x = btot
        and caclulated switch events 
        """
        cdef int Mk = int(self.params['Mk'])
        cdef int N = int(self.params['N'])
        cdef double dt = self.params['dt']
        cdef double Vk = self.params['Vk']
        
        speeds *= Vk * dt #Back to real space
        self.spbR.pos += speeds[0]
        self.spbL.pos -= speeds[0]
        if self.spbR.pos <= self.spbL.pos + 0.2:
            print "Crossing "
            self.spbR.pos = 0.1
            self.spbL.pos = - 0.1

        self.spbL.traj.append(self.spbL.pos)        
        self.spbR.traj.append(self.spbR.pos)

        for n in range(N):
            ch = self.chromosomes[n]
            ch.rightpos += speeds[self._idx[(0,n)]]
            if ch.rightpos > self.spbR.pos:
                ch.rightpos = self.spbR.pos

            ch.leftpos +=  speeds[self._idx[(1,n)]] 
            if ch.leftpos < self.spbL.pos:
                ch.leftpos = self.spbL.pos

            # if ch.dist() < 0: # undo the displacement!
            #     ch.rightpos -= dt * speeds[self._idx[(0,n)]]
            #     ch.leftpos -= dt * speeds[self._idx[(1,n)]]

            for m in range(Mk):
                ch.rplugs[m].pos +=  speeds[self._idx[(0,n,m)]]
                if self.spbL.pos > ch.rplugs[m].pos:
                   ch.rplugs[m].pos = self.spbL.pos
                if self.spbR.pos < ch.rplugs[m].pos:
                   ch.rplugs[m].pos = self.spbR.pos
                ch.rplugs[m].traj.append(ch.rplugs[m].pos)
                    
                ch.lplugs[m].pos +=  speeds[self._idx[(1,n,m)]]                
                if self.spbL.pos > ch.lplugs[m].pos:
                   ch.lplugs[m].pos = self.spbL.pos
                if self.spbR.pos < ch.lplugs[m].pos:
                   ch.lplugs[m].pos = self.spbR.pos
                ch.lplugs[m].traj.append(ch.lplugs[m].pos)

                ch.rplugs[m].state_hist.append(ch.rplugs[m].plug)
                ch.lplugs[m].state_hist.append(ch.lplugs[m].plug)                


            ch.pluged_history.append(ch.pluged())
            ch.mero_history.append(ch.mero())

            ch.lefttraj.append(ch.leftpos)
            ch.righttraj.append(ch.rightpos)

    def get_sim_duration(self):
        return (len(self.spbR.traj) - 1)*self.params['dt']


    #### ---- Not so useful code after all


    def cohesin_bound(self,n) :
        """ returns the spring force due to cohesin binding
        between siter chromatids """
        kappa_c = self.params['kappa_c']
        d0 = self.params['d0']
        ch = self.chromosomes[n]
        dist = ch.dist()
        F = - kappa_c * ( dist - d0)
        return F

    def chromatid_bound(self, n, m, side):
        kappa_k = self.params['kappa_k']

        ch = self.chromosomes[n]
        xnm = self.ppos(n,m,side)
        if side == 0 :
            xn = ch.rightpos
        else:
            xn = ch.leftpos

        F = - kappa_k * (xnm - xn) # null equilibrium distance
        return F

        
    def ppos(self, n, m, side):
        ch = self.chromosomes[n]
        if side == 0:
            kt = ch.rplugs[m]
        else:
            kt = ch.lplugs[m]
        return kt.pos

class Chromosome(object):
    """Chromosome Object """

    def __init__(self, i, spindle, L0, Mk,  d0 = 0.5, doi = 0.01, plug = None):

        """ Degrees of freedom: Position of iner plate + positions of each attachmnt site
        d0 [µm] kinetochore - kinetochore rest length
        """

        self.index = i
        self.spindle = spindle 
        center = random.gauss(0, 0.2 * (L0 - d0))
        self.leftpos = center - d0 / 2
        self.leftspeed = 0

        self.rightpos = center +  d0 / 2
        self.rightspeed = 0
        self.right_toa = 0 #Time of arrival
        self.left_toa = 0 #Time of arrival
        
        self.righttraj = [self.rightpos]
        self.lefttraj = [self.leftpos]

        #Random attachment at prometaphase : each kMT plug site
        #can be in one of three states: good (1), mero(-1), null(0)
        (nd,ng) = (0,0)
        (md,mg) = (0,0)
        self.rplugs = []
        self.lplugs = []
        
        for m in range(Mk):
            rps = PlugSite(ch = self, side = 0, doi =0.01, plug = plug)
            nd += rps.plug * (1 + rps.plug) / 2 # = 1 iff rplug = 1
            md += rps.plug * (rps.plug - 1) / 2 # = 1 iff rplug = -1
            self.rplugs.append((m, rps))
            
            lps = PlugSite(ch = self, side = 1, doi =0.01, plug = plug)
            ng += lps.plug * (1 + lps.plug) / 2 # = 1 iff rplug = 1
            mg += lps.plug * (lps.plug - 1) / 2 # = 1 iff lplug = -1
            self.lplugs.append((m, lps))

        self.rplugs = dict(self.rplugs)
        self.lplugs = dict(self.lplugs)
        
        self.pluged_history = [self.pluged()]

        self.mero_history = [self.mero()]
        self.anaphase_switch = [0,0]# = 1 once the kineto reached the pole at anaphase
        self.active_sac = 0

    def pluged(self):
        cdef int rpluged 
        cdef int lpluged 
        rpluged, lpluged = 0, 0
        for rps, lps in zip(self.rplugs.values(), self.lplugs.values()) :
            if rps.plug > 0:
                rpluged += 1
            if lps.plug > 0:
                lpluged += 1
        return rpluged, lpluged

    def mero(self):
        cdef int rmero
        cdef int lmero
        rmero, lmero = 0,0
        for rps, lps in zip(self.rplugs.values(), self.lplugs.values()) :
            if rps.plug < 0:
                rmero += 1
            if lps.plug < 0:
                lmero += 1
        return rmero, lmero

        
    def dist(self):
        return (self.rightpos - self.leftpos)

    def plug_dist(self, double plugpos):
        return abs(self.center() - plugpos)

    def center(self):
        return (self.rightpos + self.leftpos)/2

    def cen_traj(self):
        return (np.array(self.righttraj) + np.array(self.lefttraj))/2
        

    def isatrightpole(self, double tol = 0.0) :
        """tol : tolerance distance
        """
        if self.rightpos >= self.spindle.spbR.pos - tol:
            return True
        else :
            return False

    def isatleftpole (self, double tol = 0.0) :
        if abs(self.leftpos - self.spindle.spbL.pos) <= tol:
            return True
        else :
            return False
            
    def isatpole(self, side = None, double tol = 0.0):

        if side == 0 and self.isatrightpole(tol):
            return True
        elif side == 1 and self.isatleftpole(tol):
            return True
        elif side is None:
            if self.isatrightpole(tol) or self.isatleftpole(tol):
                return True
        else :
            return False


    def swap(self):
        """
        change the sides of the kinetochores
        """
        self.leftpos, self.rightpos = self.rightpos, self.leftpos
        self.lefttraj, self.righttraj = self.righttraj, self.lefttraj
        self.leftspeed, self.rightspeed =  self.rightspeed, self.leftspeed
        self.pluged_history, self.mero_history = self.mero_history, self.pluged_history

        for rps, lps in zip(self.rplugs.values(), self.lplugs.values()):
            rps.plug, lps.plug = -rps.plug, -lps.plug
            rps.state_hist, lps.state_hist = lps.state_hist, rps.state_hist
            rps.pos, lps.pos = lps.pos, rps.pos
            rps.traj, lps.traj = lps.traj, rps.traj

class PlugSite(object):
    """
    An attachment site object.

    Instanciation: PlugSite(ch, side, doi, plug).
    The plug argument can be either of the following:
    None (the default), in which case the PlugSite is randomly attached
    'null' -- -- self.plug = 0
    'amphitelic' -- self.plug = 1
    'random' -- self.plug = -1, 0, or 1 with equal probability
    'monotelic' -- self.plug = 1 for right side PlugSite and self.plug = 0 for left side ones
    'syntelic' --  self.plug = 1 for right side PlugSite and self.plug = -1 for left side ones
    
    """


    def __init__(self, ch, side, doi, plug = None):

        if side == 0:
            self.pos = ch.rightpos + doi
        else:
            self.pos = ch.leftpos - doi
        if plug == None: 
            self.plug = random.randint(-1, 1)
        elif plug == 'null':
            self.plug = 0
        elif plug == 'amphitelic':
            self.plug = 1
        elif plug == 'random':
            self.plug = random.randint(-1, 1)
        elif plug == 'monotelic':
            if side == 0:
                self.plug = 1
            else:
                self.plug = 0
        elif plug == 'syntelic':
            if side == 0:
                self.plug = 1
            else:
                self.plug = -1
        elif plug == 'merotelic':
            self.plug = random.choice([-1,1])
        else:
            self.plug = plug
        self.traj = [self.pos]
        self.state_hist = [self.plug]
        


class Spb(object) :

    """ A spindle pole object.

    Attributes:
    side: 1 or -1 (right or left. In the article, -1 correspond to the
    daggered variable)
    """
    
    def __init__(self, side, L0):
        """ side = 1 : left
        side = -1 : right
        L0 : spindle length
        """
        self.side = side
        self.pos = side * L0 / 2
        self.traj = [self.pos]
        
class Spindle(object) :

    def __init__(self, spbR, spbL ):

        self.spbR = spbR
        self.spbL = spbL
        
    def length(self):

        return self.spbR.pos - self.spbL.pos

    def length_traj(self):

        return np.array(self.spbR.traj) - np.array(self.spbL.traj)

    
