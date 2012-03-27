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
from scipy import sparse, linalg
cimport numpy as np
cimport cython

__all__ = ["KinetoDynamics", "Spb", "Chromosome",
           "Centromere", "PlugSite", "Spindle"]


class KinetoDynamics(object) :
    """ This class wraps all the simulation internals.
    
    Public methods:
    ---------------
    
    Public attributes:
    ------------------

    
    """
    
    def __init__(self, parameters, initial_plug = None):
        """KinetoDynamics instenciation method

        Paramters:
        ----------
        parameters: dict
            A dictionnary of parameters as obtained from a
            xml_handler.ParamTree instance
        initial_plug: string or None
            Defines globally the initial attachment states.
            This argument can have the following values: 
            *'null': all kinetochores are detached
            *'amphitelic': all chromosmes are amphitelic
            *'random': all attachement site can be bound to
                either pole or detached with equal prob.
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
        duration = self.params['span']
        dt = self.params['dt']
        self.num_steps = int(duration/dt)
        self.spindle = Spindle(self)
        self.spbR = Spb(self.spindle, 1, L0) #right spb
        self.spbL = Spb(self.spindle, -1, L0) #left one
        self.initial_plug = initial_plug
        self.chromosomes = [Chromosome(self.spindle)
                            for n in range(N)]
        self._idx = self._calc_idx()
        self.B_mat = self.write_B()
        self.step_count = 0
        
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

    def solve(self):
        A = self.calcA()
        b = - self.calcb()
        speeds = linalg.solve(A, b)
        return speeds

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
              self._idx[(0, n)]] = - kappa_c - Mk*kappa_k
            B[self._idx[(1, n)],
              self._idx[(1, n)]] = - kappa_c - Mk*kappa_k
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
            X[self._idx[(0, n)]] = ch.cen_A.pos
            X[self._idx[(1, n)]] = ch.cen_B.pos        
            for m in range(Mk):
                X[self._idx[(0, n, m)]] = ch.cen_A.plugsites[m].pos
                X[self._idx[(1, n, m)]] = ch.cen_B.plugsites[m].pos
        return X
        
    #Wrighting the equations (see doc/segregation_model.pdf for details)
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
        dims = 1 + 2 * N * ( Mk + 1 ) 
        A = np.zeros((dims, dims))
        cdef float p_nmA, p_nmB
        cdef float ldep_A, ldep_B
        A[0,0] = - 2 * mus - 4 * Fmz / Vmz
        # inner plate
        for n in range(N):
            ch = self.chromosomes[n]
            A[self._idx[(0,n)],
              self._idx[(0,n)]] = - Mk * muk - muc
            A[self._idx[(1,n)],
              self._idx[(1,n)]] = - Mk * muk - muc
            for m in range(Mk):
                plugsite_A = ch.cen_A.plugsites[m]
                ldep_A = plugsite_A.calc_ldep()
                p_nmA = plugsite_A.plug_state * ldep_A

                plugsite_B = ch.cen_B.plugsites[m]
                ldep_B = plugsite_B.calc_ldep()
                p_nmB = plugsite_B.plug_state * ldep_B

                #spbs diag terms:
                A[0,0] += - p_nmA - p_nmB
                #A side
                A[self._idx[(0,n,m)],
                  self._idx[(0,n,m)]] = - muk + p_nmA#- alphaR - betaR#sign?
                A[0, self._idx[(0,n,m)]] = p_nmA#alphaR - betaR
                A[self._idx[(0,n,m)], 0] = A[0, self._idx[(0,n,m)]]
                A[self._idx[(0,n,m)],
                  self._idx[(0,n)]] = muk
                A[self._idx[(0,n)],
                  self._idx[(0,n,m)]] = muk
                #B side
                A[self._idx[(1,n,m)],
                  self._idx[(1,n,m)]] = - muk + p_nmB# - alphaL - betaL
                A[0, self._idx[(1,n,m)]] = p_nmB#- alphaL + betaL
                A[self._idx[(1,n,m)], 0] = A[0, self._idx[(1,n,m)]]
                A[self._idx[(1,n,m)],
                  self._idx[(1,n)]] = muk
                A[self._idx[(1,n)],
                  self._idx[(1,n,m)]] = muk
        return A#.tocsr()

    def calc_C(self):
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef float Fmz = self.params['Fmz']
        cdef float d0 = self.params['d0']
        cdef float kappa_c = self.params['kappa_c']
        cdef float ldep_A, ldep_B
        C = np.zeros(1 + N * (1 + Mk) * 2)
        C[0] = 2 * Fmz
        for n in range(N):
            ch = self.chromosomes[n]
            C[self._idx[(0,n)]] = kappa_c * d0
            C[self._idx[(1,n)]] = - kappa_c * d0
            for m in range(Mk):
                plugsite_A = ch.cen_A.plugsites[m]
                ldep_A = plugsite_A.calc_ldep()
                p_nmA = plugsite_A.plug_state * ldep_A

                plugsite_B = ch.cen_B.plugsites[m] 
                ldep_B = plugsite_B.calc_ldep()
                p_nmB = plugsite_B.plug_state * ldep_B

                C[0] += - p_nmA - p_nmB
                C[self._idx[(0,n,m)]] = p_nmA
                C[self._idx[(1,n,m)]] = p_nmB
        return C

    def calcb(self):
        B = self.B_mat
        C = self.calc_C()
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

    def plug_unplug(self, int time_point):
        """Let's play dices ...
        """
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])

        for n in range(N):
            ch = self.chromosomes[n]
            for m in range(Mk):
                ch.cen_A.plugsites[m].plug_unplug(time_point)
                ch.cen_B.plugsites[m].plug_unplug(time_point)

    def position_update(self, np.ndarray speeds, int time_point):
        """given the speeds obtained by solving Atot.x = btot
        and caclulated switch events 
        """
        cdef int Mk = int(self.params['Mk'])
        cdef int N = int(self.params['N'])
        cdef double dt = self.params['dt']
        cdef double Vk = self.params['Vk']

        speeds *= Vk * dt #Back to real space
        self.spbR.set_pos(self.spbR.pos + speeds[0], time_point)
        self.spbL.set_pos(self.spbL.pos - speeds[0], time_point)
        for n in range(N):
            ch = self.chromosomes[n]
            ch.cen_A.set_pos(ch.cen_A.pos + speeds[self._idx[(0,n)]],
                             time_point)
            ch.cen_B.set_pos(ch.cen_B.pos + speeds[self._idx[(1,n)]],
                             time_point)
            for m in range(Mk):
                plugsite = ch.cen_A.plugsites[m]
                new_pos = plugsite.pos + speeds[self._idx[(0,n,m)]]
                plugsite.set_pos(new_pos, time_point)
                plugsite = ch.cen_B.plugsites[m]
                new_pos = plugsite.pos + speeds[self._idx[(1,n,m)]]
                plugsite.set_pos(new_pos, time_point)
            ch.pluged_history[time_point] = ch.pluged()
            ch.mero_history[time_point] = ch.mero()
        self.step_count += 1

    def get_sim_duration(self):
        return (len(self.spbR.traj) - 1)*self.params['dt']


cdef class Organite(object):
    """Base class for all the physical elements of the spindle

    Parameters
    ----------
    parent : an other subclass of :class:`Organite`
        from which the parameters are inheritated.
    init_pos : float, initial position

    Attributes
    ----------
    KD : a :class:`~spindle_dynamics.KinetoDynamics` instance
    pos : float, the position
    traj : ndarrat, the trajectory
    
    Methods
    -------
    set_pos(pos, time_point) : sets the position and updates the trajectory
    get_pos(time_point): returns the position at `time_point`
    """

    cdef int num_steps
    cdef double init_pos
    
    def __init__(self, parent, double init_pos=0.):
        """
        
        """
        self.parent = parent
        self.KD = parent.KD
        self.num_steps = parent.KD.num_steps
        self.pos = init_pos
        self.traj = np.zeros(self.num_steps)
        self.traj[0] = init_pos
        
    cpdef set_pos(self, double pos, int time_point=0):
        """sets the position. If `time_point` is provided, sets
        the corresponding value in self.traj[time_point]
        """
        self.pos = pos
        if self.pos < self.KD.spbL.pos: 
            self.pos = self.KD.spbL.pos
        elif pos > self.KD.spbR.pos:
            self.pos = self.KD.spbR.pos
        if time_point > 0:
            self.traj[time_point] = pos

    cpdef double get_pos(self, int time_point=0):
        """Returns the position.
        
        If `time_point` is null (default), returns the current position 
        If `time_point` is not null, returns the position at time_point
        """
        if time_point == 0:
            return self.pos
        return self.traj[time_point]


        
class Centromere(Organite):
    """ 
    The centromere is where the plugsites are bound to the
    chromosome and where the cohesin spring restoring force, as
    well as the friction coefficient, are applied.
    This is a subclass of :class:`Organite`
    
    Parameters:
    ----------
    chromosome: a :class:`~Chromosome` instance
       the parent chromosome
    tag: {'A', 'B'}
       Side of the centromere. Note that the centromere
       side and the SPB side are not necesseraly related
    """
    def __init__(self, chromosome, tag='A', plug=None):

        self.tag = tag
        self.chromosome = chromosome
        d0 = self.chromosome.KD.params['d0']
        if tag == 'A':
            init_pos = chromosome.pos + d0 / 2.
        elif tag == 'B':
            init_pos = chromosome.pos - d0 / 2.
        else:
            raise ValueError("the `tag` attribute must be 'A' or 'B'.")
        Organite.__init__(self, chromosome, init_pos)
        Mk = self.KD.params['Mk']
        self.toa = 0 #time of arrival at pole
        self.plugsites = [PlugSite(self, initial_plug = self.KD.initial_plug)
                          for m in range(Mk)]
        self.traj = np.zeros(self.KD.num_steps)

    def is_attached(self):
        """
        Returns True if at least one plugsite is attached
        to at least one SPB
        """
        for plugsite in self.plugsites:
            if plugsite.plug_state != 0:
                return True
        return False

    def P_attachleft(self):
        orientation = self.KD.params['orientation']
        if orientation == 0: return 0.5
        cdef double lp, rp
        lp = self.left_pluged()
        rp = self.right_pluged()
        if lp + rp == 0:
            return 0.5
        cdef double P_left
        P_left = 0.5 + orientation * (lp - rp) / (2 * (lp + rp)) 
        return P_left

    def left_pluged(self):
        cdef int lp
        lp = -sum([plugsite.plug_state for plugsite
                   in self.plugsites if plugsite.plug_state ==  -1])
        return float(lp)

    def right_pluged(self):
        cdef int lp
        rp = sum([plugsite.plug_state for plugsite
                  in self.plugsites if plugsite.plug_state == 1])
        return float(rp)

    def at_rightpole(self, double tol=0.01):
        rightpole_pos = self.KD.spbR.pos
        if abs(self.pos - rightpole_pos) < tol:
            return True
        return False
    
    def at_leftpole(self, double tol=0.01):
        leftpole_pos = self.KD.spbR.pos
        if abs(self.pos - leftpole_pos) < tol:
            return True
        return False
    

class Chromosome(Organite):
    """The chromosome, containing two centromeres ('A' and 'B')

    Parameters
    ----------
    spindle: a :class:`~Spindle` instance
    
    """

    def __init__(self, spindle, plug = None):
        
        d0 = spindle.KD.params['d0']
        L0 = spindle.KD.params['L0']
        Mk = spindle.KD.params['Mk']
        center_pos = random.gauss(0, 0.2 * (L0 - d0))
        Organite.__init__(self, spindle, center_pos)
        self.cen_A = Centromere(self, 'A')
        self.cen_B = Centromere(self, 'B')
        
        
        self.pluged_history = np.zeros((self.KD.num_steps, 2))
        self.pluged_history[0] = self.pluged()
        self.mero_history = np.zeros((self.KD.num_steps, 2))
        self.mero_history[0] = self.mero()

    def is_right_A(self):
        """returns True if centromere A is pluged
        mainly to the right pole.
        """
        right_A = self.cen_A.right_pluged() + self.cen_B.left_pluged()
        left_A = self.cen_B.right_pluged() + self.cen_A.left_pluged()
        if right_A >= left_A:
            return True
        return False
        
    def pluged(self):
        """returns the number of correctly pluged MTs 
        """
        if self.is_right_A():
            return self.cen_A.right_pluged(), self.cen_B.left_pluged() 
        else:
            return self.cen_A.left_pluged(), self.cen_B.right_pluged() 

    def mero(self):
        """returns the number of erroneously pluged MTs 
        """
        if self.is_right_A():
            return self.cen_A.left_pluged(), self.cen_B.right_pluged()
        else:
            return self.cen_A.right_pluged(), self.cen_B.left_pluged()
        
    def pair_dist(self):
        return abs(self.cen_A.pos - self.cen_B.pos)

    def plug_dist(self, double plugpos):
        return abs(self.center() - plugpos)

    def center(self):
        return (self.cen_A.pos + self.cen_B.pos)/2

    def cen_traj(self):
        return (self.cen_A.traj + self.cen_B.traj)/2

    def at_rightpole(self, double tol) :
        """tol : tolerance distance
        """
        if self.cen_A.at_rightpole(tol) or self.cen_B.at_rightpole(tol):
            return True
        return False
    
    def at_leftpole (self, double tol) :
        if self.cen_A.at_leftpole(tol) or self.cen_B.at_leftpole(tol):
            return True
        return False
        
    def at_pole(self, side=None, double tol=0.0):
        if side == 1 and self.at_rightpole(tol):
            return True
        elif side == -1 and self.at_leftpole(tol):
            return True
        elif side is None:
            if self.at_rightpole(tol) and self.at_leftpole(tol):
                return True
        else :
            return False

class PlugSite(Organite):

    """An attachment site object.

    Parameters:
    -----------
    cen: a Centromere instance
    plug : string
        The `initial_plug` argument can be either of the following:
        None (the default), in which case the PlugSite is randomly attached
        * 'null': self.plug = 0
        * 'amphitelic' -- self.plug_state = -1 for A side plugsites and
            self.plug_state = 1 for B side ones 
        * 'random' -- self.plug = -1, 0, or 1 with equal probability
        * 'monotelic' -- self.plug = -1 for A side plugsites
              and self.plug = 0 for  B side ones
        * 'syntelic' --  self.plug = 1 in any case
    
    """

    def __init__(self, centromere, initial_plug):
        Organite.__init__(self, centromere, centromere.pos)
        self.centromere = centromere
        self.tag = self.centromere.tag
        if initial_plug == None: 
            self.plug_state = random.randint(-1, 1)
        elif initial_plug == 'null':
            self.plug_state = 0
        elif initial_plug == 'amphitelic':
            self.plug_state = 1 if self.tag == 'A' else -1
        elif initial_plug == 'random':
            self.plug_state = random.randint(-1, 1)
        elif initial_plug == 'monotelic':
            self.plug_state = 1 if self.tag == 'A' else 0
        elif initial_plug == 'syntelic':
            self.plug_state = 1
        elif initial_plug == 'merotelic':
            self.plug_state = random.choice([-1,1])
        else:
            self.plug_state = initial_plug
        
        self.state_hist = np.zeros(self.KD.num_steps)
        self.state_hist[:] = self.plug_state
        self.P_att = 1 - np.exp(- self.KD.params['k_a'])

    def set_plug_state(self, state, time_point):
        self.plug_state = state
        self.state_hist[time_point:] = state

    def calc_ldep(self):
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
        ld_slope = self.KD.params['ld_slope']
        if ld_slope == 0: return 1.
        ld0 = self.KD.params['ld0']
        cdef double mt_length
        cdef double pole_pos
        pole_pos = self.KD.spbR.pos * self.plug_state
        mt_length = abs(pole_pos - self.pos)
        cdef double ldep
        ldep = ld_slope * mt_length + ld0
        # TODO: investigate: introduces a first order discontinuity
        #     suspected to trigger artifacts when the plugsite
        #     is close to the pole
        if mt_length < 0.0001:
            ldep = mt_length  # No force when at pole
        return ldep

    def plug_unplug(self, time_point):
        dice = random.random()
        #Attachment
        if self.plug_state == 0 and dice < self.P_att:
            side_dice = random.random()
            P_left = self.centromere.P_attachleft()
            if dice < P_left:
                self.set_plug_state(-1, time_point)
            else:
                self.set_plug_state(1, time_point)
        #Detachment
        elif self.plug_state != 0 and dice < self.P_det():
            self.set_plug_state(0, time_point)
                    
    def P_det(self):
        cdef d_alpha = self.KD.params['d_alpha']
        cdef double k_d0 = self.KD.params['k_d0']
        if d_alpha == 0: return k_d0
        
        cdef double k_dc
        cdef dist = abs(self.pos - self.centromere.pos)
        if dist == 0: return 1.
        k_dc = k_d0  *  d_alpha / dist
        if k_dc > 1e4: return 1.
        return 1 - np.exp(-k_dc) 


class Spb(Organite) :

    """ A spindle pole object.

    Attributes:
    side: 1 or -1 (right or left. In the article, -1 correspond to the
    daggered variable)
    """
    
    def __init__(self, spindle, side, L0):
        """ side = 1 : left
        side = -1 : right
        L0 : spindle length
        """
        self.side = side
        init_pos = side * L0 / 2.
        Organite.__init__(self, spindle, init_pos)
        
class Spindle(object) :

    def __init__(self, KD):
        self.KD = KD
        
    def all_plugsites(self):
        plugsites = []
        for ch in self.KD.chromosomes:
            for plugsite in ch.cen_A.plugsites:
                plugsites.append(plugsite) 
            for plugsite in ch.cen_B.plugsites:
                plugsites.append(plugsite) 
        return plugsites
