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

RIGHT = 1
LEFT = -1
a = 0
b = 1


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
        duration = self.params['span']
        dt = self.params['dt']
        self.num_steps = int(duration/dt)
        self.spindle = Spindle(self)
        self.spbR = Spb(self.spindle, RIGHT, L0) #right spb (RIGHT = 1)
        self.spbL = Spb(self.spindle, LEFT, L0) #left one (LEFT = -1)
        self.initial_plug = initial_plug
        self.chromosomes = [Chromosome(self.spindle)
                            for n in range(N)]
        self._idx = self._calc_idx()
        cdef int dim = 1 + N * ( 1 + Mk ) * 2
        self.calc_B()
        self.A0_mat = self.time_invariantA()
        self.At_mat = np.zeros((dim, dim), dtype = float)
        self.step_count = 0
        self.anaphase = False
        self.spindle.all_plugsites = self.spindle.get_all_plugsites()

        
    def _calc_idx(self):
        """Returns the  index dictionnary
        """
        cdef int Mk, N
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        idxs = {}
        for n in range(N):
            idxs[(a,n)] =  2*n*(Mk+1) + 1
            idxs[(b,n)] = 2*n*(Mk+1) + 2 + Mk
            for m in range(Mk):
                idxs[(a,n,m)] = 2*n*(Mk + 1) + 2 + m
                idxs[(b,n,m)] =  2*n*(Mk + 1) + 3 + Mk + m
        return idxs

    def one_step(self, time_point):
        """elementary step"""
        if not self.anaphase:
            self.plug_unplug(time_point)
        self.solve()
        self.position_update(time_point)

    def solve(self):
        X = self.get_state_vector()
        A = self.calc_A()
        B = self.B_mat
        C = self.calc_C()
        pos_dep = np.dot(B,X) + C
        self.speeds = linalg.solve(A, -pos_dep)

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

    def calc_A(self):
        self.time_dependentA()
        A = self.A0_mat + self.At_mat
        return A

    def time_invariantA(self):
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef double muc = self.params['muc']
        cdef float muk = self.params['muk']
        cdef float mus = self.params['mus']
        cdef float Vmz = self.params['Vmz']
        cdef float Fmz = self.params['Fmz']
        cdef int dims
        dims = 1 + 2 * N * ( Mk + 1 ) 
        A0 = np.zeros((dims, dims))
        A0[0,0] = - 2 * mus - 4 * Fmz / Vmz
        for n in range(N):
            ch = self.chromosomes[n]
            A0[self._idx[(0,n)],
              self._idx[(0,n)]] = - Mk * muk - muc
            A0[self._idx[(1,n)],
              self._idx[(1,n)]] = - Mk * muk - muc
            for m in range(Mk):
                #A side
                A0[self._idx[(0,n,m)],
                   self._idx[(0,n,m)]] = - muk
                A0[self._idx[(0,n,m)],
                   self._idx[(0,n)]] = muk
                A0[self._idx[(0,n)],
                   self._idx[(0,n,m)]] = muk
                #B side
                A0[self._idx[(1,n,m)],
                   self._idx[(1,n,m)]] = - muk
                A0[self._idx[(1,n,m)],
                   self._idx[(1,n)]] = muk
                A0[self._idx[(1,n)],
                   self._idx[(1,n,m)]] = muk
        return A0

        
    def time_dependentA(self):

        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        self.At_mat[0,0] = 0
        cdef int pi_nmA, pi_nmB
        cdef int pluggedA, pluggedB
        #cdef float ldep_A, ldep_B
        for n in range(N):
            ch = self.chromosomes[n]
            for m in range(Mk):
                plugsite_A = ch.cen_A.plugsites[m]
                #ldep_A = plugsite_A.calc_ldep()
                pi_nmA = plugsite_A.plug_state# * ldep_A
                pluggedA = plugsite_A.plugged
                plugsite_B = ch.cen_B.plugsites[m]
                #ldep_B = plugsite_B.calc_ldep()
                pi_nmB = plugsite_B.plug_state# * ldep_B
                pluggedB = plugsite_B.plugged
                #spbs diag terms:
                self.At_mat[0,0] -= pluggedA + pluggedB
                #A side
                self.At_mat[self._idx[(0,n,m)],
                            self._idx[(0,n,m)]] = - pluggedA
                self.At_mat[0, self._idx[(0,n,m)]] = pi_nmA
                self.At_mat[self._idx[(0,n,m)], 0] = pi_nmA
                #B side
                self.At_mat[self._idx[(1,n,m)],
                  self._idx[(1,n,m)]] = -  pluggedB
                self.At_mat[0, self._idx[(1,n,m)]] = pi_nmB
                self.At_mat[self._idx[(1,n,m)], 0] = pi_nmB

    def calc_B(self):
        """Returns the matrix containing the linear terms
        of the equation set A\dot{X} + BX + C = 0
        """
        #self.time_dependentB()
        kappa_k = self.params['kappa_k']
        kappa_c = self.params['kappa_c']
        if kappa_c > 0 :
            self.B_mat = kappa_k * self.kinetochore_B() + kappa_c * self.cohesin_B()
        else:
            self.B_mat = kappa_k * self.kinetochore_B()

    def kinetochore_B(self):
        """Returns the constant part of the matrix containing the
        linear terms of the equation set A\dot{X} + BX + C = 0
        """
        cdef N, Mk
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        dim = 1 + N * ( 1 + Mk ) * 2
        Bk = np.zeros((dim, dim), dtype = float)
        for n in range(N) :
            ch = self.chromosomes[n]
            Bk[self._idx[(0, n)],
              self._idx[(0, n)]] = - Mk
            Bk[self._idx[(1, n)],
              self._idx[(1, n)]] = - Mk
            for m in range(Mk):
                Bk[self._idx[(0, n, m)],
                  self._idx[(0, n, m)]] = - 1
                Bk[self._idx[(1, n, m)],
                  self._idx[(1, n, m)]] = - 1
                Bk[self._idx[(0, n)],
                  self._idx[(0, n, m)]] = 1
                Bk[self._idx[(0, n, m)],
                  self._idx[(0, n)]] = 1
                Bk[self._idx[(1, n)],
                  self._idx[(1, n, m)]] = 1
                Bk[self._idx[(1, n, m)],
                  self._idx[(1, n)]] = 1
        return Bk

    def cohesin_B(self):
        cdef N, Mk
        Mk = int(self.params['Mk'])
        N = int(self.params['N'])
        dim = 1 + N * ( 1 + Mk ) * 2
        Bc = np.zeros((dim, dim), dtype = float)
        for n in range(N) :
            Bc[self._idx[(0, n)],
                   self._idx[(0, n)]] = - 1
            Bc[self._idx[(1, n)],
                   self._idx[(1, n)]] = - 1
            Bc[self._idx[(1, n)],
                   self._idx[(0, n)]] = 1
            Bc[self._idx[(0, n)],
                   self._idx[(1, n)]] = 1
        return Bc 
            
    def calc_C(self):
        """Returns the vector containing the constant terms
        of the equation set A\dot{X} + BX + C = 0
        """
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        cdef float Fmz = self.params['Fmz']
        cdef float d0 = self.params['d0']
        cdef float kappa_c = self.params['kappa_c']
        #cdef float ldep_A, ldep_B
        C = np.zeros(1 + N * (1 + Mk) * 2)
        C[0] = 2 * Fmz
        for n in range(N):
            ch = self.chromosomes[n]
            delta = ch.delta()
            
            C[self._idx[(0,n)]] = - delta * kappa_c * d0
            C[self._idx[(1,n)]] = delta * kappa_c * d0
            for m in range(Mk):
                plugsite_A = ch.cen_A.plugsites[m]
                #ldep_A = plugsite_A.calc_ldep()
                pi_nmA = plugsite_A.plug_state # * ldep_A
                pluggedA = plugsite_A.plugged
                plugsite_B = ch.cen_B.plugsites[m] 
                #ldep_B = plugsite_B.calc_ldep()
                pi_nmB = plugsite_B.plug_state # * ldep_B
                pluggedB = plugsite_B.plugged
                C[0] -= pluggedA + pluggedB
                C[self._idx[(0,n,m)]] = pi_nmA
                C[self._idx[(1,n,m)]] = pi_nmB
        return C

    def plug_unplug(self, int time_point):
        """Let's play dices ...
        """
        cdef int N = int(self.params['N'])
        cdef int Mk = int(self.params['Mk'])
        for plugsite in self.spindle.all_plugsites:
            plugsite.plug_unplug(time_point)

    def position_update(self, int time_point):
        """given the speeds obtained by solving Atot.x = btot
        and caclulated switch events 
        """
        cdef int Mk = int(self.params['Mk'])
        cdef int N = int(self.params['N'])
        cdef double dt = self.params['dt']
        cdef double Vk = self.params['Vk']

        speeds = self.speeds
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
            ch.plugged_history[time_point] = ch.plugged()
            ch.mero_history[time_point] = ch.mero()

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
    
    def __init__(self, parent, double init_pos):
        """
        
        """
        self.parent = parent
        self.KD = parent.KD
        self.num_steps = parent.KD.num_steps
        self.traj = np.zeros(self.num_steps)
        self.pos = init_pos
        self.traj[0] = init_pos
        
    cpdef set_pos(self, double pos, int time_point=-1):
        """sets the position. If `time_point` is provided, sets
        the corresponding value in self.traj[time_point]
        """
        self.pos = pos
        if pos > self.KD.spbR.pos:
            self.pos = self.KD.spbR.pos
        elif self.pos < self.KD.spbL.pos: 
            self.pos = self.KD.spbL.pos
        if time_point >= 0:
            self.traj[time_point] = self.pos

    cpdef double get_pos(self, int time_point=0):
        """Returns the position.
        
        If `time_point` is null (default), returns the current position 
        If `time_point` is not null, returns the position at time_point
        """
        if time_point == 0:
            return self.pos
        return self.traj[time_point]

class Spindle(object) :

    def __init__(self, KD):
        self.KD = KD
        
    def get_all_plugsites(self):
        plugsites = []
        for ch in self.KD.chromosomes:
            for plugsite in ch.cen_A.plugsites:
                plugsites.append(plugsite) 
            for plugsite in ch.cen_B.plugsites:
                plugsites.append(plugsite) 
        return plugsites

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

        self.plugged_history = np.zeros((self.KD.num_steps, 2))
        self.plugged_history[0] = self.plugged()
        self.mero_history = np.zeros((self.KD.num_steps, 2))
        self.mero_history[0] = self.mero()

    def is_right_A(self):
        """returns True if centromere A is plugged
        mainly to the right pole.
        """
        right_A = self.cen_A.right_plugged() + self.cen_B.left_plugged()
        left_A =  self.cen_A.left_plugged() + self.cen_B.right_plugged()
        if right_A >= left_A:
            return True
        return False

    def delta(self):
        """In case the centromeres swap (exchange side), the direction
        of the cohesin restoring force needs to be changed
        """
        return 1 if self.cen_A.pos < self.cen_B.pos else -1
    

        
    def plugged(self):
        """returns the number of correctly plugged MTs 
        """
        if self.is_right_A():
            return self.cen_A.right_plugged(), self.cen_B.left_plugged() 
        else:
            return self.cen_A.left_plugged(), self.cen_B.right_plugged() 

    def mero(self):
        """returns the number of erroneously plugged MTs 
        """
        if self.is_right_A():
            return self.cen_A.left_plugged(), self.cen_B.right_plugged()
        else:
            return self.cen_A.right_plugged(), self.cen_B.left_plugged()
        
    def pair_dist(self):
        return abs(self.cen_A.pos - self.cen_B.pos)

    def signed_dist(self):
        return self.cen_A.pos - self.cen_B.pos

    def plug_dist(self, double plugpos):
        return abs(self.center() - plugpos)

    def center(self):
        return (self.cen_A.pos + self.cen_B.pos)/2

    def center_traj(self):
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
            init_pos = chromosome.pos - d0 / 2.
        elif tag == 'B':
            init_pos = chromosome.pos + d0 / 2.
        else:
            raise ValueError("the `tag` attribute must be 'A' or 'B'.")
        Organite.__init__(self, chromosome, init_pos)
        Mk = self.KD.params['Mk']
        self.toa = 0 #time of arrival at pole
        self.state_vector = np.zeros(Mk)
        self.plugsites = [PlugSite(self) for m in range(Mk)]
        self.calc_state_vector()
        
    def is_attached(self):
        """
        Returns True if at least one plugsite is attached
        to at least one SPB
        """
        for plugsite in self.plugsites:
            if plugsite.plug_state != 0:
                return True
        return False

    def calc_state_vector(self):
        cdef np.ndarray state
        state = np.array([plugsite.plug_state for plugsite
                          in self.plugsites])
        self.state_vector = state
        
    def P_attachleft(self):

        orientation = self.KD.params['orientation']
        if orientation == 0: return 0.5
        cdef float lp, rp
        self.calc_state_vector()
        lp = self.left_plugged()
        rp = self.right_plugged()
        if lp + rp == 0:
            return 0.5
        cdef float P_left
        P_left = 0.5 + orientation * (lp - rp) / (2 * (lp + rp)) 
        return P_left

    def left_plugged(self):
        cdef float lp 
        cdef np.ndarray left_plugged
        left_plugged = self.state_vector * (self.state_vector - 1) / 2.
        lp = left_plugged.sum()
        return lp
        
    def right_plugged(self):
        cdef float rp
        right_plugged = self.state_vector * (1 + self.state_vector) / 2.
        rp = right_plugged.sum()
        return rp

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
        

class PlugSite(Organite):

    """An attachment site object.

    Parameters:
    -----------
    centromere: a Centromere instance
    """
    def __init__(self, centromere):
        cdef double init_pos = centromere.pos
        Organite.__init__(self, centromere, init_pos)
        initial_plug = self.KD.initial_plug
        self.centromere = centromere
        self.tag = self.centromere.tag
        if initial_plug == None: 
            self.plug_state = random.randint(-1, 1)
        elif initial_plug == 'null':
            self.plug_state = 0
        elif initial_plug == 'amphitelic':
            self.plug_state = - 1 if self.tag == 'A' else 1
        elif initial_plug == 'random':
            self.plug_state = random.randint(-1, 1)
        elif initial_plug == 'monotelic':
            self.plug_state = - 1 if self.tag == 'A' else 0
        elif initial_plug == 'syntelic':
            self.plug_state = 1
        elif initial_plug == 'merotelic':
            self.plug_state = random.choice([-1,1])
        else:
            self.plug_state = initial_plug
        self.set_pos(init_pos)
        self.state_hist = np.zeros(self.KD.num_steps)
        self.state_hist[:] = self.plug_state
        self.P_att = 1 - np.exp(- self.KD.params['k_a'])

    def set_plug_state(self, state, int time_point = -1):
        self.plug_state = state
        self.plugged = 0 if state == 0 else 1
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
            if side_dice < P_left:
                self.set_plug_state(-1, time_point)
            else:
                self.set_plug_state(1, time_point)
        #Detachment
        elif dice < self.P_det():
            self.set_plug_state(0, time_point)
                    
    def P_det(self):
        cdef d_alpha = self.KD.params['d_alpha']
        cdef double k_d0 = self.KD.params['k_a']
        if d_alpha == 0: return k_d0
        
        cdef double k_dc
        cdef dist = np.abs(self.pos - self.centromere.pos)
        if dist == 0: return 1.
        k_dc = k_d0  *  d_alpha / dist
        if k_dc > 1e4: return 1.
        return 1 - np.exp(-k_dc) 
