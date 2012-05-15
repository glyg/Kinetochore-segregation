class Organite(Object):
    """Base class for all the physical elements of the spindle

    """
    def __init__(self, parent, init_pos=0):
        """

        """
        self.parent = parent
        self.num_steps = parent.num_steps
        self.KD = parent.KD
        self.traj = np.zeros(self.num_steps)
        self.traj[0] = init_pos
        
    def set_pos(double pos, int time_point):
        if pos < self.KD.spbL.pos:
            pos = self.KD.spbL.pos
        elif pos < self.KD.spbR.pos:
            pos = self.KD.spbR.pos
        self.pos = pos
        self.traj[time_point] = pos

    def get_pos(time_point):
        return self.traj[time_point]

        
class Centromere(Organite):
    """The centromere is where each plugsite is attached to the
    chromosome and where the cohesin spring restoring force, as
    well as the friction coefficient are applied
    
    """
    def __init__(self, d0, chromosome, Mk, tag='A', plug=None):

        self.tag = tag
        if tag == 'A':
            init_pos = self.chromosome.pos - d0 / 2.
        elif tag == 'B':
            init_pos = self.chromosome.pos + d0 / 2.
        else:
            raise ValueError('the `tag` attribute must be 'A' or 'B' ')
        Organite.__init__(chromosome, init_pos)

        self.toa = 0 #time of arrival at pole
        self.plugsites = [PlugSite(cen=self, plug=plug)
                          for m in range(Mk)]
        self.traj = np.zeros(self.num_steps)

    def is_attached(self):
        for plugsite in self.plugsites:
            if plugsite.plugstate != 0:
                return True
        return False

    def P_attachleft(self):
        orientation = self.params['orientation']
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
        return lp

    def right_pluged(self):
        cdef int lp
        rp = sum([plugsite.plug_state for plugsite
                  in self.plugsites if plugsite.plug_state == 1])
        return rp

    def at_rightpole(self, double tol=0.001):
        rightpole_pos = self.chromosome.KD.spbR.pos
        is abs(self.pos - rightpole_pos) < tol:
            return True
        return False
    
    def at_leftpole(self, double tol=0.001):
        leftpole_pos = self.chromosome.KD.spbR.pos
        is abs(self.pos - rightpole_pos) < tol:
            return True
        return False
    

class Chromosome(Organite):
    """Chromosome Object """

    def __init__(self, KD, plug = None):
        d0 = KD.params['d0']
        L0 = KD.params['L0']
        Mk = KD.params['Mk']
        center_pos = random.gauss(0, 0.2 * (L0 - d0))

        Organite.__init__(KD, center_pos)
        self.cen_A = Centromere(d0, self, Mk, 'A', plug)
        self.cen_B = Centromere(d0, self, Mk, 'B', plug)
        
        
        self.pluged_history = np.zeros(self.num_steps)
        self.pluged_history[0] = self.pluged()
        self.mero_history = np.zeros(self.num_steps)
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
            return self.cen_A.left_pluged, self.cen_B.right_pluged 
        else:
            return self.cen_A.right_pluged, self.cen_B.left_pluged 
        
    def pair_dist(self):
        return abs(self.cen_A.pos - self.cen_B.pos)

    def plug_dist(self, double plugpos):
        return abs(self.center() - plugpos)

    def center(self):
        return (self.cen_A.pos + self.cen_B.pos)/2

    def cen_traj(self):
        return (self.cen_A.traj + self.cen_B.traj)/2

    def at_rightpole(self, double tol=0.0) :
        """tol : tolerance distance
        """
        if self.cen_A.at_rightpole(tol) or self.cen_B.at_rightpole(tol):
            return True
        return False
    
    def at_leftpole (self, double tol=0.0) :
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

    Instanciation: PlugSite(ch, side, doi, plug).
    The plug argument can be either of the following:
    None (the default), in which case the PlugSite is randomly attached
    * 'null': self.plug = 0
    * 'amphitelic' -- self.plug_state = -1 for A side PlugSite and
        self.plug_state = 1 for Bside PlugSite
    * 'random' -- self.plug = -1, 0, or 1 with equal probability
    * 'monotelic' -- self.plug = -1 for A side PlugSite
          and self.plug = 0 for  B side ones
    * 'syntelic' --  self.plug = 1 
    
    """

    def __init__(self, cen, plug = None):
        Organite.__init__(cen, cen.pos)

        if plug == None: 
            self.plug_state = random.randint(-1, 1)
        elif plug == 'null':
            self.plug_state = 0
        elif plug == 'amphitelic':
            self.plug_state = -1 if self.tag == 'A' else 1
        elif plug == 'random':
            self.plug_state = random.randint(-1, 1)
        elif plug == 'monotelic':
            self.plug_state = -1 if self.tag == 'A' else 0
        elif plug == 'syntelic':
            self.plug_state = 1
        elif plug == 'merotelic':
            self.plug_state = random.choice([-1,1])
        else:
            self.plug_state = plug
        
        self.state_hist = np.zeros(self.num_steps)
        self.state_hist[:] = self.plug_state
        self.P_att = 1 - np.exp(self.params['k_a'])

    def set_plug_state(state, time_point):
        self.plug_state = state
        self.state_hist[time_point:] = state
        
    ### This should be a method of PlugSite
    #@cython.profile(False)
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
        cdef double dist
        cdef double pole_pos
        pole_pos = self.KD.spbR.pos * self.plug_state
        dist = abs(pole_pos - self.pos)
        cdef double ldep
        ldep = ld_slope * dist + ld0
        # TODO: investigate: introduces a first order discontinuity
        #     suspected to trigger artifacts when the plugsite
        #     is close to the pole
        if dist < 0.0001:
            ldep = dist  # No force when at pole
        return ldep

    def plug_unplug(self):
        dice = random.random()
        #Attachment
        if self.plug_state == 0 and dice < self.P_att:
            side_dice = random.random()
            P_left = self.cen.P_attachleft()
            if dice < P_left:
                self.set_plug_state(-1, time_point)
            else:
                self.set_plug_state(1, time_point)
        #Detachment
        elif self.plug_state != 0 and dice < self.P_det():
            self.set_plug_state(0, time_point)
                    
    def P_det(self):
        
        cdef d_alpha = self.params['d_alpha']
        cdef double k_d0 = self.params['k_d0']
        if d_alpha == 0: return k_d0
        
        cdef double k_dc
        cdef dist = abs(self.pos - self.cen.pos)
        k_dc = k_d0  *  d_alpha / dist
        if k_dc > 1e4: return 1.
        
        return 1 - np.exp(-k_dc) 


class Spb(object) :

    """ A spindle pole object.

    Attributes:
    side: 1 or -1 (right or left. In the article, -1 correspond to the
    daggered variable)
    """
    
    def __init__(self, KD, side, L0):
        """ side = 1 : left
        side = -1 : right
        L0 : spindle length
        """
        self.side = side
        init_pos = side * L0 / 2
        Organite.__init__(KD, init_pos)
        
class Spindle(object) :

    def __init__(self, spbR, spbL ):
        self.spbR = spbR
        self.spbL = spbL
        
    def length(self):
        return self.spbR.pos - self.spbL.pos

    def length_traj(self):
        return self.spbR.traj - self.spbL.traj

    
