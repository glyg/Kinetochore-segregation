#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module provides the core simulation functionalities.

See Gay et al. J. Cell Biol., 2012 http://dx.doi.org/10.1083/jcb.201107124
The original framework was adapted from:
Civelekoglu-Scholey et al. Biophys. J. 90(11), 2006
http://dx.doi.org/10.1529/biophysj.105.078691

"""

import time
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal
from xml.etree.ElementTree import Element, SubElement, tostring
from Image import fromarray as image_fromarray

import pyximport
pyximport.install()
## local imports
from .spindle_dynamics import KinetoDynamics
from .xml_handler import ParamTree, indent, ResultTree
from ..analysis.eval_simul import evaluations # as eval_simul


__all__ = ["Metaphase", "reduce_params", "PARAMFILE",
           "MEASUREFILE", "get_fromfile"]

CURRENT_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.dirname(CURRENT_DIR)
PARAMFILE = os.path.join(ROOT_DIR, 'default', 'params.xml')
MEASUREFILE = os.path.join(ROOT_DIR, 'default', 'measures.xml')
MEASURETREE = ParamTree(MEASUREFILE, adimentionalized = False)
MEASURES = MEASURETREE.absolute_dic


class Metaphase(object):
    """An instance of the Metaphase class is a wrapper around
    the whole simulation.

    Typical usage :
    ---------------

    >>> from kt_simul.simul_spindle import Metaphase
    >>> m = Metaphase()
    >>> m.simul()
    >>> m.show_trajs()
    >>> m.write_results('examples/docstring_results.xml',
                        'examples/docstring_data.npy')

    From an already runned simulation:
    ----------------------------------
    
    >>> from kt_simul.simul_spindle import Metaphase
    >>> m1 = get_fromfile('examples/docstring_results.xml')
    >>> m1.show_one(1) #This shows the trajactory of the chromosome 1
    >>> m2 = Metaphase(m1.paramtree, m1.measuretree) #A new simulation
    >>> m2.simul(ablat = 600) #this time with spindle ablation    

    Public methods:
    ---------------
    simul() : runs the simulation
    show_trajs() : displays the trajectories
    show_one() : displays one of the chromosomes' trajectories
    show_one_article() : another display style
    evaluate() : calculates the characteristics of the simulation
    write_results() : saves the trajectories, evaluations and parameters
    write_asoctave() : save the trajectories with a different data layout
    get_ch(n) : returns chromosome n object
    get_2Dtraj_list : returns 2D trajectories of the spindle elements
    
    Public attributes:
    ------------------
    report : a paragraph containing various output from the simulation
    KD : the spindle_dynamics.KinetoDynamics instance
    paramtree : a xml_handle.ParamTree instance
    measuretree : a xml_handle.MeasureTree instance
    num_steps : int
        Total number of time points
    timelapse : ndarray
        An array of num_steps time points with dt increment in seconds

    See Also:
    ---------
    spindle_dynamics.KinetoDynamics and xml_handler.ParamTree  
    
    """


    def __init__(self,  paramtree=None, measuretree=None,
                 paramfile=PARAMFILE, measurefile=MEASUREFILE,
                 initial_plug='random', reduce_p=True):

        """Metaphase instanciation method
        
        Key-word arguments:
        -------------------
        
        duration  : a float
            the duration of the mitosis in seconds (defaults to 900)

        paramtree : a ParamTree instance or None
            The paramtree contains the parameters for the simulation
            if paramtree is None, the parameters are read
            from the file paramfile. Defaults to None.
 
        measuretree : a ParamTree instance or None
            The measuretree contains the observed characteristics
            of the mitosis e.g. metaphase spindle elongation rate, etc.
            if measuretree is None, the measures are read from the file
            indicated by the measurefile argument. Defaults to None.

        paramfile : string
            Path to a xml file to read the parameters from. Defaults to the
            file params.xml in the module's default directory. Other parameter
            files can be produced by editing and changing the default one.
            If the paramtree argument is not None,  paramfile is ignored

        measurefile : string
            Path to a xml file to read the measures from. Defaults to the
            file measures.xml in the module's default directory.
            Other measure files can be produced by editing and changing
            the default one.
            If the measuretree argument is not None, measurefile is ignored

        initial_plug : string or None
            Defines globally the initial attachment states.
            This argument can have the following values: 
            - 'null': all kinetochores are detached
            - 'amphitelic': all chromosmes are amphitelic
            - 'random': all attachement site can be bound to
                        either pole or deteched with equal prob.
            - 'monotelic': right kinetochores are attached to the same pole,
                           left ones are detached
            - 'syntelic' : all kinetochores are attached to the same pole
            
        reduce_p  : bool
            If True, changes the parameters according to the measures
            so that the simulation average behaviour complies with
            the data in the measures dictionary

        """
        if paramtree is None:
            self.paramtree = ParamTree(paramfile)
        else:
            self.paramtree = paramtree
        if measuretree is None:
            self.measuretree = ParamTree(measurefile, adimentionalized=False)
        else:
            self.measuretree = measuretree
        if reduce_p:
            reduce_params(self.paramtree, self.measuretree)

        params = self.paramtree.relative_dic
        # Reset explicitely the unit parameters to their
        # dimentionalized value
        params['Vk'] = self.paramtree.absolute_dic['Vk']
        params['Fk'] = self.paramtree.absolute_dic['Fk']
        params['dt'] = self.paramtree.absolute_dic['dt']
            
        self.KD = KinetoDynamics(params, initial_plug=initial_plug)
        dt = self.paramtree.absolute_dic['dt']
        duration = self.paramtree.absolute_dic['span']
        self.num_steps = int(duration/dt)
        self.KD.anaphase = False
        self.timelapse = np.arange(0, duration, dt)
        self.report = []
        self.delay = -1
        self.observations = {}

    def __str__(self):
        lines = []
        lines.append('Metaphase class')
        try:
            lines.append('Parameters:')
            for line in str(self.paramtree.relative_dic).split(','):
                lines.append(line)
        except AttributeError:
            pass
        try:
            lines.append('Measures:' )
            for line in str(self.measuretree.absolute_dic).split(','):
                lines.append(line)
        except AttributeError:
            pass
        lines.append('')
        lines.append('Observations:')
        if len(self.observations.keys()) > 0:
            for line in str(self.observations).split(','):
                lines.append(line)
        else:
            lines.append('Not yet evaluated')
        return '\n'.join(lines)
    

    def simul(self, ablat=None, ablat_pos=0.): 
        """ The simulation main loop. 
        
        Keyword arguments:
        ------------------
        
        movie: bool, optional
            If True, runs _make_movie (default False) during the simulation.
            TODO: implement a public make_movie method that can be runned after
            the simulation 
        ablat: float, optional
            Timepoint at which ablation takes place. If None (default)
            no ablation is performed.

        """
        dt = self.KD.params['dt']
        kappa_c = self.KD.params['kappa_c']

        for time_point in range(1, self.num_steps):
            # Ablation test
            if ablat == time_point:
                self._ablation(time_point, pos=ablat_pos)
            # Anaphase transition ?
            self._anaphase_test(time_point)
            self.KD.one_step(time_point)
        self.KD.params['kappa_c'] = kappa_c
        delay_str = "delay = %2d seconds" % self.delay
        self.report.append(delay_str)
        for ch in self.KD.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()
            ch.cen_A.calc_toa()
            ch.cen_B.calc_toa()

    def _anaphase_test(self, time_point):
        """returns True if anaphase has been executed.
        At anaphase onset, set the cohesin spring constent to 0 and
        self.KD.anaphase to True.
        """
        t_A = int(self.KD.params['t_A'])
        dt = self.KD.params['dt']
        t = time_point * dt
        if self.KD.anaphase :
            return True
        if t >= t_A and self._plug_checkpoint():
            if self.delay == -1 :
                self.delay = t - t_A 
                #Then we just get rid of cohesin
                self.KD.params['kappa_c'] = 0.
                self.KD.calc_B()
                nb_mero = self._mero_checkpoint()
                if nb_mero:
                    s = ("There were %d merotelic MT at anaphase onset"
                         % nb_mero)
                    self.report.append(s)
                self.KD.anaphase = True
                return True
        return False

    def _ablation(self, time_point, pos = None):
        """Simulates a laser ablation: detaches all the kinetochores
        and sets the midzone stall force Fmz to 0

        Parameter:
        ----------
            pos: float or None, optional
                position of the laser beam within the spindle
        """
        if pos == None: pos = self.KD.spbR.pos
        if not self.KD.spbL.pos <= pos <= self.KD.spbR.pos:
            print 'Missed shot, same player play again!'
            return 
        self.KD.params['Fmz'] = 0.
        self.KD.params['k_a'] = 0.
        self.KD.params['k_d0'] = 0.

        for plugsite in self.KD.spindle.all_plugsites():
            if pos < plugsite.pos and plugsite.plug_state == - 1:
                plugsite.set_plug_state(0, time_point)
            elif pos > plugsite.pos and plugsite.plug_state == 1:
                plugsite.set_plug_state(0, time_point)

    def _plug_checkpoint(self):
        """If the spindle assembly checkpoint is active, returns True
        if all chromosomes are plugged by at least one kMT, False
        otherwise.

        """
        sac = self.KD.params['sac']
        if sac == 0:
            return True
        for ch in self.KD.chromosomes :
            if not ch.cen_A.is_attached() or not ch.cen_B.is_attached():
                ch.active_sac = 1
                return False
        return True

    def _mero_checkpoint(self):
        '''returns the total number of merotellic kT
        '''
        nb_mero = 0
        for ch in self.KD.chromosomes :
            if np.any(ch.erroneous()) :
                nb_mero += sum(ch.erroneous())
                #print "active checkpoint"
        return nb_mero

    def _mplate_checkpoint(self):
        """returns True if each kinetochore is in the proper half
        of the spindle

        """
        for ch in self.KD.chromosomes:
            ktR = ch.cen_A.pos
            ktL = ch.cen_B.pos
            if min(ktR, ktL) <= 0 and max(ktR, ktL) >= 0:
                return True 
        return True

    def evaluate(self):
        """ passes all the evaluations in eval_simul.py
        results are stored in the self.observations dictionnary

        """
        if len(self.KD.step_count) < 2:
            print "No simulation was runned ... exiting"
            return 0
        for name, function in evaluations.items():
            self.observations[name] = function(self.KD)

    def show_trajs(self, axes = None): 
        """ Plot the different trajectories

        """
        N = int(self.KD.params['N'])
        if axes == None:
            fig = plt.figure()
            axes = fig.gca()
        spbRtraj = self.KD.spbR.traj
        spbLtraj = self.KD.spbL.traj
        axes.plot(self.timelapse, spbRtraj, color='r', ls='-', lw=1)
        axes.plot(self.timelapse, spbLtraj, color='r', ls='-', lw=1)

        fmt_list = ['g-', 'b-', 'm-']
        for n in range(N):
            ch = self.KD.chromosomes[n]
            right_traj = ch.cen_A.traj
            left_traj = ch.cen_B.traj
            fmt = fmt_list[np.mod(n, 3)]
            line1 = axes.plot(self.timelapse, right_traj, fmt, alpha=0.5)
            line2 = axes.plot(self.timelapse, left_traj, fmt, alpha=0.5)
            if n == 0:
                line1[0].set_alpha(1.)
                line2[0].set_alpha(1.)
        axes.set_xlabel('Time (seconds)', fontsize = 'small')
        axes.set_ylabel(u'Distance from center (um)', fontsize = 'small')
        plt.show()

    def write_asoctave(self, xmlfname = "results.xml",
                       datafname = "oct_data.txt"):
        """Save results in a file compatible with the tools developped
        to caracterize real life movies under octave.

        Parameters
        ----------
        xmlfname : string
            name of the xml file where the paramters and observations
            are stored, defaults to ``results.xml``
        datafname : string
            name of the xml file where the trajectories coordinates
            are stored, defaults to ``oct_data.txt``

        Note that this will also call self.write_results(xmlfname).
        
        """
        self.write_results(xmlfname = xmlfname)
        spbRx = self.KD.spbR.traj
        ys = np.zeros(spbRx.shape)
        ts = np.arange(spbRx.shape[0])
        idxs = np.ones(spbRx.shape)
        out_array = np.vstack((spbRx, ys, ts, idxs))
        
        spbLx = self.KD.spbL.traj
        idxs += 1
        vstack = np.vstack((spbLx, ys, ts, idxs))
        out_array = np.append(out_array, vstack, axis=0)

        for ch in self.KD.chromosomes:
            idxs += 1
            chRx = ch.cen_A.traj
            vstack = np.vstack((chRx, ys, ts, idxs))
            out_array = np.append(out_array, vstack, axis=0)

            idxs += 1            
            chLx = ch.cen_B.traj
            vstack = np.vstack((chLx, ys, ts, idxs))
            out_array = np.append(out_array, vstack, axis = 0)
            
        dataout = file(datafname, 'w+')
        dataout.write("# desciptor: "+xmlfname+"\n")
        np.savetxt(dataout, out_array, delimiter=' ')
        
    def write_results(self, xmlfname = "results.xml", datafname = "data.npy"):
        """ Saves the results of the simulation in two files
        with the parameters, measures and observations in one file
        and the trajectories in the other.

        Keyword arguments:
        ------------------
        xmlfname : string, optional
            name of the xml file where parameters and observations
            will be written
        datafname : string, optional
            name of the file where the trajectories will be written
            file type is determined by the file suffix:
                 - *.npy : data are stored in numpy's binary format
                           (less portable but quite efficient)
                 - *.txt : simple text
                 - *.txt.gz : text files compressed transparently

        Any other suffix will be saved as plain text. Column index for
        each trajectory is an attribute of the corresponding element
        in the xml file.

        TODO : This function should be chopped off some how, it's messy

        """
        if not hasattr(self, 'observations'):
            self.evaluate()

        chromosomes = self.KD.chromosomes
        wavelist = []
        out = file(xmlfname, 'w+')
        out.write('<?xml version="1.0"?>\n')
        today = time.asctime()
        experiment = Element("experiment", date=today, datafile=datafname)
        experiment.append(self.paramtree.root)
        experiment.append(self.measuretree.root)

        #right SPB
        spbR = SubElement(experiment, "trajectory", name = "rightspb",
                          column='0', units='mu m')
        SubElement(spbR, "description").text="right spb trajectory"
        spbRtraj = np.array(self.KD.spbR.traj)
        wavelist.append(spbRtraj)

        #left SPB
        spbL = SubElement(experiment, "trajectory", name = "leftspb",
                          column='1', units='mu m')
        SubElement(spbL, "description").text="left spb trajectory"
        spbLtraj = np.array(self.KD.spbL.traj)
        wavelist.append(spbLtraj)

        col_num = 2
        #chromosomes
        for n, ch in enumerate(chromosomes):
            rch = SubElement(experiment, "trajectory", name="centromereA",
                             index = str(n), column=str(col_num), units='mu m')
            text = "chromosome %i centromere A trajectory" % n
            SubElement(rch, "description").text = text
            wavelist.append(ch.cen_A.traj)
            col_num += 1
            
            SubElement(experiment, "numbercorrect", name="centromereA",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.correct_history[:, 0])
            col_num += 1
            SubElement(experiment, "numbererroneous",
                       name="centromereA", index = str(n),
                       column=str(col_num))
            wavelist.append(ch.erroneous_history[:, 0])
            col_num += 1
            
            lch = SubElement(experiment, "trajectory", index=str(n), 
                             column=str(col_num), units='mu m')
            text = "chromosome %s left kinetochore trajectory" % n
            SubElement(lch, "description").text = text
            wavelist.append(np.array(ch.cen_B.traj))
            col_num += 1
            SubElement(experiment, "numbercorrect", name="centromereB",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.correct_history[:, 1])
            col_num += 1

            SubElement(experiment, "numbererroneous", name="centromereB",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.erroneous_history[:, 1])
            col_num += 1

            #Plug Sites
            for m, plugsite in enumerate(ch.cen_A.plugsites):
                SubElement(experiment, "trajectory", name="plugsite",
                           index=str((n, m)), cen_tag='A',
                           column=str(col_num), units='mu m')
                wavelist.append(plugsite.traj)
                col_num += 1

                SubElement(experiment, "state", name="plugsite",
                           index=str((n, m)), cen_tag='A',
                           column=str(col_num), units='')
                wavelist.append(plugsite.state_hist)
                col_num += 1
            
            for m, plugsite in enumerate(ch.cen_B.plugsites):

                SubElement(experiment, "trajectory", name="plugsite",
                           index = str((n, m)), cen_tag='B',
                           column=str(col_num), units='mu m')
                wavelist.append(plugsite.traj)
                col_num += 1

                SubElement(experiment, "state", name="plugsite",
                           index = str((n, m)), cen_tag='B',
                           column=str(col_num), units='')
                wavelist.append(plugsite.state_hist)
                col_num += 1

        #Observations
        obs_elem = SubElement(experiment, "observations")
        if not hasattr(self, 'observations'):
            self.evaluate()
        for key, val in self.observations.items():
            SubElement(obs_elem, key).text = str(val)

        #Now we write down the whole experiment XML element
        indent(experiment)
        out.write(tostring(experiment))
        out.close()
        #And the numbers in the file datafname
        dataout = file(datafname, 'w+')
        data = np.vstack(wavelist).T
        if datafname.endswith('.npy'):
            np.save(dataout, data)
        else:
            dataout.write("# desciptor: "+xmlfname+"\n")
            np.savetxt(dataout, data, delimiter=' ')
        print "Simulation saved to file %s " % xmlfname

    def _make_movie(self, t, imsize):
        """somehow deprecated"""
        try:
            os.mkdir("movie")
        except OSError:
            pass
        
        imname = "movie/spindle%03i.tif" % t
        if os.path.isfile(imname):
            os.remove(imname)

        imarray = np.zeros(imsize, np.uint8)

        xspbR = scale(self.KD.spbR.pos, imsize[1])
        yspbR = scale(0, imsize[0])
        xspbL = scale(self.KD.spbL.pos, imsize[1])
        yspbL = scale(0, imsize[0])
        imarray[yspbR, xspbR, 0] += 128
        imarray[yspbL, xspbL, 0] += 128

        chromosomes = self.KD.chromosomes
        for n, ch in enumerate(chromosomes):
            xchR = scale(ch.cen_A.pos, imsize[1])
            ychR = scale(0.1 * n - 0.1, imsize[0])
            xchL = scale(ch.cen_B.pos, imsize[1])
            ychL = scale(0.1 * n - 0.1, imsize[0])
            imarray[ychR, xchR, 1] += 42
            imarray[ychL, xchL, 1] += 42

        if 80 < t < 100:
            self.testim =  imarray
            print xchR, ychR, t
        im = image_fromarray(imarray, 'RGB')
        im.save(imname)
        del im

    def show_one(self, n = 0, fig = None):
        """ Shows chromosome n trajectory and plug state """
        dt = self.KD.params['dt']
        ch = self.KD.chromosomes[n]

        if fig == None:
            fig = plt.figure()
        fig.clear()
        
        #fig.add_subplot(312)
        gridspec = plt.GridSpec(5,1)
        subplotspec = gridspec.new_subplotspec((1,0), rowspan=3)
        traj_ax = fig.add_subplot(subplotspec)
        traj_ax.plot(self.timelapse, ch.cen_A.traj, 'g', lw=2, alpha=0.5)
        traj_ax.plot(self.timelapse, ch.cen_B.traj, 'purple', lw=2, alpha=0.5)
        traj_ax.plot(self.timelapse, self.KD.spbR.traj, 'k')
        traj_ax.plot(self.timelapse, self.KD.spbL.traj, 'k')
        for plugsite in ch.cen_A.plugsites:
            traj_ax.plot(self.timelapse, plugsite.traj, 'g')
        for plugsite in ch.cen_B.plugsites:
            traj_ax.plot(self.timelapse, plugsite.traj, 'purple')
        traj_ax.set_xticks([], '')

        erroneous_hist = ch.erroneous_history
        correct_hist = ch.correct_history

        subplotspec = gridspec.new_subplotspec((0,0), rowspan=1)
        ax = fig.add_subplot(subplotspec, sharex=traj_ax)
        ax.plot(self.timelapse, erroneous_hist[:, 0], 'r',
                label='number of erroneoustellic MTs')
        ax.plot(self.timelapse, correct_hist[:, 0], 'g',
                label='number of correct MTs')
        ax.axis((0, self.num_steps*dt, -0.5, 4.5))
        ax.set_xticks([], '')

        subplotspec = gridspec.new_subplotspec((4,0), rowspan=1)
        ax = fig.add_subplot(subplotspec, sharex=traj_ax)
        ax.plot(self.timelapse, erroneous_hist[:, 1], 'r',
                label='number of erroneous MTs')
        ax.plot(self.timelapse, correct_hist[:, 1], 'purple',
                label='number of correct MTs')
        ax.axis((0, self.num_steps*dt, -0.5, 4.5))
        plt.show()
        
    def get_ch(self, n = 0):
        return self.KD.chromosomes[n]

    def get_2Dtraj_list(self, ang_noise=5e-3, pos_noise=1e-2,
                        perp_distance=0.3,
                        cen=None, cen_sigma=0.6):
        
        """returns a list of 2D trajectories, adding a random mouvement
        to the whole spindle.

        Keyword arguments:
        ------------------
        ang_noise : float
            angular variation in radians per seconds
        pos_noise : float
            center of mass displacement in um/s
        perp_distance: float
            perpendicular distance between the kinetochore pairs 
        cen: int or None
            chromosome index for which a centromere trajectory is
            simulated. Be aware that this adds a random coiled coil
            movement around the kinetochore position.
            If cen is None returns all 6 trajectories plus the SPBs
        cen_sigma: float
            sets the amplitude of the simulated random coiled coil
            movement of the cen marker with respect to the chromosome

        Returns:
        --------
        trajectories: ndarray
            Trajectories is an array of shape (n, 2, num_steps),
            where n is the number of trajectories
        """

        dt = self.KD.params['dt']
        N =  self.KD.params['N']
        ang_noise *= dt # rad.s^-1
        pos_noise *=  dt # um.s^-1

        n_traj = 2 + 2 * N if cen is None else 4
        trajs_shape = (n_traj, 2, self.num_steps)
        trajectories = np.zeros(trajs_shape)

        trajectories[0, 0, :] = self.KD.spbR.traj
        trajectories[1, 0, :] = self.KD.spbL.traj

        if cen is not None:
            ch = self.get_ch(cen)
            trajectories[2, 0, :] = ch.cen_A.traj
            trajectories[3, 0, :] = ch.cen_B.traj
            trajectories[2:, ...] += normal(0, scale=cen_sigma,
                                           size=(2, 2, self.num_steps))
        else:
            for n, ch in enumerate(self.KD.chromosomes):
                trajectories[2 + 2 * n, 0, :] = ch.cen_A.traj
                trajectories[2 + 2 * n, 1, :] += (1 - n) * perp_distance
                trajectories[2 + 2 * n + 1, 0, :] = ch.cen_B.traj
                trajectories[2 + 2 * n + 1, 1, :] += (1 - n) * perp_distance
        # TODO: block needs to be vectorized GG march 2012
        xcs = []
        ycs = []
        rots = []
        xd, yd, thetad = (0.,)*3
        for n in range(self.num_steps):
            xd += normal(0, scale = pos_noise)
            yd += normal(0, scale = pos_noise)
            thetad += normal(0, scale = ang_noise)
            xcs.append(xd)
            ycs.append(yd)
            rots.append([[np.cos(thetad), np.sin(thetad)],
                         [np.cos(thetad), - np.sin(thetad)]])
        # TODO: block needs to be vectorized GG march 2012        
        for traj in trajectories:
            traj += np.vstack((np.array(ycs), np.array(xcs)))
            n = 0
            for pos, rot in zip(traj, rots):
                new_pos = np.dot(pos, rot)
                traj[n] = new_pos
                n += 1
        # TODO: implement a general vectorized Brownian motion
        return trajectories

    def get_3Dtraj_list(self, ang_noise=5e-2, pos_noise=3e-2,
                        radial_distance=0.3, cen=None, cen_sigma=0.6):
        """returns a list of 3D trajectories, adding a random mouvement [1]_
        to the whole spindle. 

        Keyword arguments:
        ------------------
        ang_noise : float
            angular variation in radians per seconds, sets the standard
            deviation of the spindle axis angle in um/s.

        pos_noise : float
            std. dev of the spindle center displacement in um/s 
        radial_distance: float
            radial distance of the centromeres to the spindle axis.
            Here, chromosomes are distributed evenly around the spindle
            axis. 
        cen: int or None
            chromosome index for which a centromere trajectory is
            simulated. Be aware that this adds a random coiled coil
            movement around the kinetochore position.
            If cen is None returns all 6 trajectories plus the SPBs
        cen_sigma: float or None
            if cen is not None, sets the 

        Returns:
        --------

        trajectories: a list of ndarrays


        .. [1] By adding a normaly distributed noise with a Gaussian
               distribution see numpy.random.normal for further details
        """

        dt = self.KD.params['dt']
        N = self.KD.params['N']
        trajectories = []
        ang_noise *= dt # rad.s^-1
        pos_noise *=  dt # um.s^-1
        
        n_traj = 2 + 2 * N if cen is None else 4
        trajs_shape = (n_traj, 3, self.num_steps)
        trajectories = np.zeros(trajs_shape)
        trajectories[0, 0, :] = self.KD.spbR.traj
        trajectories[1, 0, :] = self.KD.spbL.traj

        if cen is not None:
            ch = self.get_ch(cen)
            trajectories[2, 0, :] = ch.cen_A.traj
            trajectories[3, 0, :] = ch.cen_B.traj
            trajectories[2:, ...] += normal(0, scale=cen_sigma,
                                           size=(2, 2, self.num_steps))
        else:
            for n, ch in enumerate(self.KD.chromsomes):
                trajectories[2 + 2 * n, 0, :] = ch.cen_A.traj
                phi = (n - 1) * 2 * np.pi / N
                radial_y = radial_distance * np.cos(phi)
                trajectories[2 + 2 * n, 1, :] += radial_y
                trajectories[2 + 2 * n + 1, 1, :] += radial_y
                radial_z = radial_distance * np.sin(phi)                
                trajectories[2 + 2 * n + 1, 2, :] += radial_z
        
        xcs, ycs, zcs = np.zeros((3, self.num_steps))
        #Fix the x axis, rotate around the two others
        xy_rots = np.zeros((3, 3, self.num_steps))
        zx_rots = np.zeros((3, 3, self.num_steps))

        thetad, phid = 0., 0.
        for n in range():
            thetad += normal(0, scale = ang_noise)
            phid += normal(0, scale = ang_noise)
            xcs[n] += normal(0, scale = pos_noise)
            ycs[n] += normal(0, scale = pos_noise)
            zcs[n] += normal(0, scale = pos_noise)
            xy_rots[:, :, n] = [[np.cos(thetad), - np.sin(thetad), 0], 
                                [np.sin(thetad), np.cos(thetad), 0],
                                [0, 0, 1]]
            zx_rots[:, :, n] = [[np.cos(phid), 0, - np.sin(phid)], 
                                [0, 1, 0],
                                [np.sin(phid), 0, np.cos(phid)]]
        for traj in trajectories:
            traj += np.vstack((xcs, ycs, zcs))
            for n, pos in enumerate(traj.T):
                tmp_pos = np.dot(pos, xy_rots[:, :, n])
                new_pos = np.dot(tmp_pos, zx_rots[:, :, n])
                traj[:, n] = new_pos
        return trajectories


def reduce_params(paramtree, measuretree):
    """ This functions changes the paramters so that 
    the dynamical characteristics complies with the measures [1]_.

    input:
    -----
    paramtree : xml_handler.ParamTree instance
    measuretree : xml_handler.MeasureTree

    paramtree is modified in place


    .. [1] G. Gay, T.Courthéoux, C. Reyes, S. Tournier, Y. Gachet.
           J. Cell Biol 2012 http://dx.doi.org/10.1083/jcb.201107124

    """
    params = paramtree.absolute_dic
    measures = measuretree.absolute_dic
    try:
        poleward_speed = measures['poleward_speed'] 
        metaph_rate = measures['metaph_rate']  
        anaph_rate = measures['anaph_rate']   
        mean_metaph_k_dist = measures['mean_metaph_k_dist']
        max_metaph_k_dist = measures['max_metaph_k_dist']        
        outer_inner_dist = measures['oi_dist']
        tau_k = measures['tau_k'] 
        tau_c = measures['tau_c'] 
        obs_d0 = measures['obs_d0']
    except KeyError:
        print ("The measures dictionary should contain"
               "at least the following keys: ")
        print MEASURES.keys()
        return 0

    k_a = params['k_a'] # 'free' attachement event frequency
    k_d0 = params['k_d0'] # 'free' detachement event frequency
    d_alpha = params['d_alpha']
    N = int(params['N'])
    Mk = int(params['Mk'])
    kappa_k = params['kappa_k']
    Fk = params['Fk']

    #Let's go for the direct relations
    d0 = params['d0'] = obs_d0
    Vk = params['Vk'] = poleward_speed
    Vmz = params['Vmz'] = anaph_rate
    #Aurora modifies fd
    if d_alpha != 0:
        k_d_eff = k_a * d_alpha / mean_metaph_k_dist
    else:
        print "Warning; things don't go well without Aurora "
        k_d_eff = k_d0
    
    # alpha_mean = float(mean_attachment(k_a/fd_eff) / Mk)
    alpha_mean = 1/(1 + k_d_eff/k_a)
    #Take metaphase kt pair distance as the maximum one
    kappa_c = Fk * Mk / ( max_metaph_k_dist - d0 )
    params['kappa_c'] = kappa_c

    #kop = alpha_mean * ( 1 + metaph_rate/2 ) / ( outer_inner_dist )
    kappa_k = Fk * Mk / ( 2 * outer_inner_dist )
    params['kappa_k'] = kappa_k
    #Ensure we have sufficientely small time steps
    dt = params['dt']
    params['dt'] = min( tau_c/4., tau_k/4., params['dt'])
    if params['dt'] != dt:
        print 'Time step changed' 

    mus = params['mus']
    Fmz =  ( Fk * N * Mk * alpha_mean * (1 +  metaph_rate / ( 2 * Vk ))
             + mus * metaph_rate / 2.  ) / (1 -  metaph_rate / Vmz )
    params['Fmz'] = Fmz
    muc = ( tau_c * kappa_c )
    params['muc'] = muc
    muk = ( tau_k * kappa_k )
    params['muk'] = muk 
    for key, val in params.items():
        paramtree.change_dic(key, val, write = False, verbose = False)

def get_fromfile(xmlfname = "results.xml"):
    """Creates a simul_spindle.Metaphase from a XML file.

    Returns:
    --------
        metaphase: a Metaphase instance

    """
    restree = ResultTree(xmlfname)
    param_root = restree.root.find('parameters')
    paramtree = ParamTree(root = param_root)
    params = paramtree.relative_dic
    measure_root = restree.root.find('measures')
    measuretree = ParamTree(root = measure_root,
                            adimentionalized = False)
    metaphase = Metaphase(paramtree = paramtree,
                          measuretree = measuretree)
    
    traj_matrix = restree.get_all_trajs()
    correct_matrix = restree.get_all_correct()
    erroneous_matrix = restree.get_all_erroneous()
    state_hist_matrix = restree.get_all_plug_state()
    KD = KinetoDynamics(params)
    KD.spbR.traj = traj_matrix[:, 0]
    KD.spbL.traj = traj_matrix[:, 1]
    Mk = int(params['Mk'])
    col_num = 2
    state_num = 0
    for n, ch in enumerate(KD.chromosomes) :
        ch.cen_A.traj = traj_matrix[:, col_num]
        col_num += 1
        ch.cen_B.traj = traj_matrix[:, col_num]
        col_num += 1
        ch.erroneous_history = (erroneous_matrix[:, n*2 : n*2 + 2])
        ch.correct_history = (correct_matrix[:, n*2 : n*2 + 2])
        for plugsite in ch.cen_A.plugsites:
            plugsite.traj = traj_matrix[:, col_num]
            col_num += 1
            plugsite.state_hist = state_hist_matrix[:, state_num]
            state_num += 1
        for plugsite in ch.cen_B.plugsites:
            plugsite.traj = traj_matrix[:, col_num]
            col_num += 1
            plugsite.state_hist = state_hist_matrix[:, state_num]
            state_num += 1
    metaphase.KD = KD

    return metaphase

def scale(x, size, pix_size = 0.0645):
    '''
    Scale the position x on a line of size "size" from microns to pixels
    origin is put in the center of the line
    '''
    return int(x / pix_size + size / 2)
    
def circ_box(xy, rad):
    return (xy[0] - rad, xy[1] - rad , xy[0] + rad , xy[1] + rad)
    
    
