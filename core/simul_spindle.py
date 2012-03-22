#!/usr/bin/python
# -*- coding: utf-8 -*-

### Adapted from Civelekoglu-Scholey et al. Biophys.J 90(11) 2006
### doi: 10.1529/biophysj.105.078691
### More details in Gay et al. JCB 2012
import numpy as np

import time
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy import linalg
from numpy.random import normal
from xml.etree.ElementTree import Element, SubElement, tostring
from Image import fromarray as image_fromarray

import pyximport
pyximport.install()

## local imports
from .spindle_dynamics import KinetoDynamics
from .xml_handler import ParamTree, indent, ResultTree
import ..analysis.eval_simul as eval_simul

__all__ = ["Metaphase", "reduce_params", "paramfile", "measurefile", "get_fromfile"]

try:
    paramfile = os.path.join(os.path.dirname(__file__),
                             'default', 'params.xml')
    measurefile = os.path.join(os.path.dirname(__file__),
                             'default', 'measures.xml')
except NameError:
    paramfile = 'default/params.xml'
    measurefile = 'default/measures.xml'
    
measuretree = ParamTree(measurefile, adimentionalized = False)
MEASURES = measuretree.absolute_dic


class Metaphase(object):
    """An instance of the Metaphase class is a wrapper around a whole simulation.

    Typical usage :
    ---------------

    >>> from kt_simul.simul_spindle import Metaphase
    >>> m = Metaphase()
    >>> m.simul()
    >>> m.show_trajs()
    >>> m.write_results('examples/docstring_results.xml',
                        'examples/docstring_data.npy')

    From already ran simulations:
    -----------------------------
    
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
    write_asoctave() : save the trajectories with a different txt organisation
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
                 paramfile=paramfile, measurefile=measurefile,
                 plug='random', reduce_p=True):

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

        plug : string or None
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
        params['Vk'] = self.paramtree.absolute_dic['Vk']
        params['Fk'] = self.paramtree.absolute_dic['Fk']
        params['dt'] = self.paramtree.absolute_dic['dt']
            
        self.KD = KinetoDynamics(params, plug=plug)
        dt = self.paramtree.absolute_dic['dt']
        duration = self.paramtree.absolute_dic['span']
        self.num_steps = int(duration/dt)
        self.KD.anaphase = False
        self.timelapse = np.arange(0, duration + dt, dt)
        self.report = []
        self.delay = -1

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
        try:
            lines.append('')
            lines.append('Observations:')
            for line in str(self.observations).split(','):
                lines.append(line)
        except AttributeError:
            lines.append('Not yet evaluated')
            pass
        return '\n'.join(lines)
    
    def _one_step(self):
        if not self.KD.anaphase:
            self.KD.plug_unplug()
        A = self.KD.calcA()
        b = - self.KD.calcb()
        speeds = linalg.solve(A,b)
        self.KD.position_update(speeds)

    def simul(self, movie = False, ablat = None): 
        """ The simulation main loop. 
        
        Keyword arguments:
        ------------------
        
        movie: bool, optional
            If True, runs _make_movie (default False) during the simulation.
            TODO: implement a public make_movie method that can be runned after the simulation 
        ablat: float, optional
            Timepoint at which ablation takes place. If None (default)
            no ablation is performed.

        """
        dt = self.KD.params['dt']
        kappa_c = self.KD.params['kappa_c']

        for t in self.timelapse[1:]:
            # Ablation test
            if ablat is not None and ablat == t:
                self._ablation(pos = self.KD.spbL.pos)
            # Anaphase transition ?
            if self._anaphase_test(t):
                self._toa_test(t)
            self._one_step()
            
            if movie and t/dt % 15 == 0: #One picture every 15 time point
                self._make_movie(t, (50, 100,  3))

        self.KD.params['kappa_c'] = kappa_c
        self.KD.delay = self.delay - 1
        s = "delay = %2d seconds" % self.delay
        self.report.append(s)

    def _toa_test(self, t):
        for n, ch in self.KD.chromosomes.items():    
            if ch.anaphase_switch[0] > 0 and ch.right_toa == 0:
                ch.right_toa = t
                s = ("Right kt of chromosome %d reached the pole\n"
                     "at time %.2f" %(n,t))
                self.report.append(s)
            if ch.anaphase_switch[1] > 0 and ch.left_toa == 0:
                ch.left_toa = t
                s = ("Left kt of chromosome %d reached the pole\n"
                     "at time %.2f" %(n,t))
                self.report.append(s)

    def _anaphase_test(self, t):
        t_A = int(self.KD.params['t_A'])
        if self.KD.anaphase:
            return True
        if t_A <= t and self._plug_checkpoint():
            if self.delay == -1 :
                self.delay = t - t_A 
                #Then we just get rid of cohesin
                self.KD.params['kappa_c'] = 0.
                self.KD.B_mat = self.KD.write_B()
                self.nb_mero = self._mero_checkpoint()
                if self.nb_mero:
                    s = ("There were %d merotelic MT at anaphase onset"
                         %self.nb_mero)
                    self.report.append(s)
                self.KD.anaphase = True
                self.KD.test_anaphase_switch()
                return True
        return False

    def _ablation(self, pos = None):
        if pos == None: pos = self.KD.spbR.pos
        if not self.KD.spbL.pos <= pos <= self.KD.spbR.pos:
            print 'Missed shot, same player play again!'
            return 

        self.KD.params['Fmz'] = 0.
        self.KD.params['k_a'] = 0.
        self.KD.params['k_d0'] = 0.
        for ch in self.KD.chromosomes.values():
            (right_pluged, left_pluged) = ch.pluged()
            (right_mero, left_mero) = ch.mero()
            if pos > ch.rightpos:
                for rplug in ch.rplugs.values():
                    rplug.plug = min(0, rplug.plug)
                for lplug in ch.lplugs.values():
                    lplug.plug = max(0, lplug.plug)
                #ch.pluged = (right_pluged, left_pluged)
                #ch.mero = (0, 0)

            elif ch.rightpos >= pos > ch.leftpos:
                for rplug in ch.rplugs.values():
                    rplug.plug = max(0, rplug.plug)
                for lplug in ch.lplugs.values():
                    lplug.plug = max(0, lplug.plug)
                #ch.pluged = (0, left_pluged)
                #ch.mero = (right_mero, 0)

            elif pos < ch.leftpos:
                for lplug in ch.lplugs.values():
                    lplug.plug = min(0, lplug.plug)
                for rplug in ch.rplugs.values():
                    rplug.plug = max(0, rplug.plug)

            (right_pluged, left_pluged) = ch.pluged()
            (right_mero, left_mero) = ch.mero()
            print (right_pluged, left_pluged)
            print (right_mero, left_mero)
                #ch.pluged = (right_pluged, 0)
                #ch.mero = (0,left_mero)

    def _plug_checkpoint(self):
        """If the spindle assembly checkpoint is active, returns True
        if all chromosomes are pluged by at least one kMT, False
        otherwise.

        """
        sac = self.KD.params['sac']
        if sac == 0:
             return True
        
        for ch in self.KD.chromosomes.values() :
            if not np.all(ch.pluged()) :
                ch.active_sac = 1
                #print "active checkpoint"
                return False
            else:
                ch.active_sac = 0
                
        return True

    def _mero_checkpoint(self):
        '''returns the total number of merotellic kT
        '''
        nb_mero = 0
        for ch in self.KD.chromosomes.values() :
            if np.any(ch.mero) :
                nb_mero += sum(ch.mero())
                #print "active checkpoint"
        return nb_mero

    def _mplate_checkpoint(self):
        """returns True if each kinetochore is in the proper half
        of the spindle
        """
        for ch in self.KD.chromosomes.values():
            ktR = ch.rightpos
            ktL = ch.leftpos
            if ktR <= 0 or ktL >= 0:
                return True 
        return True

    def evaluate(self):
        """ passes all the evaluations in eval_simul.py
        results are stored in the self.observations dictionnary

        """
        if len(self.KD.spbR.traj) < 2:
            print "No simulation was ran ... exiting"
            return 0
        
        observations = {'anaphase_rate' : eval_simul.anaphase_rate(self.KD),
                        'metaph_rate' : eval_simul.metaph_rate(self.KD),
                        'mean_metaph_k_dist': eval_simul.metaph_kineto_dist(self.KD),
                        'pitch' : eval_simul.auto_corel(self.KD, smooth = 5.),
                        'poleward_speed' : eval_simul.poleward_speed(self.KD),
                        'kt_rms_speed' : eval_simul.kt_rms_speed(self.KD),
                        'times_of_arrival' : eval_simul.time_of_arrival(self.KD),
                        'pluged_stats' : eval_simul.pluged_stats(self.KD)}


    def show_trajs(self, axes = None): 
        """ Plot the different trajectories

        """
        N = int(self.KD.params['N'])
        if axes == None:
            fig = plt.figure()
            axes = fig.gca()

        spbRtraj = np.array(self.KD.spbR.traj)
        spbLtraj = np.array(self.KD.spbL.traj)
        try:
            axes.plot(self.timelapse, spbRtraj, color='r', ls='-', lw=1)
            axes.plot(self.timelapse, spbLtraj, color='r', ls='-', lw=1)
        except ValueError:
            print 'please run the simulation before'
            return 0

        fmt_list = ['g-', 'b-', 'm-']
        for n in range(N):
            ch = self.KD.chromosomes[n]
            right_traj = np.array(ch.righttraj)
            left_traj = np.array(ch.lefttraj)
            fmt = fmt_list[np.mod(n, 3)]
            line1 = axes.plot(self.timelapse,right_traj, fmt, alpha=0.5)
            line2 = axes.plot(self.timelapse,left_traj, fmt, alpha=0.5)
            if n == 0:
                line1[0].set_alpha(1.)
                line2[0].set_alpha(1.)
            # if n == 0:
            #     axes.plot(self.timelapse[::int(15/dt)], right_traj[::int(15/dt)],
            #               mec, mfc = 'None')
            #     axes.plot(self.timelapse[::int(15/dt)], left_traj[::int(15/dt)],
            #               'go', mfc = 'None')
                
        #axes.axis([0, self.num_steps * dt, -2.2, 2.2 ])
        axes.set_xlabel('Time (seconds)', fontsize = 'small')
        axes.set_ylabel(u'Distance from center (um)', fontsize = 'small')
        plt.show()

    def write_asoctave(self, xmlfname = "results.xml",
                       datafname = "oct_datas.txt"):
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

        spbRx = np.array(self.KD.spbR.traj)
        ys = np.zeros(spbRx.shape)
        ts = np.arange(spbRx.shape[0])
        idxs = np.ones(spbRx.shape)
        out_array = np.vstack((spbRx, ys, ts, idxs))
        
        spbLx = np.array(self.KD.spbL.traj)
        idxs += 1
        vstack = np.vstack((spbLx, ys, ts, idxs))
        out_array = np.append(out_array, vstack, axis=0)

        for ch in self.KD.chromosomes.values():
            idxs += 1
            chRx = np.array(ch.righttraj)
            vstack = np.vstack((chRx, ys, ts, idxs))
            out_array = np.append(out_array, vstack, axis=0)

            idxs += 1            
            chLx = np.array(ch.lefttraj)
            vstack = np.vstack((chLx, ys, ts, idxs))
            out_array = np.append(out_array, vstack, axis = 0)
            
        dataout = file(datafname, 'w+')
        dataout.write("# desciptor: "+xmlfname+"\n")
        np.savetxt(dataout, out_array, delimiter=' ')
        
    def write_results(self, xmlfname = "results.xml", datafname = "datas.npy"):
        """ Saves the results of the simulation in two files with the parameters,
        the measures and observations in one file and the trajectories in
        an other file

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

        """
        if not hasattr(self, 'observations'):
            self.evaluate()

        chromosomes = self.KD.chromosomes.values()
        wavelist=[]
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
            # Adding pluged_history and mero_history
            chp = np.array(ch.pluged_history)
            rp = chp[:, 0]
            lp = chp[:, 1]
            chm = np.array(ch.mero_history)
            rm = chm[:, 0]
            lm = chm[:, 1]

            rch = SubElement(experiment, "trajectory", name="rightkineto",
                             index = str(n), column=str(col_num), units='mu m')
            text = "chromosome %i right kinetochore trajectory" %n
            SubElement(rch, "description").text = text
            wavelist.append(np.array(ch.righttraj))
            col_num += 1
            
            SubElement(experiment, "numberpluged", name="rightkineto",
                       index=str(n), column=str(col_num))
            wavelist.append(rp)
            col_num += 1
            SubElement(experiment, "numbermero",
                       name="rightkineto", index = str(n),
                       column=str(col_num))
            wavelist.append(rm)
            col_num += 1
            
            lch = SubElement(experiment, "trajectory", index=str(n), 
                             column=str(col_num), units='mu m')
            text = "chromosome %s left kinetochore trajectory" %n
            SubElement(lch, "description").text = text
            wavelist.append(array(ch.lefttraj))
            col_num += 1
            SubElement(experiment, "numberpluged", name="leftkineto",
                       index=str(n), column=str(col_num))
            wavelist.append(lp)
            col_num += 1
            SubElement(experiment, "numbermero", name="leftkineto",
                       index=str(n), column=str(col_num))
            wavelist.append(lm)
            col_num += 1

            #Plug Sites
            for m, rplug in enumerate(ch.rplugs.values()):
                rpt_traj = np.array(rplug.traj)
                SubElement(experiment, "trajectory",
                           name="rightplugsite", index = str((n, m)),
                           column=str(col_num), units='mu m')
                wavelist.append(rpt_traj); col_num += 1
                rpt_plug = np.array(rplug.state_hist)
                SubElement(experiment, "state",
                           name="rightplugsite", index = str((n, m)),
                           column=str(col_num), units='')
                wavelist.append(rpt_plug); col_num += 1
            
            for m, lplug in enumerate(ch.lplugs.values()):
                lpt_traj = np.array(lplug.traj)
                SubElement(experiment, "trajectory",
                           name="leftplugsite", index = str((n, m)),
                           column=str(col_num), units='mu m')
                wavelist.append(lpt_traj); col_num += 1
                lpt_plug = np.array(lplug.state_hist)
                SubElement(experiment, "state",
                           name="leftplugsite", index = str((n, m)),
                           column=str(col_num), units='')
                wavelist.append(lpt_plug); col_num += 1

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
        data = np.vstack(wavelist)
        if datafname.endswith('.npy'):
            np.save(dataout, data)
        else:
            dataout.write("# desciptor: "+xmlfname+"\n")
            np.savetxt(dataout, data, delimiter=' ')
        print "Simulation saved to file %s " %datafname

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

        chromosomes = self.KD.chromosomes.values()
        for n, ch in enumerate(chromosomes):
            xchR = scale(ch.rightpos, imsize[1])
            ychR = scale(0.1 * n - 0.1, imsize[0])
            xchL = scale(ch.leftpos, imsize[1])
            ychL = scale(0.1 * n - 0.1, imsize[0])
            imarray[ychR, xchR, 1] += 42
            imarray[ychL, xchL, 1] += 42

        if 80< t <100:
            self.testim =  imarray
            print xchR, ychR, t
        im = image_fromarray(imarray, 'RGB')
        im.save(imname)
        del im

    def show_one(self, n = 0, fig = None):
        """ Shows one chromosome trajectory and plug state """
        dt = self.KD.params['dt']
        if fig == None:
            fig = plt.figure()

        ch = self.KD.chromosomes[n]
        fig.clear()
        fig.add_subplot(312)
        ax = fig.gca()
        ax.plot(self.timelapse, ch.righttraj, 'b', lw=1.5)
        ax.plot(self.timelapse, ch.lefttraj, 'k', lw=1.5)
        rps = ch.rplugs
        for ps in rps.values():
            ax.plot(self.timelapse, ps.traj, 'g')
        lps = ch.lplugs
        for ps in lps.values():
            ax.plot(self.timelapse, ps.traj, 'purple')
        mero_hist = np.asarray(ch.mero_history)
        pluged_hist = np.asarray(ch.pluged_history)
        
        fig.add_subplot(311)
        ax = fig.gca()   
        ax.plot(self.timelapse, mero_hist[:,0], 'r',
                label='number of merotellic MTs')
        ax.plot(self.timelapse, pluged_hist[:,0], 'g',
                label='number of pluged MTs')
        ax.axis((0, self.num_steps*dt, -0.5, 4.5))

        fig.add_subplot(313)
        ax = fig.gca()
        ax.plot(self.timelapse, mero_hist[:,1], 'r',
                label='number of merotelic MTs')
        ax.plot(self.timelapse, pluged_hist[:,1], 'g',
                label='number of pluged MTs')
        ax.axis((0, self.num_steps*dt, -0.5, 4.5))
        plt.show()

    def show_one_article(self, n = 0, fig = None):
        """ Shows one chromosome trajectory and plug state,
        alternate layout"""        

        dt = self.KD.params['dt']
        if fig == None:
            fig = plt.figure()

        ch = self.KD.chromosomes[n]
        fig.clear()
        ax_traj = fig.add_subplot(312)
        spbRtraj = np.array(self.KD.spbR.traj)
        spbLtraj = np.array(self.KD.spbL.traj)

        ax_traj.plot(self.timelapse, spbRtraj, 'r-')
        ax_traj.plot(self.timelapse, spbLtraj, 'r-')
        ax_traj.plot(self.timelapse, ch.righttraj, 'g-', lw=1., alpha=1.)
        ax_traj.plot(self.timelapse, ch.lefttraj, 'g-', lw=1., alpha=1.)
        rps = ch.rplugs
        for ps in rps.values():
            ax_traj.plot(self.timelapse, ps.traj,
                         'purple', lw=0.5, alpha=0.5)
        lps = ch.lplugs
        for ps in lps.values():
            ax_traj.plot(self.timelapse, ps.traj,
                         'purple', lw=0.5, alpha=0.5)

        ax_traj.axis((0, self.num_steps*dt, -2, 2))
        ax_traj.yticks(range(-2,4,2))
        ax_traj.xticks(range(0, 850 , 240), ['0', '4', '8', '12'])
        ax_traj.ylabel(u'Position (um)', fontsize='small')
        ax_traj.xlabel('Time (min)', fontsize='small')
        
        mero_hist = np.array(ch.mero_history)
        pluged_hist = np.array(ch.pluged_history)
        ax_plug = fig.add_subplot(311)
        ax_plug.plot(self.timelapse, mero_hist[:,0], 'r-',
                     label= 'Number of merotellic kMTs')
        ax_plug.plot(self.timelapse, pluged_hist[:,0], 'g-',
                     label= 'Number of synthelic kMTs')
        ax_plug.axis((0, self.num_steps*dt, -0.5, 4.5))
        ax_plug.xticks(range(0, 850 , 240), (''))
        ax_plug.ylabel('kMTs state', fontsize = 'small')
        ax_plug.yticks(np.arange(0,5,2))

        ax_plug2 = fig.add_subplot(313)
        ax_plug2.plot(self.timelapse, mero_hist[:,1], 'r-',
                      label= 'number of merotellic MTs')
        ax_plug2.plot(self.timelapse, pluged_hist[:,1], 'g-',
                      label= 'number of synthelic MTs')
        ax_plug2.axis((0, self.num_steps*dt, -0.5, 4.5))
        ax_plug2.xticks(range(0, 850 , 240), (''))
        ax_plug2.ylabel('kMTs state', fontsize = 'small')
        plt.show()
        return fig
        
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
        trajectories = np.zeros(t_shape)

        trajectories[0, 0,:] = self.KD.spbR.traj
        trajectories[1, 0,:] = self.KD.spbL.traj

        if cen is not None:
            ch = self.get_ch(cen)
            trajectories[2, 0,:] = ch.righttraj
            trajectories[3, 0,:] = ch.lefttraj
            trajectories[2:,...] += normal(0, scale=cen_sigma,
                                           size=(2, 2, traj_length))
        else:
            for n, ch in enumerate(self.KD.chromosomes):
                trajectories[2 + 2 * n, 0,:] = ch.righttraj
                trajectories[2 + 2 * n, 1,:] += (1 - n) * perp_distance
                trajectories[2 + 2 * n + 1, 0,:] = ch.lefttraj
                trajectories[2 + 2 * n + 1, 1,:] += (1 - n) * perp_distance

        # This block needs to be vectorized GG march 2012
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
        # This block needs to be vectorized GG march 2012        
        for traj in trajectories:
            traj += np.vstack((np.array(ycs), np.array(xcs)))
            n = 0
            for pos, rot in zip(traj, rots):
                new_pos = np.dot(pos, rot)
                traj[n] = new_pos
                n += 1
        # TODO: implement a general vectorized Brownian motion

        return trajectories

    def get_3Dtraj_list(self, ang_noise=5e-2,
                        pos_noise=3e-2, cen=None
                        radial_distance=0.3):
        
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
        trajectories = []
        ang_noise *= dt # rad.s^-1
        pos_noise *=  dt # um.s^-1

        
        n_traj = 2 + 2 * N if cen is None else 4
        trajs_shape = (n_traj, 3, self.num_steps)
        trajectories = np.zeros(t_shape)
        

        trajectories[0, 0,:] = self.KD.spbR.traj
        trajectories[1, 0,:] = self.KD.spbL.traj

        if cen is not None:
            ch = self.get_ch(cen)
            trajectories[2, 0,:] = ch.righttraj
            trajectories[3, 0,:] = ch.lefttraj
            trajectories[2:,...] += normal(0, scale=cen_sigma,
                                           size=(2, 2, traj_length))
        else:
            for n, ch in enumerate(self.KD.chromsomes):
                trajectories[2 + 2 * n, 0,:] = ch.righttraj
                phi = (n - 1) * 2 * np.pi / N
                radial_y = radial_distance * np.cos(phi)
                trajectories[2 + 2 * n, 1,:] += radial_y
                trajectories[2 + 2 * n + 1, 1,:] += radial_y
                radial_z = radial_distance * np.sin(phi)                
                trajectories[2 + 2 * n + 1, 2,:] += radial_z
        
        xcs, ycs, zcs = np.zeros((3, xspbR.size))

        #Fix the x axis, rotate around the two others
        xy_rots = np.zeros((3,3,xspbR.size))
        zx_rots = np.zeros((3,3,xspbR.size))

        xd, yd, zd, thetad, phid = (0.,)*5
        for n in range():
            thetad += normal(0, scale = ang_noise)
            phid += normal(0, scale = ang_noise)
            xcs[n] += normal(0, scale = pos_noise)
            ycs[n] += normal(0, scale = pos_noise)
            zcs[n] += normal(0, scale = pos_noise)
            xy_rots[:,:,n] = [[cos(thetad), - np.sin(thetad), 0], 
                              [np.sin(thetad), np.cos(thetad), 0],
                              [0, 0, 1]]

            zx_rots[:,:,n] = [[cos(phid), 0, - np.sin(phid)], 
                              [0, 1, 0],
                              [np.sin(phid), 0, np.cos(phid)]]
            
        
        for traj in trajectories:
            traj += np.vstack((xcs, ycs, zcs))

            for n, pos in enumerate(traj.T):
                tmp_pos = np.dot(pos, xy_rots[:,:,n])
                new_pos = np.dot(tmp_pos, zx_rots[:,:,n])
                traj[:,n] = new_pos

        return trajectories

def reduce_params(paramtree, measuretree):

    ''' This functions changes the paramters so that the average behaviour
    complies with the measures.

    input:
    paramtree : ParamTree instance
    measures : a dictionary of measures
    paramtree is modified in place
    '''
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
        print "The measures dictionary should contain at least the following "
        print "keys "
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

    '''re-creates the Kinetochores Dynamics datas structure
    i.e a dictionnary of "waves" for each element
    returns a Metaphase instance
    '''
    restree = ResultTree(xmlfname)
    param_root = restree.root.find('parameters')
    paramtree = ParamTree(root = param_root)
    params = paramtree.relative_dic
    measure_root = restree.root.find('measures')
    measuretree = ParamTree(root = measure_root, adimentionalized = False)
    metaphase = Metaphase(paramtree = paramtree,
                          measuretree = measuretree)
    
    traj_matrix = restree.get_all_trajs()
    pluged_matrix = restree.get_all_pluged()
    mero_matrix = restree.get_all_mero()
    state_hist_matrix = restree.get_all_plug_state()
    KD = KinetoDynamics(params)
    KD.spbR.traj = traj_matrix[:,0]
    KD.spbL.traj = traj_matrix[:,1]
    N = int(params['N'])
    Mk = int(params['Mk'])
    col_num = 2
    state_num = 0
    for n in range(N):
        KD.chromosomes[n].righttraj = traj_matrix[:, col_num]
        col_num += 1
        KD.chromosomes[n].lefttraj = traj_matrix[:, col_num]
        col_num += 1
        KD.chromosomes[n].pluged_history = (pluged_matrix[:, n*2: n*2 + 2])
        KD.chromosomes[n].mero_history = (mero_matrix[:, n*2 : n*2 + 2])
        for m in range(Mk):
            KD.chromosomes[n].rplugs[m].traj = traj_matrix[:, col_num]
            col_num += 1
            KD.chromosomes[n].rplugs[m].state_hist = state_hist_matrix[:, state_num]
            state_num += 1
        for m in range(Mk):
            KD.chromosomes[n].lplugs[m].traj = traj_matrix[:, col_num]
            col_num += 1
            KD.chromosomes[n].lplugs[m].state_hist = state_hist_matrix[:, state_num]
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
    
    
