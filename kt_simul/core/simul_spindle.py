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
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal
from xml.etree.ElementTree import Element, SubElement, tostring
from Image import fromarray as image_fromarray

import pyximport
pyximport.install()

# Local imports
from kt_simul.core.spindle_dynamics import KinetoDynamics
from kt_simul.core.xml_handler import ParamTree, indent, ResultTree
from kt_simul.analysis.evaluate import evaluations
import parameters
import utils


__all__ = ["Metaphase", "PARAMFILE",
           "MEASUREFILE", "get_fromfile"]

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
PARAMFILE = os.path.join(ROOT_DIR, 'default', 'params.xml')
MEASUREFILE = os.path.join(ROOT_DIR, 'default', 'measures.xml')
MEASURETREE = ParamTree(MEASUREFILE, adimentionalized = False)
MEASURES = MEASURETREE.absolute_dic


class Metaphase(object):
    """
    An instance of the Metaphase class is a wrapper around
    the whole simulation.

    **Typical usage**

    >>> from kt_simul.simul_spindle import Metaphase
    >>> m = Metaphase()
    >>> m.simul()
    >>> m.show_trajs()
    >>> m.write_results('examples/docstring_results.xml',
                        'examples/docstring_data.npy')

    **From an already runned simulation**

    >>> from kt_simul.simul_spindle import Metaphase
    >>> m1 = get_fromfile('examples/docstring_results.xml')
    >>> m1.show_one(1) #This shows the trajactory of the chromosome 1
    >>> m2 = Metaphase(m1.paramtree, m1.measuretree) #A new simulation
    >>> m2.simul(ablat = 600) #this time with spindle ablation

    """


    def __init__(self,  paramtree=None, measuretree=None,
                 paramfile=PARAMFILE, measurefile=MEASUREFILE,
                 initial_plug='random', reduce_p=True,
                 verbose=False):

        """
        Metaphase instanciation method


        :param duration: The duration of the mitosis in seconds (defaults to 900)
        :type duration: float

        :param paramtree: The paramtree contains the parameters for the simulation
            if paramtree is None, the parameters are read
            from the file paramfile. Defaults to None.
        :type paramtree: ParamTree instance or None

        :param measuretree: The measuretree contains the observed characteristics
            of the mitosis e.g. metaphase spindle elongation rate, etc.
            if measuretree is None, the measures are read from the file
            indicated by the measurefile argument. Defaults to None.
        :type measuretree: ParamTree instance or None

        :param paramfile: Path to a xml file to read the parameters from. Defaults to the
            file params.xml in the module's default directory. Other parameter
            files can be produced by editing and changing the default one.
            If the paramtree argument is not None,  paramfile is ignored
        :type paramfile: string

        :param measurefile: Path to a xml file to read the measures from. Defaults to the
            file measures.xml in the module's default directory.
            Other measure files can be produced by editing and changing
            the default one. If the measuretree argument is not None, measurefile
            is ignored
        :type measurefile: string

        :param initial_plug: Defines globally the initial attachment states.
            This argument can have the following values:
                * 'null': all kinetochores are detached
                * 'amphitelic': all chromosmes are amphitelic
                * 'random': all attachement site can be bound to
                        either pole or deteched with equal prob.
                * 'monotelic': right kinetochores are attached to the same pole,
                           left ones are detached
                * 'syntelic' : all kinetochores are attached to the same pole
        :type initial_plug: string or None

        :param reduce_p: If True, changes the parameters according to the measures
            so that the simulation average behaviour complies with
            the data in the measures dictionary
        :type reduce_p: bool

        """

        # Enable or disable log console
        self.verbose = verbose
        logger = logging.getLogger()
        if not self.verbose:
            logger.disabled = True
        else:
            logger.disabled = False

        if paramtree is None:
            self.paramtree = ParamTree(paramfile)
        else:
            self.paramtree = paramtree
        if measuretree is None:
            self.measuretree = ParamTree(measurefile, adimentionalized=False)
        else:
            self.measuretree = measuretree
        if reduce_p:
            parameters.reduce_params(self.paramtree, self.measuretree)

        logging.info('Parameters loaded')

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

        logging.info('Simulation initialized')

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
        """
        The simulation main loop.

        :param ablat: Timepoint at which ablation takes place. If None (default)
            no ablation is performed.
        :type ablat: float, optional

        """

        #dt = self.KD.params['dt']
        kappa_c = self.KD.params['kappa_c']

        logging.info('Running simulation')
        bef = 0
        log_anaphase_onset = False

        for time_point in range(1, self.num_steps):

            progress = int((time_point * 100.0) / self.num_steps)

            if self.verbose and progress != bef:
                utils.progress(int(progress))
                bef = progress

            # Ablation test
            if ablat == time_point:
                logging.info("Performing ablation")
                self._ablation(time_point, pos=ablat_pos)

            # Anaphase transition ?
            if self._anaphase_test(time_point):
                if not log_anaphase_onset:
                    utils.progress(-1)
                    logging.info("Anaphase onset at %i / %i" %
                        (time_point, self.num_steps))
                    log_anaphase_onset = True

            self.KD.one_step(time_point)

        if self.verbose:
            utils.progress(-1)

        logging.info('Simulation done')
        self.KD.params['kappa_c'] = kappa_c
        delay_str = "delay = %2d seconds" % self.delay
        self.report.append(delay_str)

        for ch in self.KD.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()
            ch.cen_A.calc_toa()
            ch.cen_B.calc_toa()

    def evaluate(self):
        """
        Passes all the evaluations in eval_simul.py
        results are stored in the self.observations dictionnary
        """
        if not self.KD.simulation_done:
            logging.info("No simulation was runned")
            return False

        for name, function in evaluations().iteritems():
            self.observations[name] = function(self.KD)

        return True

    def _anaphase_test(self, time_point):
        """
        Returns True if anaphase has been executed.
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
        """
        Simulates a laser ablation: detaches all the kinetochores
        and sets the midzone stall force Fmz to 0

        :param pos: Position of the laser beam within the spindle
        :type pos: float or None, optional

        """
        if pos == None: pos = self.KD.spbR.pos
        if not self.KD.spbL.pos <= pos <= self.KD.spbR.pos:
            logging.warning('Missed shot, same player play again!')
            return
        self.KD.params['Fmz'] = 0.
        self.KD.params['k_a'] = 0.
        self.KD.params['k_d0'] = 0.
        self.KD.A0_mat = self.KD.time_invariantA()

        for plugsite in self.KD.spindle.all_plugsites:
            if pos < plugsite.pos and plugsite.plug_state == - 1:
                plugsite.set_plug_state(0, time_point)
            elif pos > plugsite.pos and plugsite.plug_state == 1:
                plugsite.set_plug_state(0, time_point)

    def _plug_checkpoint(self):
        """
        If the spindle assembly checkpoint is active, returns True
        if all chromosomes are plugged by at least one kMT, False
        otherwise.
        """
        sac = self.KD.params['sac']
        if sac == 0:
            return True
        for ch in self.KD.chromosomes :
            if not ch.cen_A.is_attached() or not ch.cen_B.is_attached():
                return False
        return True

    def _mero_checkpoint(self):
        """
        :return: The total number of merotellic kT
        """
        nb_mero = 0
        for ch in self.KD.chromosomes :
            if np.any(ch.erroneous()) :
                nb_mero += sum(ch.erroneous())
                #print "active checkpoint"
        return nb_mero

    def _mplate_checkpoint(self):
        """
        :return: Returns True if each kinetochore is in the proper half
        of the spindle
        """
        for ch in self.KD.chromosomes:
            ktR = ch.cen_A.pos
            ktL = ch.cen_B.pos
            if min(ktR, ktL) <= 0 and max(ktR, ktL) >= 0:
                return True
        return True
