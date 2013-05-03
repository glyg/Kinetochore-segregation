# -*- coding: utf-8 -*-
"""
This module provides the core simulation functionalities.

See Gay et al. J. Cell Biol., 2012 http://dx.doi.org/10.1083/jcb.201107124
The original framework was adapted from:
Civelekoglu-Scholey et al. Biophys. J. 90(11), 2006
http://dx.doi.org/10.1529/biophysj.105.078691

"""

import logging
import numpy as np
import collections

# Local imports
from kt_simul.core.spindle_dynamics import KinetoDynamics
from kt_simul.io.xml_handler import ParamTree
from kt_simul.analysis import evaluations
from kt_simul.core import parameters
from kt_simul.utils.progress import print_progress
from kt_simul.utils.format import pretty_dict

logger = logging.getLogger(__name__)

__all__ = ["Metaphase", "PARAMFILE", "MEASUREFILE"]

CURRENT_DIR = parameters.CURRENT_DIR
ROOT_DIR = parameters.ROOT_DIR
PARAMFILE = parameters.PARAMFILE
MEASUREFILE = parameters.MEASUREFILE
MEASURETREE = parameters.MEASURETREE
MEASURES = parameters.MEASURES


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
        logger = logging.getLogger(__name__)
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

        logger.info('Parameters loaded')

        params = self.paramtree.relative_dic
        # Reset explicitely the unit parameters to their
        # dimentionalized value
        params['Vk'] = self.paramtree.absolute_dic['Vk']
        params['Fk'] = self.paramtree.absolute_dic['Fk']
        params['dt'] = self.paramtree.absolute_dic['dt']

        self.KD = KinetoDynamics(params, initial_plug=initial_plug)
        dt = self.paramtree.absolute_dic['dt']
        duration = self.paramtree.absolute_dic['span']
        self.num_steps = int(duration / dt)
        self.KD.anaphase = False
        self.timelapse = np.arange(0, duration, dt)
        self.report = []
        self.delay = -1
        self.observations = {}

        logger.info('Simulation initialized')
        logger.disabled = False

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
            lines.append('Measures:')
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

        # Check is simulation has already be done
        if self.KD.simulation_done:
            raise SimulationAlreadyDone("""A simulation is already done on this
instance. Please create another Metaphase instance to launch a new simulation.""")

        kappa_c = self.KD.params['kappa_c']

        if self.verbose:
            logger.info('Running simulation')
        bef = 0
        log_anaphase_onset = False

        for time_point in range(1, self.num_steps):

            progress = int((time_point * 100.0) / self.num_steps)

            if self.verbose and progress != bef:
                print_progress(int(progress))
                bef = progress

            # Ablation test
            if ablat == time_point:
                if self.verbose:
                    logger.info("Performing ablation")
                self._ablation(time_point, pos=ablat_pos)

            # Anaphase transition ?
            if self._anaphase_test(time_point):
                if not log_anaphase_onset:
                    print_progress(-1)
                    if self.verbose:
                        logger.info("Anaphase onset at %i / %i" %
                                        (time_point, self.num_steps))
                    log_anaphase_onset = True

            self.KD.one_step(time_point)
            # if time_point % 100 == 0:
                # print self.KD.At_mat

        if self.verbose:
            print_progress(-1)

        if self.verbose:
            logger.info('Simulation done')
        self.KD.params['kappa_c'] = kappa_c
        delay_str = "delay = %2d seconds" % self.delay
        self.report.append(delay_str)

        for ch in self.KD.chromosomes:
            ch.calc_correct_history()
            ch.calc_erroneous_history()
            ch.cen_A.calc_toa()
            ch.cen_B.calc_toa()

    def evaluate(self, name=None, groups=[],
                 debug=False,
                 verbose=False,
                 draw=False,
                 run_all=False):
        """
        Passes all the evaluations in kt_simul.analysis.valuations module
        results are stored in the self.observations dictionnary

        TODO: most of the code of this method should be moved to
        kt_simul.analysis.evaluations.__init__
        """
        if not self.KD.simulation_done:
            logger.info("No simulation was runned")
            return False

        if not name and verbose:
            logger.info("Starting evaluations")
        all_evaluations = evaluations.find_evaluations(name=name, groups=groups, run_all=run_all)

        if not all_evaluations:
            if verbose:
                logger.info("No evaluations found")
            return False

        for evaluation in all_evaluations:
            if verbose:
                logger.info("Running %s" % evaluation.name)
            if debug:
                result = evaluation().run(self.KD, draw)
                if verbose:
                    logger.info("%s done" % evaluation.name)
            else:
                try:
                    result = evaluation().run(self.KD, draw)
                    if verbose:
                        logger.info("%s done" % evaluation.name)
                except Exception as e:
                    result = np.nan
                    if verbose:
                        logger.info("%s returns errors : %s" % (evaluation.name, e))

            if name and not run_all:
                return result
            else:
                current_name = evaluation.name.replace(" ", "_")
                self.observations[current_name] = result

        if verbose:
            logger.info("All evaluations processed")

        return self.observations

    def get_report(self, time = 0):
        """
        Print simulation state about a specific time point
        """
        params = self.paramtree.relative_dic

        report = collections.OrderedDict()

        report["Total time (span)"] = params["span"]
        report["Time precision (dt)"] = params["dt"]
        report["Chromosomes (N)"] = params["N"]
        report["Kinetochores (Mk)"] = params["Mk"]
        report["separate"] = ""

        report["Current time"] = "%i\n" % ((time + 1) * float(params["dt"]))

        report["separate"] = ""
        report["spbR pos"] = round(self.KD.spbR.traj[time], 3)
        report["spbL pos"] = round(self.KD.spbL.traj[time], 3)
        report["separate"] = ""

        for ch in self.KD.chromosomes:
            chdict = collections.OrderedDict()
            chdict["correct"] = ch.correct_history[time]
            chdict["erroneous"] = ch.erroneous_history[time]
            for cent in [ch.cen_A, ch.cen_B]:
                cen_dict = collections.OrderedDict()
                cen_dict["position"] = round(cent.traj[time], 3)
                for site in cent.plugsites:
                    site_dict = collections.OrderedDict()
                    site_dict['position'] = round(site.traj[time], 3)
                    site_dict['Plug state'] = site.state_hist[time]

                    cen_dict["PlugSite %i" % site.site_id] = site_dict

                chdict["Centromere %s" % cent.tag] = cen_dict

            report["Chromosome %i" % ch.ch_id] = chdict
            report["separate"] = ""

        return pretty_dict(report)

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
            logger.warning('Missed shot, same player play again!')
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

class SimulationAlreadyDone(Exception):
    pass
