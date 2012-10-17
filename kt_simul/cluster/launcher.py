"""
Launcher can run several simulation with same parameters.
It handle data storing.
"""

import logging
import multiprocessing
from multiprocessing import Pool

from kt_simul.core.xml_handler import ParamTree
from kt_simul.core.simul_spindle import Metaphase, PARAMFILE, MEASUREFILE
from kt_simul.core import parameters

from kt_simul.io import SimuIO
from kt_simul.draw import Drawer


class Launcher:
    """
    """

    def __init__(self, results_path,
                 nsimu,
                 ncore=None,
                 paramtree=None, measuretree=None,
                 paramfile=PARAMFILE, measurefile=MEASUREFILE,
                 verbose=True, use_multi_process=True):
        """

        :results_path: The path where simulation results are stored
        :type results_path: string

        :nsimu: Number of simulation to run
        :type nsimu: int

        :ncore: Number of process to launch at the same time
        :type ncore: int

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

        # Reduce parameters
        parameters.reduce_params(self.paramtree, self.measuretree)

        self.results_path = results_path
        self.nsimu = nsimu
        if not ncore:
            self.ncore = multiprocessing.cpu_count() + 1
        else:
            self.ncore = ncore
        self.use_multi_process = use_multi_process

    def run(self):
        """
        """

        logging.info("Starting %i simulations on %i core" % \
            (self.nsimu, self.ncore))

        logging.info("Building parameters set")
        args = []
        i = 0
        for job in range(self.nsimu):
            args.append({"verbose" : self.verbose,
                    "paramtree" : self.paramtree,
                    "measuretree" : self.measuretree,
                    "process_number" : i,
                    "total_process_number" : self.nsimu})
            i += 1

        logging.info("Launching simulations")

        p = Pool(processes=self.ncore)
        p.map(run_one, args)
        p.close()
        p.join()

        logging.info("Simulations are done")

def run_one(*args):
    """
    This function need to be outisde Launcher class to allow
    multiprocessing module to work
    """

    # Retrieving parameters
    for k, v in args[0].iteritems():
        globals()[k]=v

    logging.info("Launching simulation %i" % process_number)

    meta = Metaphase(verbose=verbose,
        paramtree=paramtree,
        measuretree=measuretree)
    meta.simul()

    logging.info("Simulation %i done" % process_number)

