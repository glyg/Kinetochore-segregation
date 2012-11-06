import logging
import os

import numpy as np

from kt_simul.analysis.pool_evaluations import find_pool_evaluations


class Process:
    """
    """

    def __init__(self, results_path):
        """
        """
        self.results_path = results_path

        # Where simulations raw data are stored
        self.raw_path = os.path.join(self.results_path, 'raw')

        # Where to store analysis files (plot, etc)
        self.eval_results = os.path.join(self.results_path, 'analysis')

        self.observations = {}

        logging.info("Processor initialised")

    def evaluate(self, name=None, groups=[], debug=False, run_all=False, verbose=True):
        """
        """

        logging.info("Starting pool evaluations")
        all_pool_evaluations = find_pool_evaluations(name=name, groups=groups, run_all=run_all)

        if not all_pool_evaluations:
            logging.info("No pool evaluations found")
            return False

        for pool_evaluation in all_pool_evaluations:
            logging.info("Running %s" % pool_evaluation.name)
            if debug:
                result = pool_evaluation().run(self.results_path,
                                                self.raw_path,
                                                self.eval_results,
                                                verbose=verbose)
                logging.info("%s done" % pool_evaluation.name)
            else:
                try:
                    result = pool_evaluation().run(self.results_path,
                                                self.raw_path,
                                                self.eval_results,
                                                verbose=verbose)
                    logging.info("%s done" % pool_evaluation.name)
                except Exception as e:
                    result = np.nan
                    logging.info("%s returns errors : %s" % (pool_evaluation.name, e))

            name = pool_evaluation.name
            self.observations[name] = result

        logging.info("All pool evaluations processed")

        del result
        return self.observations
