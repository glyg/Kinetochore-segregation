from kt_simul.analysis.evaluations import Evaluation

import numpy as np


class KinetoAttachment(Evaluation):
    """
    """

    name = "Kineto Attachment"
    description = """Return the naive attachment state of kinetochore
    (correct, erroneous or unattached)"""
    group = "attachment_state"
    enable = True

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        Mk = int(KD.params['Mk'])
        num_steps = len(KD.spbR.traj)

        balance = np.vsplit(np.zeros((2 * Mk - 1, num_steps)), 2 * Mk - 1)

        merotelic_types = {'corrected': np.zeros(num_steps),
                           'cut': np.zeros(num_steps),
                           'monotelic': np.zeros(num_steps),
                           'syntelic': np.zeros(num_steps)}

        for ch in KD.chromosomes:

            mh = np.array(ch.erroneous_history)
            ph = np.array(ch.correct_history)

        return mh
