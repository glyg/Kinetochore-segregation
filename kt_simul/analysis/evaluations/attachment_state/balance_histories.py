from kt_simul.analysis.evaluations import Evaluation

import numpy as np


class BalanceHistories(Evaluation):
    """
    input: a kt_simul.spindle_dynamics.KinetochoreDynamics()
    after a simulation instance.

    Returns the difference between the number of correctly and erroneously
    attached plug sites when those two are != 0. This is  stored in an 2D array
    for which each line gives the number of kt in each of the possible cases
    (i.e balance = -Mk+2, -Mk+3, .., Mk-2) for each time point. """

    name = "Balance Histories"
    description = """Returns the difference between the number of correctly and
                    erroneously attached plug sites when those two are != 0."""
    group = "attachement_state"
    enable = True

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        Mk = int(KD.params['Mk'])
        num_steps = len(KD.spbR.traj)

        # The number of cases is (Mk-1) * 2 + 1

        balance = np.vsplit(np.zeros((2 * Mk - 1, num_steps)), 2 * Mk - 1)

        merotelic_types = {'corrected': np.zeros(num_steps),
                           'cut': np.zeros(num_steps),
                           'monotelic': np.zeros(num_steps),
                           'syntelic': np.zeros(num_steps)}

        for ch in KD.chromosomes:

            mh = np.array(ch.erroneous_history)
            ph = np.array(ch.correct_history)

            rbalance = np.zeros(num_steps)
            lbalance = np.zeros(num_steps)
            for j, (m, p) in enumerate(zip(mh, ph)):
                if m[0] * p[0] == 0:  # So the kt is not merotelic
                    rbalance[j] = np.nan
                else:
                    rbalance[j] = p[0] - m[0]
                if m[1] * p[1] == 0:
                    lbalance[j] = np.nan
                else:
                    lbalance[j] = p[1] - m[1]

                if np.isfinite(lbalance[j]) or np.isfinite(rbalance[j]):
                    if p[0] == m[0]:
                        if p[0] != 0:
                            merotelic_types['cut'][j] += 1
                        else:
                            merotelic_types['monotelic'][j] += 1
                    if p[1] == m[1]:
                        if p[1] != 0:
                            merotelic_types['cut'][j] += 1
                        else:
                            merotelic_types['monotelic'][j] += 1

                    # They are pulled the same way
                    if (p[1] - m[1]) * (p[0] - m[0]) < 0:
                        merotelic_types['syntelic'][j] += 1

                    # they are pulled appart
                    if (p[1] - m[1]) * (p[0] - m[0]) > 0:
                        merotelic_types['corrected'][j] += 1

            for i, bal in enumerate(balance):
                bal += np.array([rbalance == i - Mk + 2]).flatten().astype(int)
                bal += np.array([lbalance == i - Mk + 2]).flatten().astype(int)

        return np.vstack(balance), merotelic_types
