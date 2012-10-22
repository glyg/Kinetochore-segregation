from kt_simul.analysis.evaluations import Evaluation

import numpy as np


class ChromoAttachment(Evaluation):
    """
    """

    name = "Chromosome Attachment"
    description = """Return the chromosomes attachment state"""
    group = "attachment_state"
    enable = False

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        num_steps = KD.spbR.traj.size
        N = int(KD.params['N'])
        Mk = int(KD.params['Mk'])

        value_max = N * Mk * 2

        attach_state = {'monotelic': np.zeros(num_steps, dtype='float'),
                        'incorrect_attached': np.zeros(num_steps, dtype='float'),
                        'unattached': np.zeros(num_steps, dtype='float')
                       }

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
