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

        num_steps = KD.spbR.traj.size
        N = int(KD.params['N'])
        Mk = int(KD.params['Mk'])

        value_max = N * Mk * 2

        attach_state = {'correct_attached': np.zeros(num_steps, dtype='float'),
                        'incorrect_attached': np.zeros(num_steps, dtype='float'),
                        'unattached': np.zeros(num_steps, dtype='float')
                       }

        for ch in KD.chromosomes:

            ch_incorrect = np.array(ch.erroneous_history).sum(1)
            ch_correct = np.array(ch.correct_history).sum(1)

            attach_state['correct_attached'] += ch_correct
            attach_state['incorrect_attached'] += ch_incorrect

        attach_state['unattached'] = np.ones(num_steps, dtype='float') * value_max
        attach_state['unattached'] -= attach_state['correct_attached']
        attach_state['unattached'] -= attach_state['incorrect_attached']

        # Scale values on 1
        attach_state['correct_attached'] /= value_max
        attach_state['incorrect_attached'] /= value_max
        attach_state['unattached'] /= value_max

        return attach_state
