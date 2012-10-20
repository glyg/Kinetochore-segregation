from kt_simul.analysis.pool_evaluations import PoolEvaluation

import numpy as np


class KinetoAttachment(PoolEvaluation):
    """
    """

    name = "Kineto Attachment"
    description = """Return the naive attachment state of kinetochore
    (correct, erroneous or unattached)"""
    group = "attachment_state"
    enable = True

    def __init__(self,):
        pass

    def run(self, raw_path, eval_results):
        """
        """

        return True
