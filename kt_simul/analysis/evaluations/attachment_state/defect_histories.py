from kt_simul.analysis.evaluations import Evaluation

import numpy as np


class DefectHistories(Evaluation):
    """
    Returns the history of the various attachment states.
    Takes a KinetochoreDynamics instance as unique argument

    output:
    num_defects, were_defects: dectionnaries

    The following outputs are established on a per chromosome
    (i.e. kt pairs) basis

    Those outputs are returned in a dictionnary called num_defects:
    num_defects = {"amphitelic":the number of correctely plugged chromosomes,
               "merotelic":the number of merotelic chromosomes,
               "monotelic":the number of monotelic chromosomes,
               "syntelic":the number of syntelic chromosomes
               "unattached":the number of unpluged chromosomes}
    Normaly num_plugs+num_synt+num_mono+num_mero+num_unplugs = N

    were_plug, were_mero, etc. : 2D np.arrays of booleans giving the state
    of each chromosome
    """

    name = "Defect Histories"
    description = "Returns the history of the various attachment states. Takes a KinetochoreDynamics instance as unique argument"
    group = "attachement_state"
    enable = True

    def __init__(self,):
        pass

    def run(self, KD):
        """
        """

        meros = [np.array(ch.erroneous_history) for ch in KD.chromosomes]
        plugs = [np.array(ch.correct_history) for ch in KD.chromosomes]

        meros = np.hstack(meros)
        plugs = np.hstack(plugs)
        mh_b = meros.astype(bool)
        ph_b = plugs.astype(bool)

        # Put all this in a dictionary
        were_defects = {'amphitelic': self.were_plug(mh_b, ph_b),
                       'merotelic': self.were_mero(mh_b, ph_b),
                       'monotelic': self.were_mono(mh_b, ph_b),
                       'syntelic': self.were_synt(mh_b, ph_b),
                       'unattached': self.were_unat(mh_b, ph_b)}

        num_defects = {}

        for key, val in were_defects.items():
            num_defects[key] = val.sum(axis=1)
        return num_defects, were_defects

    def were_plug(self, mh_b, ph_b):
        """
        Test wether chromosomes were  amphitelic
                    Correct  Merotelic
        -------------------------------
             Left  |  True     False
        AND  Right |  True     False
        OR
             Left  |  False    True
        AND  Right |  False    True
        """

        correct_r = ph_b[:, ::2] & np.logical_not(mh_b[:, ::2])
        correct_l = ph_b[:, 1::2] & np.logical_not(mh_b[:, 1::2])
        and_x = correct_r & correct_l

        mero_r = mh_b[:, ::2] & np.logical_not(ph_b[:, ::2])
        mero_l = mh_b[:, 1::2] & np.logical_not(ph_b[:, 1::2])
        and_x_mero = mero_r & mero_l

        return (and_x | and_x_mero)

    def were_synt(self, mh_b, ph_b):
        """
        #The number of syntelic chromosomes
        #        Correct  Merotelic OR  Correct  Merotelic
        #------------------------------------------------
        # Left  |  True     False        False    True
        # Right | False     True         True     False
        """
        xor_l = mh_b[:, ::2] ^ ph_b[:, ::2]   # left
        xor_r = mh_b[:, 1::2] ^ ph_b[:, 1::2]  # right
        xor_x = ph_b[:, ::2] ^ ph_b[:, 1::2]  # crossed
        xor_t = xor_l & xor_r & xor_x
        return xor_t

    def were_mero(self, mh_b, ph_b):
        """
        The number of merotelic chromosomes
                  Correct  Merotelic
             Left   True     True
         OR  Right  True     True
        """
        and_r = mh_b[:, ::2] & ph_b[:, ::2]
        and_l = mh_b[:, 1::2] & ph_b[:, 1::2]
        or_x = and_r | and_l
        return or_x

    # On a phenotype point of view, one should futher distinguish between merotelic chromosomes :
    # The ones that are amphitelic - merotelic (i.e kts tend to segregate correctly)
    # The ones that are syntelic - merotelic (both kts go to the same poles)
    # The ones that are monotelic - merotelic (one kt is merotelic and the other is unattached)
    # The ones that are cut (one or both kts are equaly attached to both poles)

    def were_amphi_mero(self, mh_b, ph_b):

        and_r = mh_b[:, ::2] & ph_b[:, ::2]
        and_l = mh_b[:, 1::2] & ph_b[:, 1::2]

        nor_r = np.logical_not(mh_b[:, ::2] | ph_b[:, ::2])
        nor_l = np.logical_not(mh_b[:, 1::2] | ph_b[:, 1::2])

        return (and_r & nor_l) ^ (and_l & nor_r)

    def were_mono(self, mh_b, ph_b):
        """
                    Correct  Merotelic
              Left   False     False
         XOR  Right  False     False
        """
        xor_r = mh_b[:, 1::2] ^ ph_b[:, 1::2]
        xor_l = mh_b[:, ::2] ^ ph_b[:, ::2]

        # The number of monotelic chromosomes
        nand_r = np.logical_not(mh_b[:, ::2] | ph_b[:, ::2])
        nand_l = np.logical_not(mh_b[:, 1::2] | ph_b[:, 1::2])
        xor_x = (nand_r & xor_r) ^ (nand_l & xor_l)
        return xor_x

    def were_unat(self, mh_b, ph_b):
        """
        The number of unattached chromosomes
                     Correct  Merotelic
         -------------------------------
              Left  |  False     False
         AND  Right |  False     False
        """
        # xor_r = mh_b[:, 1::2] ^ ph_b[:, 1::2]
        # xor_l = mh_b[:, ::2] ^ ph_b[:, ::2]
        nand_r = np.logical_not(mh_b[:, ::2] | ph_b[:, ::2])
        nand_l = np.logical_not(mh_b[:, 1::2] | ph_b[:, 1::2])
        and_x = nand_r & nand_l
        return and_x
