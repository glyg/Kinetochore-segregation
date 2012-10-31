# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from kt_simul.analysis.evaluations import Evaluation
from kt_simul.draw.plot import dic_plot


class MitoticPlate(Evaluation):
    """
    """

    name = "Mitotic Plate"
    description = "Retrieve mitotic plate (Kt dispersion) over time"
    group = "plate"
    enable = True

    def __init__(self,):
        pass

    def run(self, KD, draw=False):
        """
        """

        N = int(KD.params["N"])
        ana_onset = int(KD.params["t_A"])

        timelapse = np.arange(0, KD.duration, KD.dt)
        num_steps = KD.duration / KD.dt

        spindle_length = abs(KD.spbR.traj - KD.spbL.traj)

        max_kt = np.zeros((N * 2, num_steps), dtype="float")
        min_kt = np.zeros((N * 2, num_steps), dtype="float")

        for i, ch in enumerate(KD.chromosomes):
            id1 = (i * 2) -1
            id2 = (i * 2)

            max_kt[id1] = ch.cen_A.traj
            max_kt[id2] = ch.cen_B.traj

            min_kt[id1] = ch.cen_A.traj
            min_kt[id2] = ch.cen_B.traj

        max_kt_distance = max_kt.max(axis=0)
        min_kt_distance = min_kt.min(axis=0)

        dispersion = (abs(max_kt_distance - min_kt_distance) / spindle_length)

        # if not draw:
        #     return trajectories

        # Draw attachment state with matplotlib
        timelapse = np.arange(0, KD.duration, KD.dt)

        plot_data = {}
        plot_data['title'] = "Mitotic plate formation"
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time (second)'}
        plot_data['yaxis'] = {'label': "Dispersion measure (relative to the spindle length)", 'axis': []}
        plot_data['logx'] = False
        plot_data['legend'] = False
        plot_data['limit_y_min'] = 0

        # Add annotation about anaphase onset
        plot_data["annotations"] = []
        plot_data["annotations"].append({'s': 'Anaphase onset: %i' % ana_onset,
                                       'xy': (ana_onset, 0),
                                       'xytext': (0, 20),
                                       'textcoords': 'offset points',
                                       'ha': 'center',
                                       'va': 'bottom',
                                       'bbox': dict(boxstyle='round,pad=0.2',
                                                    fc='yellow',
                                                    alpha=0.3),
                                       'arrowprops': dict(arrowstyle='->',
                                                          # connectionstyle='arc3,rad=0.5',
                                                          color='black')
                                       })

        plot_data['yaxis']['axis'].append({'data': dispersion,
                                          'color': 'blue',
                                          })

        dic_plot(plot_data, fname="../noldep+.svg")

        # return trajectories
