import numpy as np
import os
import logging

import matplotlib.pyplot as plt

from kt_simul.analysis.explo_pool_evaluations import ExploPoolEvaluation
from kt_simul.draw.plot import dic_plot
from kt_simul.cluster.process import Process
from kt_simul.utils.color import color


class MitoticPlate(ExploPoolEvaluation):
    """
    """

    name = "Mitotic Plate"
    description = "Retrieve mitotic plate (Kt dispersion) over time"
    group = "plate"
    enable = True

    def __init__(self,):
        pass

    def run(self, simu_path, pool_folder, eval_results_path):
        """
        """

        nparams = len(pool_folder)
        plog = self.get_param_from_log(simu_path)
        nsimu = int(plog['nsimu'])
        param_name = plog["parameter to explore name"]

        kt_plate = []

        logging.info("Running processes on each pool")
        for i, simudir in enumerate(sorted(pool_folder)):

            logging.info(color("Running processes on pool called: %s" % os.path.split(simudir)[-1], 'BOLD'))

            p = Process(results_path=simudir)
            resu = p.evaluate(name="Mitotic Plate", draw=False, verbose=False)

            kt_plate.append(resu)

        # Configure and plot the graph
        logging.info("Plotting results")
        timelapse = np.arange(0, plog['duration'], plog['dt'])

        plot_data = {}
        plot_data['title'] = "Mitotic plate formation with %s variable" % param_name
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time (second)'}
        plot_data['yaxis'] = {'label': 'Dispersion measure (relative to the spindle length)', 'axis': []}
        plot_data['error'] = False
        plot_data['legend'] = True
        # plot_data['limit_y_min'] = 0

        # Draw parameters box
        plot_data["params_box"] = [{'name': "Name", 'data': plog["name"]},
                                   {'name': "Simulations number", 'data': nsimu},

                             ]

        # User matplotlib to get color gradient
        cm = plt.get_cmap('gist_rainbow')

        for i, resu in enumerate(kt_plate):

            plot_color = cm(1. * i / nparams)

            plot_data['yaxis']['axis'].append({'data': resu['dispersion'],
                                               'color': plot_color,
                                               'legend': "%s" % resu['params'][param_name]
                                    })


        dic_plot(plot_data, os.path.join(eval_results_path, "Mitotic_Plate_Formation.svg"))

        return kt_plate
