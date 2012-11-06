import numpy as np
import os
import logging

from kt_simul.analysis.explo_pool_evaluations import ExploPoolEvaluation
from kt_simul.draw.plot import dic_plot
from kt_simul.cluster.process import Process


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
        nsimu = int(self.get_nsimu(simu_path))

        for simudir in pool_folder:
            p = Process(results_path=simudir)

            print p.evaluate(name="Mitotic Plate", verbose=False)
            print "###################################"

        # kt_plate = {'dispersion': np.zeros((nsimu, num_steps))
        #            }

        # logging.info("Loading data from simulations files")
        # for i, (simu_id, meta) in enumerate(self.iter_simulations(raw_path,
        #                                                 nsimu=nsimu,
        #                                                 print_progress=True)):
        #     results = meta.evaluate(name="Mitotic Plate", verbose=False)

        #     kt_plate['dispersion'][i] = results['dispersion']

        # logging.getLogger().disabled = False

        # kt_plate['dispersion_std'] = kt_plate['dispersion'].std(axis=0)
        # kt_plate['dispersion'] = kt_plate['dispersion'].mean(axis=0)

        # # Configure and plot the graph
        # logging.info("Plotting results")
        # timelapse = meta_infos.timelapse

        # plot_data = {}
        # plot_data['title'] = "Mitotic plate formation"
        # plot_data['xaxis'] = {'data': timelapse, 'label': 'Time (second)'}
        # plot_data['yaxis'] = {'label': 'Dispersion measure (relative to the spindle length)', 'axis': []}
        # plot_data['error'] = True
        # plot_data['legend'] = False
        # plot_data['limit_y_min'] = 0

        # # Draw parameters box
        # plot_data["params_box"] = [{'name': "Name", 'data': name},
        #                            {'name': "Simulations number", 'data': nsimu},
        #                            {'name': "Lenght Dependance factor", 'data': params["ld_slope"]},
        #                            {'name': 'Lenght Dependance base', 'data': params['ld0']}
        #                      ]

        # # Add annotation about anaphase onset
        # plot_data["annotations"] = []
        # plot_data["annotations"].append({'s': 'Anaphase onset: %i' % ana_onset,
        #                                  'xy': (ana_onset, 0),
        #                                  'xytext': (0,50),
        #                                  'textcoords': 'offset points',
        #                                  'ha': 'center',
        #                                  'va': 'bottom',
        #                                  'bbox': dict(boxstyle='round,pad=0.2',
        #                                               fc='yellow',
        #                                               alpha=0.3),
        #                                  'arrowprops': dict(arrowstyle='->',
        #                                                     # connectionstyle='arc3,rad=0.5',
        #                                                     color='black')
        #                                  })

        # plot_data['yaxis']['axis'].append({'data': kt_plate['dispersion'],
        #                                     'color': 'blue',
        #                                     'error': kt_plate['dispersion_std']
        #                                     })

        # dic_plot(plot_data, os.path.join(eval_results_path, "Mitotic_Plate_Formation.svg"))

        # return kt_plate
