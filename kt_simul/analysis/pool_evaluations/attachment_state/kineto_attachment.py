import numpy as np
import os
import logging

from kt_simul.analysis.pool_evaluations import PoolEvaluation
from kt_simul.draw.plot import dic_plot


class KinetoAttachment(PoolEvaluation):
    """
    """

    name = "Kineto Attachment"
    description = """Process and plot mean and std of the returned values of
    Evaluation called 'Kineto Attachment'"""
    group = "attachment_state"
    enable = True

    def __init__(self,):
        pass

    def run(self, simu_path, raw_path, eval_results_path):
        """
        """

        meta_infos = self.get_simus_params(simu_path)
        num_steps = meta_infos.num_steps
        nsimu = int(self.get_nsimu(simu_path))

        attach_state = {'correct_attached': np.zeros((nsimu, num_steps)),
                        'incorrect_attached': np.zeros((nsimu, num_steps)),
                        'unattached': np.zeros((nsimu, num_steps))
                       }

        logging.info("Loading data from simulations files")
        for i, (simu_id, meta) in enumerate(self.iter_simulations(raw_path,
                                                        nsimu=nsimu,
                                                        print_progress=True)):
            results = meta.evaluate(name="Kineto Attachment", verbose=False)

            attach_state['correct_attached'][i] = results['correct_attached']
            attach_state['incorrect_attached'][i] = results['incorrect_attached']
            attach_state['unattached'][i] = results['unattached']

        logging.getLogger().disabled = False
        # Mean data
        logging.info("Processing data")
        correct_attached = attach_state['correct_attached'].mean(axis=0)
        incorrect_attached = attach_state['incorrect_attached'].mean(axis=0)
        unattached = attach_state['unattached'].mean(axis=0)

        # Configure and plot the graph
        logging.info("Plotting results")
        timelapse = meta_infos.timelapse

        plot_data = {}
        plot_data['title'] = "Kinetochores attachment rate on %i simulations" % nsimu
        plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
        plot_data['yaxis'] = {'label': 'Attachment rate', 'axis': []}

        correct_attached_axis = {'data': correct_attached,
                                 'color': 'green',
                                 'legend': 'Correct attached kinetochore'}

        incorrect_attached_axis = {'data': incorrect_attached,
                                   'color': 'red',
                                   'legend': 'Incorrect attached kinetochore'}

        unattached_axis = {'data': unattached,
                           'color': 'blue',
                           'legend': 'Unattached kinetochore'}

        plot_data['yaxis']['axis'].append(correct_attached_axis)
        plot_data['yaxis']['axis'].append(incorrect_attached_axis)
        plot_data['yaxis']['axis'].append(unattached_axis)

        dic_plot(plot_data, os.path.join(simu_path, "Kineto_Attachment.svg"))

        return True
