# -*- coding: utf-8 -*-

from kt_simul.draw.plot import dic_plot


def draw_kineto_attachment(kt_attach, meta):
    """
    """

    timelapse = meta.timelapse

    # Configure and plot the graph

    plot_data = {}
    plot_data['title'] = "Kinetochores attachment rate"
    plot_data['xaxis'] = {'data': timelapse, 'label': 'Time'}
    plot_data['yaxis'] = {'label': 'Attachment rate', 'axis': []}

    correct_attached_axis = {'data': kt_attach['correct_attached'],
                             'color': 'green',
                             'legend': 'Correct attached kinetochore'}

    incorrect_attached_axis = {'data': kt_attach['incorrect_attached'],
                               'color': 'red',
                               'legend': 'Incorrect attached kinetochore'}

    unattached_axis = {'data': kt_attach['unattached'],
                       'color': 'blue',
                       'legend': 'Unattached kinetochore'}

    plot_data['yaxis']['axis'].append(correct_attached_axis)
    plot_data['yaxis']['axis'].append(incorrect_attached_axis)
    plot_data['yaxis']['axis'].append(unattached_axis)

    dic_plot(plot_data)
