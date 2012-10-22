# -*- coding: utf-8 -*-

"""
Usefull functions to plot with matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt


def dic_plot(plot_data, fname=None):
    """
    Plot curve with dictionnary description
    """

    d = plot_data

    fig, ax = plt.subplots(1)

    # Add label and title
    if 'label' in d['xaxis'].keys():
        ax.set_xlabel(d['xaxis']['label'])
    if 'label' in d['yaxis'].keys():
        ax.set_ylabel(d['yaxis']['label'])
    if 'title' in d.keys():
        ax.set_title(d['title'])

    # Set plot configuration
    plotter = ax.plot
    if 'logx' in d.keys() and d['logx'] == True and\
        'logy' in d.keys() and d['logy'] == True:
        plotter = ax.loglog
    if 'logx' in d.keys() and d['logx'] == True:
        plotter = ax.semilogx
    if 'logy' in d.keys() and d['logy'] == True:
        plotter = ax.semilogy

    xaxis = d['xaxis']['data']
    legends = []
    legends_label = []

    for yaxis in d['yaxis']['axis']:
        mu = yaxis['data']
        color = yaxis['color']

        yplot, = plotter(xaxis, mu, color, label=yaxis['legend'])
        if 'legend' in yaxis.keys():
            legends.append(yplot)
            legends_label.append(yaxis['legend'])

        if 'error' in yaxis.keys():
            sigma = yaxis['error']
            ax.fill_between(xaxis, mu + sigma, mu - sigma,
                facecolor=color, alpha=0.5)

    lgd = ax.legend(tuple(legends), tuple(legends_label),
        loc='upper center', bbox_to_anchor=(0.5,-0.1))
    ax.grid()

    if fname:
        plt.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()
