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

    plt.figure()

    # Add label and title
    if 'label' in d['xaxis'].keys():
        plt.xlabel(d['xaxis']['label'])
    if 'label' in d['yaxis'].keys():
        plt.ylabel(d['yaxis']['label'])
    if 'title' in d.keys():
        plt.title(d['title'])

    xaxis = d['xaxis']['data']
    legends = []

    for yaxis in d['yaxis']['axis']:
        #plt.semilogx(xaxis, yaxis['data'], yaxis['color'])
        plt.plot(xaxis, yaxis['data'], yaxis['color'])
        if 'legend' in yaxis.keys():
            legends.append(yaxis['legend'])

    plt.legend(tuple(legends), 'best')

    if fname:
        plt.savefig(fname)
    else:
        plt.show()
