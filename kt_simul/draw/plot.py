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

    if "params_box" in d.keys():
        fig = plt.figure(figsize=(20, 6))
        ax = plt.subplot2grid((1, 6), (0, 0), colspan=4)
        box = plt.subplot2grid((1, 6), (0, 4), colspan=2)

        box.set_axis_off()
    else:
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)

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

    if "params_box" in d.keys():
        txtstr = "Simulation parameters\n\n"

        for p in plot_data["params_box"]:
            line = "%s : %s" % (p['name'], p['data'])
            txtstr += line + "\n"

        # these are matplotlib.patch.Patch properies
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        box.text(0, 0.9, txtstr,
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=14,
                 bbox=props)

    # plt.tight_layout()

    if fname:
        plt.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()
