# -*- coding: utf-8 -*-
"""
This module is used to lighten Metaphase class
"""

import time
import logging

import numpy as np
import matplotlib.pyplot as plt

from xml.etree.ElementTree import Element, SubElement, tostring

from kt_simul.core.xml_handler import ParamTree, indent, ResultTree
from kt_simul.analysis.eval_simul import evaluations

class Processor():
    """
    Processor class can performs various things with simulation results
    such as writing results in files, plotting graph, making movie (not yet),
    etc
    """

    def __init__(self, meta_instance):
        """
        Parameters
        ----------

        meta_instance : Metaphase instance
            Should have been already perform a simulation
        """

        self.meta = meta_instance
        self.KD = self.meta.KD
        self.timelapse = self.meta.timelapse
        self.num_steps = self.meta.num_steps
        self.paramtree = self.meta.paramtree
        self.measuretree = self.meta.measuretree
        self.observations = {}

    def evaluate(self):
        """
        Passes all the evaluations in eval_simul.py
        results are stored in the self.observations dictionnary
        """
        if not self.KD.simulation_done:
            logging.info("No simulation was runned")
            return False

        for name, function in evaluations().iteritems():
            self.observations[name] = function(self.KD)

        return True

    def write_results(self, xmlfname = "results.xml",
                        datafname = "data.npy"):
        """ Saves the results of the simulation in two files
        with the parameters, measures and observations in one file
        and the trajectories in the other.

        Keyword arguments:
        ------------------
        xmlfname : string, optional
            name of the xml file where parameters and observations
            will be written
        datafname : string, optional
            name of the file where the trajectories will be written
            file type is determined by the file suffix:
                 - *.npy : data are stored in numpy's binary format
                           (less portable but quite efficient)
                 - *.txt : simple text
                 - *.txt.gz : text files compressed transparently

        Any other suffix will be saved as plain text. Column index for
        each trajectory is an attribute of the corresponding element
        in the xml file.

        TODO : This function should be chopped off some how, it's messy

        """
        if not hasattr(self, 'observations'):
            self.evaluate()

        chromosomes = self.KD.chromosomes
        wavelist = []
        out = file(xmlfname, 'w+')
        out.write('<?xml version="1.0"?>\n')
        today = time.asctime()
        experiment = Element("experiment", date=today, datafile=datafname)
        experiment.append(self.paramtree.root)
        experiment.append(self.measuretree.root)

        #right SPB
        spbR = SubElement(experiment, "trajectory", name = "rightspb",
                          column='0', units='mu m')
        SubElement(spbR, "description").text="right spb trajectory"
        spbRtraj = np.array(self.KD.spbR.traj)
        wavelist.append(spbRtraj)

        #left SPB
        spbL = SubElement(experiment, "trajectory", name = "leftspb",
                          column='1', units='mu m')
        SubElement(spbL, "description").text="left spb trajectory"
        spbLtraj = np.array(self.KD.spbL.traj)
        wavelist.append(spbLtraj)

        col_num = 2
        #chromosomes
        for n, ch in enumerate(chromosomes):
            rch = SubElement(experiment, "trajectory", name="centromereA",
                             index = str(n), column=str(col_num), units='mu m')
            text = "chromosome %i centromere A trajectory" % n
            SubElement(rch, "description").text = text
            wavelist.append(ch.cen_A.traj)
            col_num += 1

            SubElement(experiment, "numbercorrect", name="centromereA",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.correct_history[:, 0])
            col_num += 1
            SubElement(experiment, "numbererroneous",
                       name="centromereA", index = str(n),
                       column=str(col_num))
            wavelist.append(ch.erroneous_history[:, 0])
            col_num += 1

            lch = SubElement(experiment, "trajectory", index=str(n),
                             column=str(col_num), units='mu m')
            text = "chromosome %s left kinetochore trajectory" % n
            SubElement(lch, "description").text = text
            wavelist.append(np.array(ch.cen_B.traj))
            col_num += 1
            SubElement(experiment, "numbercorrect", name="centromereB",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.correct_history[:, 1])
            col_num += 1

            SubElement(experiment, "numbererroneous", name="centromereB",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.erroneous_history[:, 1])
            col_num += 1

            #Plug Sites
            for m, plugsite in enumerate(ch.cen_A.plugsites):
                SubElement(experiment, "trajectory", name="plugsite",
                           index=str((n, m)), cen_tag='A',
                           column=str(col_num), units='mu m')
                wavelist.append(plugsite.traj)
                col_num += 1

                SubElement(experiment, "state", name="plugsite",
                           index=str((n, m)), cen_tag='A',
                           column=str(col_num), units='')
                wavelist.append(plugsite.state_hist)
                col_num += 1

            for m, plugsite in enumerate(ch.cen_B.plugsites):

                SubElement(experiment, "trajectory", name="plugsite",
                           index = str((n, m)), cen_tag='B',
                           column=str(col_num), units='mu m')
                wavelist.append(plugsite.traj)
                col_num += 1

                SubElement(experiment, "state", name="plugsite",
                           index = str((n, m)), cen_tag='B',
                           column=str(col_num), units='')
                wavelist.append(plugsite.state_hist)
                col_num += 1

        #Observations
        obs_elem = SubElement(experiment, "observations")
        if not hasattr(self, 'observations'):
            self.evaluate()
        for key, val in self.observations.items():
            SubElement(obs_elem, key).text = str(val)

        #Now we write down the whole experiment XML element
        indent(experiment)
        out.write(tostring(experiment))
        out.close()
        #And the numbers in the file datafname
        dataout = file(datafname, 'w+')
        data = np.vstack(wavelist).T
        if datafname.endswith('.npy'):
            np.save(dataout, data)
        else:
            dataout.write("# desciptor: "+xmlfname+"\n")
            np.savetxt(dataout, data, delimiter=' ')
        logging.info("Simulation saved to file %s and %s "
            % (xmlfname, datafname))

    def write_asoctave(self, xmlfname = "results.xml",
                       datafname = "oct_data.txt"):
        """
        Save results in a file compatible with the tools developped
        to caracterize real life movies under octave.

        Parameters
        ----------
        xmlfname : string
            name of the xml file where the paramters and observations
            are stored, defaults to ``results.xml``
        datafname : string
            name of the xml file where the trajectories coordinates
            are stored, defaults to ``oct_data.txt``

        Note that this will also call self.write_results(xmlfname).

        """
        self.write_results(xmlfname = xmlfname)
        spbRx = self.KD.spbR.traj
        ys = np.zeros(spbRx.shape)
        ts = np.arange(spbRx.shape[0])
        idxs = np.ones(spbRx.shape)
        out_array = np.vstack((spbRx, ys, ts, idxs))

        spbLx = self.KD.spbL.traj
        idxs += 1
        vstack = np.vstack((spbLx, ys, ts, idxs))
        out_array = np.append(out_array, vstack, axis=0)

        for ch in self.KD.chromosomes:
            idxs += 1
            chRx = ch.cen_A.traj
            vstack = np.vstack((chRx, ys, ts, idxs))
            out_array = np.append(out_array, vstack, axis=0)

            idxs += 1
            chLx = ch.cen_B.traj
            vstack = np.vstack((chLx, ys, ts, idxs))
            out_array = np.append(out_array, vstack, axis = 0)

        dataout = file(datafname, 'w+')
        dataout.write("# desciptor: "+xmlfname+"\n")
        np.savetxt(dataout, out_array, delimiter=' ')

    def show_trajs(self, axes = None):
        """
        Plot the different trajectories
        """
        N = int(self.KD.params['N'])
        if axes == None:
            fig = plt.figure()
            axes = fig.gca()
        spbRtraj = self.KD.spbR.traj
        spbLtraj = self.KD.spbL.traj
        axes.plot(self.timelapse, spbRtraj, color='r', ls='-', lw=1)
        axes.plot(self.timelapse, spbLtraj, color='r', ls='-', lw=1)

        fmt_list = ['g-', 'b-', 'm-']
        for n in range(N):
            ch = self.KD.chromosomes[n]
            right_traj = ch.cen_A.traj
            left_traj = ch.cen_B.traj
            fmt = fmt_list[np.mod(n, 3)]
            line1 = axes.plot(self.timelapse, right_traj, fmt, alpha=0.5)
            line2 = axes.plot(self.timelapse, left_traj, fmt, alpha=0.5)
            if n == 0:
                line1[0].set_alpha(1.)
                line2[0].set_alpha(1.)
        axes.set_xlabel('Time (seconds)', fontsize = 'small')
        axes.set_ylabel(u'Distance from center (um)', fontsize = 'small')
        plt.show()

    def show_one(self, n = 0, fig = None):
        """
        Shows chromosome n trajectory and plug state
        """
        dt = self.KD.params['dt']
        ch = self.KD.chromosomes[n]

        if fig == None:
            fig = plt.figure()
        fig.clear()

        #fig.add_subplot(312)
        gridspec = plt.GridSpec(5,1)
        subplotspec = gridspec.new_subplotspec((1,0), rowspan=3)
        traj_ax = fig.add_subplot(subplotspec)
        traj_ax.plot(self.timelapse, ch.cen_A.traj, 'g', lw=2, alpha=0.5)
        traj_ax.plot(self.timelapse, ch.cen_B.traj, 'purple', lw=2, alpha=0.5)
        traj_ax.plot(self.timelapse, self.KD.spbR.traj, 'k')
        traj_ax.plot(self.timelapse, self.KD.spbL.traj, 'k')
        for plugsite in ch.cen_A.plugsites:
            traj_ax.plot(self.timelapse, plugsite.traj, 'g')
        for plugsite in ch.cen_B.plugsites:
            traj_ax.plot(self.timelapse, plugsite.traj, 'purple')
        traj_ax.set_xticks([], '')

        erroneous_hist = ch.erroneous_history
        correct_hist = ch.correct_history

        subplotspec = gridspec.new_subplotspec((0,0), rowspan=1)
        ax = fig.add_subplot(subplotspec, sharex=traj_ax)
        ax.plot(self.timelapse, erroneous_hist[:, 0], 'r',
                label='number of erroneoustellic MTs')
        ax.plot(self.timelapse, correct_hist[:, 0], 'g',
                label='number of correct MTs')
        ax.axis((0, self.num_steps*dt, -0.5, 4.5))
        ax.set_xticks([], '')

        subplotspec = gridspec.new_subplotspec((4,0), rowspan=1)
        ax = fig.add_subplot(subplotspec, sharex=traj_ax)
        ax.plot(self.timelapse, erroneous_hist[:, 1], 'r',
                label='number of erroneous MTs')
        ax.plot(self.timelapse, correct_hist[:, 1], 'purple',
                label='number of correct MTs')
        ax.axis((0, self.num_steps*dt, -0.5, 4.5))
        plt.show()