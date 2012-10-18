# -*- coding: utf-8 -*-
"""
"""

import time
import logging

import numpy as np

from xml.etree.ElementTree import Element, SubElement, tostring

from kt_simul.io.xml_handler import ParamTree, indent, ResultTree


class SimuIO():
    """
    """

    def __init__(self, meta_instance=None):
        """
        :param meta_instance: Should have been already perform a simulation
        :type meta_instance: Metaphase instance
        """

        if meta_instance:
            self.meta = meta_instance
            self.KD = self.meta.KD
            self.timelapse = self.meta.timelapse
            self.num_steps = self.meta.num_steps
            self.paramtree = self.meta.paramtree
            self.measuretree = self.meta.measuretree
            self.observations = self.meta.observations

    def save(self, xmlfname="results.xml",
                datafname="data.npy"):
        """
        Saves the results of the simulation in two files
        with the parameters, measures and observations in one file
        and the trajectories in the other.

        :param xmlfname: name of the xml file where parameters and observations will be written
        :type xmlfname: string, optional

        :param datafname: name of the file where the trajectories will be written
            file type is determined by the file suffix:
                 - *.npy : data are stored in numpy's binary format
                           (less portable but quite efficient)
                 - *.txt : simple text
                 - *.txt.gz : text files compressed transparently
        :type datafname: string, optional

        Any other suffix will be saved as plain text. Column index for
        each trajectory is an attribute of the corresponding element
        in the xml file.

        .. todo:: This function should be chopped off some how, it's messy

        """
        if not self.observations:
            self.meta.evaluate()
            self.observations = self.meta.observations

        chromosomes = self.KD.chromosomes
        wavelist = []
        out = file(xmlfname, 'w+')
        out.write('<?xml version="1.0"?>\n')
        today = time.asctime()
        experiment = Element("experiment", date=today,
                            datafile=datafname)
        experiment.append(self.paramtree.root)
        experiment.append(self.measuretree.root)

        #right SPB
        spbR = SubElement(experiment, "trajectory", name="rightspb",
                          column='0', units='mu m')
        SubElement(spbR, "description").text="right spb trajectory"
        spbRtraj = np.array(self.KD.spbR.traj)
        wavelist.append(spbRtraj)

        #left SPB
        spbL = SubElement(experiment, "trajectory", name="leftspb",
                          column='1', units='mu m')
        SubElement(spbL, "description").text="left spb trajectory"
        spbLtraj = np.array(self.KD.spbL.traj)
        wavelist.append(spbLtraj)

        col_num = 2
        #chromosomes
        for n, ch in enumerate(chromosomes):
            rch = SubElement(experiment, "trajectory",
                            name="centromereA",
                            index = str(n),
                            column=str(col_num),
                            units='mu m')
            text="chromosome %i centromere A trajectory" % n
            SubElement(rch, "description").text = text
            wavelist.append(ch.cen_A.traj)
            col_num += 1

            SubElement(experiment, "numbercorrect", name="centromereA",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.correct_history[:, 0])
            col_num += 1
            SubElement(experiment, "numbererroneous",
                       name="centromereA", index=str(n),
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

    def read(self, xmlfname="results.xml"):
        """
        Creates a simul_spindle.Metaphase from a XML file.

        :param xmlfname: The xml file where results from the existing
                            simulation are.
        :type xmlfname: string, optional
        :return: a Metaphase instance
        """

        restree = ResultTree(xmlfname)
        param_root = restree.root.find('parameters')
        paramtree = ParamTree(root = param_root)
        params = paramtree.relative_dic
        measure_root = restree.root.find('measures')
        measuretree = ParamTree(root = measure_root,
                                adimentionalized = False)
        metaphase = Metaphase(paramtree = paramtree,
                              measuretree = measuretree)

        traj_matrix = restree.get_all_trajs()
        correct_matrix = restree.get_all_correct()
        erroneous_matrix = restree.get_all_erroneous()
        state_hist_matrix = restree.get_all_plug_state()
        KD = KinetoDynamics(params)
        KD.spbR.traj = traj_matrix[:, 0]
        KD.spbL.traj = traj_matrix[:, 1]
        Mk = int(params['Mk'])
        col_num = 2
        state_num = 0
        for n, ch in enumerate(KD.chromosomes) :
            ch.cen_A.traj = traj_matrix[:, col_num]
            col_num += 1
            ch.cen_B.traj = traj_matrix[:, col_num]
            col_num += 1
            ch.erroneous_history = (erroneous_matrix[:, n*2 : n*2 + 2])
            ch.correct_history = (correct_matrix[:, n*2 : n*2 + 2])
            for plugsite in ch.cen_A.plugsites:
                plugsite.traj = traj_matrix[:, col_num]
                col_num += 1
                plugsite.state_hist = state_hist_matrix[:, state_num]
                state_num += 1
            for plugsite in ch.cen_B.plugsites:
                plugsite.traj = traj_matrix[:, col_num]
                col_num += 1
                plugsite.state_hist = state_hist_matrix[:, state_num]
                state_num += 1
        metaphase.KD = KD

        return metaphase
