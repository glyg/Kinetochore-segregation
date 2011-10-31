#!/usr/bin/python
# -*- coding: utf-8 -*-


### Adapted from Civelekoglu-Scholey et al. Biophys.J 90(11) 2006
### doi: 10.1529/biophysj.105.078691
try:
    from pylab import *
except RuntimeError:
    print 'No display..'
    pass

from scipy import linalg, column_stack, sparse
from numpy import any, all, arange, array, uint8, save, savetxt, mod, asarray, zeros
from xml.etree.ElementTree import Element, SubElement, tostring
import time
from Image import fromarray as image_fromarray
import os
import sys

import pyximport
pyximport.install()


## local imports
#from kt_simul import *
from kt_simul.spindle_dynamics import KinetoDynamics
from kt_simul.xml_handler import ParamTree, indent, ResultTree
from kt_simul.eval_simul import *
from kt_simul.dyninst import mean_attachment, delta_Mk

__all__ = ["Metaphase", "reduce_params", "paramfile", "measurefile", "get_fromfile"]

try:
    paramfile = os.path.join(os.path.dirname(__file__), 'params.xml')
    measurefile = os.path.join(os.path.dirname(__file__), 'measures.xml')
except NameError:
    paramfile = 'params.xml'
    measurefile = 'measures.xml'

    
measuretree = ParamTree(measurefile)
measuretree.create_dic(adimentionalized = False)
MEASURES = measuretree.dic


def reduce_params(paramtree, measuretree):

    ''' This functions changes the paramters so that the average behaviour
    complies with the measures.

    input:
    paramtree : ParamTree instance
    measures : a dictionary of measures
    paramtree is modified in place

    '''
    try:    
        params = paramtree.dic
    except AttributeError:
        paramtree.create_dic(adimentionalized = True)
        params = paramtree.dic
    try:
        measures = measuretree.dic
    except AttributeError:
        measuretree.create_dic(adimentionalized = False)
        measures = measuretree.dic

    try:
        poleward_speed = measures['poleward_speed'] 
        metaph_rate = measures['metaph_rate']  
        metaph_rate /= poleward_speed 
        anaph_rate = measures['anaph_rate']   
        anaph_rate /= poleward_speed
        mean_metaph_k_dist = measures['mean_metaph_k_dist']
        max_metaph_k_dist = measures['max_metaph_k_dist']        
        outer_inner_dist = measures['oi_dist']
        tau_o = measures['tau_o'] 
        tau_i = measures['tau_i'] 
        obs_d0 = measures['obs_d0']

    except KeyError:
        
        print "The measures dictionary should contain at least the following "
        print "keys "
        print MEASURES.keys()
        return 0


    Fk = params['Fk']
    fa = params['fa'] # 'free' attachement event frequency
    fd = params['fd'] # 'free' detachement event frequency
    aurora = params['aurora']
    N = params['N']
    Mk = params['Mk']
    kop = params['kop']

    #Let's go for the direct relations
    d0 = params['d0'] = obs_d0
    Vk = params['Vk'] = poleward_speed
    Vmz = params['Vmz'] = anaph_rate
    #Aurora modifies fd
    if aurora != 0:
        fd_eff = fa * aurora / mean_metaph_k_dist
    else:
        fd_eff = fd
    
    # alpha_mean = float(mean_attachment(fa/fd_eff) / Mk)
    alpha_mean = float(1/(1 + fd_eff/fa))
    #Take metaphase kt pair distance as the maximum one
    kappa = 2 * Mk / ( max_metaph_k_dist - d0 )
    params['kappa'] = kappa

    kop = alpha_mean * ( 1 + metaph_rate/2 ) / ( outer_inner_dist )
    params['kop'] = kop
    #Ensure we have sufficientely small time steps
    dt = params['dt']
    params['dt'] = min( tau_i/4 , tau_o/4. , params['dt'])
    if params['dt'] != dt:
        print 'Time step changed' 

    mus = params['mus']
    Fmz = N * Mk * alpha_mean * (1 + ( metaph_rate/2 ) *
                                 (1 + mus / ( N * Mk * alpha_mean ))
                                 / (1 -  metaph_rate / Vmz ))
    params['Fmz'] = Fmz
    mui = ( tau_i * kappa ) * Vk #Due to adimentionalization ### THIS SHOULD BE DONE A BIT MORE CARREFULY 
    params['mui'] = mui
    muo = ( tau_o * kop ) * Vk  #Due to adimentionalization 
    params['muo'] = muo 
    for key, val in params.items():
        paramtree.change_dic(key, val, write = False, verbose = False)




class Metaphase(object):
    
    """
    An instance of the Metaphase class is a wrapper around a whole simulation.



    Typical usage :
    >>> m = Metaphase()
    >>> m.simul()
    >>> m.show_trajs()
    >>> m.write_results('docstring_params.xml', 'docstring_results.xml')
     
    Public methods:
    simul() : runs the simulation
    show_trajs() : displays the trajectories
    show_one() : displays one of the chromosomes' trajectories
    show_one_article() : another display style
    evaluate() : calculates the characteristics of the simulation
    write_results() : saves the trajectories, evaluations and parameters
    write_asoctave() : save the trajectories with a different txt organisation
    get_ch(n) : returns chromosome n object
    get_2Dtraj_list : returns 2D trajectories of the spindle elements
    """

    def __init__(self,  duration=900, paramtree=None, measuretree=None,
                 paramfile=paramfile, measurefile=measurefile, plug=None, 
                 reduce_p=True):

        """
        Metaphase instanciation method
        
        Key-word arguments:
        duration  : the duration of the mitosis in seconds (defaults to 900)
        paramtree : a ParamTree instance containing the parameters for the
            simulation if paramtree is None, the parameters are read
            from the file paramfile. Defaults to None
        paramfile : a xml file to read the parameters from. Defaults to the
            file params.xml in the module's directory. Other parameter files
            can be produced by editing and changing the default one
        reduce_p  : if True, changes the parameters according to the measures
            so that the simulation average behaviour complies with
            the data in the measures dictionary
        measurefile : a xml file to read the measures from. Defaults to the
            file measures.xml in the module's directory. Other parameter files
            can be produced by editing and changing the default one
        measures  : a dictionary containing the observed characteristics
            of the mitosis e.g. metaphase spindle elongation rate, etc.
            Defaults to the data provided by the file measures.xml in
            the module's directory
        """

        if paramtree is None:
            self.paramtree = ParamTree(paramfile)
            self.paramtree.create_dic()
        else:
            self.paramtree = paramtree

        if measuretree is None:
            self.measuretree = ParamTree(measurefile)
            self.measuretree.create_dic(adimentionalized = False)
        else:
            self.measuretree = measuretree
            self.measuretree.create_dic(adimentionalized = False)
            
        if reduce_p:
            reduce_params(self.paramtree, self.measuretree)
            
        self.KD = KinetoDynamics(self.paramtree.dic, plug = plug)
        dt = self.KD.params['dt']
        self.nb_steps = int(duration/dt)
        self.KD.anaphase = False
        self.timelapse = arange(0, duration + dt, dt)
        self.report = []
        self.delay = -1



    def __str__(self):

        lines = []
        lines.append('Metaphase class')
        try:
            lines.append('Parameters:')
            for line in str(self.paramtree.dic).split(','):
                lines.append(line)
        except AttributeError:
            pass

        try:
            lines.append('Measures:' )
            for line in str(self.measuretree.dic).split(','):
                lines.append(line)
        except AttributeError:
            pass

        try:
            lines.append('')
            lines.append('Observations:')
            for line in str(self.observations).split(','):
                lines.append(line)
        except AttributeError:
            lines.append('Not yet evaluated')
            pass

        return '\n'.join(lines)
            
    
    def _one_step(self):

        if not self.KD.anaphase:
            self.KD.plug_unplug()
        A = self.KD.calcA()
        b = - self.KD.calcb()
        #Linalg solves A.x = b and not A.x + b = 0 !!!!! (F**cking shit)
        speeds = linalg.solve(A,b)

        # sparse.linalg.use_solver( useUmfpack = False )
        # A = A.astype('f')
        # speeds = sparse.linalg.spsolve(A,b)


        self.KD.position_update(speeds)#[0])

    def simul(self, movie = False, ablat = None): 

        ''' The simulation main loop. 
        
        Keyword arguments:
        movie: if True, runs make_movie (default False)
        ablat: timepoint at which ablation takes place. If None (default)
        no ablation is performed

        '''

        N = self.KD.params['N']
        Mk = self.KD.params['Mk']
        dt = self.KD.params['dt']
        kappa = self.KD.params['kappa']
        mus = self.KD.params['mus']


        for t in self.timelapse[1:]:
            # Ablation test
            if ablat is not None and ablat == t:
                self.ablation(pos = self.KD.spbL.pos)
            # Anaphase transition ?
            if self.anaphase_test(t):
                self.toa_test(t)

            self._one_step()
            
            if movie and t/dt % 15 == 0: #One picture every 15 time point
                self._make_movie(t, (50, 100,  3))



        self.KD.params['kappa'] = kappa
        self.KD.delay = self.delay - 1
        s = "delay = %2d seconds" % self.delay
        self.report.append(s)
        k_dist = metaph_kineto_dist(self.KD)
        s = "Mean Kt - Kt distance: %.3f m" % k_dist[0]
        self.report.append(s)
        #self.evaluate()
        #         for l in self.report:
        #             print l

    def toa_test(self, t):
        
        for n, ch in self.KD.chromosomes.items():    
            if ch.anaphase_switch[0] > 0 and ch.right_toa == 0:
                ch.right_toa = t
                s = "Right kt of chromosome %d reached the pole at time %.2f" %(n,t)
                self.report.append(s)
            if ch.anaphase_switch[1] > 0 and ch.left_toa == 0:
                ch.left_toa = t
                s = "Left kt of chromosome %d reached the pole at time %.2f" %(n,t)
                self.report.append(s)


    def anaphase_test(self, t):
        
        transition = int(self.KD.params['trans'])
        if self.KD.anaphase:
            return True
            
        if transition <= t and self._plug_checkpoint():
            if self.delay == -1 :
                self.delay = t - transition 
                #Then we just get rid of cohesin
                self.KD.params['kappa'] = 0.
                self.KD.B_mat = self.KD.write_B()
                self.nb_mero = self._mero_checkpoint()
                if self.nb_mero:
                    s = 'There were %d merotelic MT at anaphase onset' %self.nb_mero
                    self.report.append(s)
                self.KD.anaphase = True
                self.KD.test_anaphase_switch()
                return True

        return False


    def ablation(self, pos = None):
        if pos == None: pos = self.KD.spbR.pos
        if not self.KD.spbL.pos <= pos <= self.KD.spbR.pos:
            print 'Bad shot, same player play again!'
            return 

        self.KD.params['Fmz'] = 0.
        self.KD.params['fa'] = 0.
        self.KD.params['fd'] = 0.
        for ch in self.KD.chromosomes.values():
            (right_pluged, left_pluged) = ch.pluged()
            (right_mero, left_mero) = ch.mero()
            if pos > ch.rightpos:
                for rplug in ch.rplugs.values():
                    rplug.plug = min(0, rplug.plug)
                for lplug in ch.lplugs.values():
                    lplug.plug = max(0, lplug.plug)
                #ch.pluged = (right_pluged, left_pluged)
                #ch.mero = (0, 0)

            elif ch.rightpos >= pos > ch.leftpos:
                for rplug in ch.rplugs.values():
                    rplug.plug = max(0, rplug.plug)
                for lplug in ch.lplugs.values():
                    lplug.plug = max(0, lplug.plug)
                #ch.pluged = (0, left_pluged)
                #ch.mero = (right_mero, 0)

            elif pos < ch.leftpos:
                for lplug in ch.lplugs.values():
                    lplug.plug = min(0, lplug.plug)
                for rplug in ch.rplugs.values():
                    rplug.plug = max(0, rplug.plug)


            (right_pluged, left_pluged) = ch.pluged()
            (right_mero, left_mero) = ch.mero()
            print (right_pluged, left_pluged)
            print (right_mero, left_mero)
                #ch.pluged = (right_pluged, 0)
                #ch.mero = (0,left_mero)

    def _plug_checkpoint(self):
        '''returns True if all chromosomes are pluged by at least
        one kMT, False otherwise
        '''
        
        N = self.KD.params['N']
        sac = self.KD.params['sac']
        # if sac == 0:
        #     return True
        
        for ch in self.KD.chromosomes.values() :
            if not all(ch.pluged()) :
                ch.active_sac = 1
                #print "active checkpoint"
                return False
            else:
                ch.active_sac = 0
                
        return True

    def _mero_checkpoint(self):
        '''returns the total number of merotellic kT
        '''
        N = self.KD.params['N']
        nb_mero = 0
        for ch in self.KD.chromosomes.values() :
            if any(ch.mero) :
                nb_mero += sum(ch.mero())
                #print "active checkpoint"
        return nb_mero


    def _mplate_checkpoint(self):
        '''
        returns True if each kinetochore is in the proper half of the spindle
        '''
        
        N = self.KD.params['N']
        nb_out = 0
        for ch in self.KD.chromosomes.values() :
            ktR = ch.rightpos
            ktL = ch.leftpos
            if ktR <= 0 or ktL >= 0:
                return True#False
        return True


    def evaluate(self):

        '''
        passes all the evaluations in eval_simul.py
        results are stored in the self.observations dictionnary
        '''
        if len(self.KD.spbR.traj) < 2:
            print "No simulation was ran ... exiting"
            return 0

        self.observations = {'anaphase_rate':anaphase_rate(self.KD),
                             'metaph_rate': metaph_rate(self.KD),
                             'mean_metaph_k_dist': metaph_kineto_dist(self.KD),
                             'pitch': auto_corel(self.KD, smooth = 5.),
                             'poleward_speed': poleward_speed(self.KD),
                             'kt_rms_speed':kt_rms_speed(self.KD),
                             'times_of_arrival': time_of_arrival(self.KD),
                             'pluged_stats': pluged_stats(self.KD)}


    def show_trajs(self, axes = None): 
        ''' Plot the different trajectories
        '''
        N = self.KD.params['N']
        Mk = self.KD.params['Mk']
        dt = self.KD.params['dt']

        if axes == None:
            fig = figure(1)
            axes = gca()

        spbRtraj = array(self.KD.spbR.traj)
        spbLtraj = array(self.KD.spbL.traj)
        try:
            axes.plot(self.timelapse, spbRtraj, color = 'r', ls = '-', lw=1)
            axes.plot(self.timelapse, spbLtraj, color = 'r', ls = '-', lw=1)
        except AttributeError:
            print 'please run the simulation before'
            axes.close()
            return 0

        fmt_list = ['g-', 'b-', 'm-']
        for n in range(N):
            ch = self.KD.chromosomes[n]
            right_traj = array(ch.righttraj)
            left_traj = array(ch.lefttraj)
            fmt = fmt_list[mod(n,3)]
            line1 = axes.plot(self.timelapse,right_traj,fmt, alpha = 0.5)
            line2 = axes.plot(self.timelapse,left_traj,fmt, alpha = 0.5)

            if n == 0:
                line1[0].set_alpha(1.)
                line2[0].set_alpha(1.)
            # if n == 0:
            #     axes.plot(self.timelapse[::int(15/dt)], right_traj[::int(15/dt)],
            #               mec, mfc = 'None')
            #     axes.plot(self.timelapse[::int(15/dt)], left_traj[::int(15/dt)],
            #               'go', mfc = 'None')
                
        #axes.axis([0, self.nb_steps * dt, -2.2, 2.2 ])
        axes.set_xlabel('Time (seconds)', fontsize = 'small')
        axes.set_ylabel(u'Distance from center (µm)', fontsize = 'small')
        
        show()

    def write_asoctave(self, xmlfname = "results.xml",
                       datafname = "oct_datas.txt"):

        '''
        Results file compatible with the tools developped to caracterize
        real life movies under octave
        '''
        
        N = self.KD.params['N']
        self.write_results(xmlfname = xmlfname)

        spbRx = array(self.KD.spbR.traj)
        ys = zeros(spbRx.shape)
        ts = arange(spbRx.shape[0])
        idxs = ones(spbRx.shape)

        out_array = column_stack((spbRx, ys, ts, idxs))
        
        spbLx = array(self.KD.spbL.traj)
        idxs += 1 
        out_array = append(out_array, column_stack((spbLx, ys, ts, idxs)),
                           axis = 0)

        for ch in self.KD.chromosomes.values():

            idxs += 1
            chRx = array(ch.righttraj)
            out_array = append(out_array, column_stack((chRx, ys, ts, idxs)),
                               axis = 0)
            idxs += 1            
            chLx = array(ch.lefttraj)
            out_array = append(out_array, column_stack((chLx, ys, ts, idxs)),
                               axis = 0)
            
        dataout = file(datafname, 'w+')
        dataout.write("# desciptor: "+xmlfname+"\n")
        write_array(dataout, out_array, separator=' ', linesep='\n')


        
        
    def write_results(self, xmlfname = "results.xml", datafname = "datas.npy"):
        '''
        Saves the results of the simulation in two files with the parameters,
        the measures and observations in one file and the trajectories in
        an other file

        Keyword arguments:
        xmlfname : name of the xml file where parameters and observations
        will be written
        datafname : name of the file where the trajectories will be written
        file type is determined by the file suffix:
         - *.npy : data are stored in numpy"s binary format
                   (less portable but quite efficient)
         - *.txt : simple text
         - *.txt.gz : text files compressed transparently
        Any other suffix will be saved as plain text
         
        Column index for each trajectory is an
        attribute of the corresponding element in the xml file.
        '''

        N = self.KD.params['N']
        Mk = self.KD.params['Mk']
        dt = self.KD.params['dt']
        if not hasattr(self, 'observations'):
            self.evaluate()
        

        chromosomes = self.KD.chromosomes.values()
        wavelist=[]
        out = file(xmlfname, 'w+')
        out.write('<?xml version="1.0"?>\n')
        today = time.asctime()
        experiment = Element("experiment", date=today, datafile=datafname)
        #the parameters used
        experiment.append(self.paramtree.root)
        #the measures used
        experiment.append(self.measuretree.root)
        #right SPB
        spbR = SubElement(experiment, "trajectory", name = "rightspb",
                          column='0', units='mu m')
        SubElement(spbR, "description").text="right spb trajectory"
        #we'll systematicaly append the corresponding trajectory to the list, 
        #so that it's placed according to the "column" attribute
        spbRtraj = array(self.KD.spbR.traj)
        wavelist.append(spbRtraj)
        #left SPB
        spbL = SubElement(experiment, "trajectory", name = "leftspb",
                          column='1', units='mu m')
        SubElement(spbL, "description").text="left spb trajectory"
        spbLtraj = array(self.KD.spbL.traj)
        wavelist.append(spbLtraj)

        col_num = 2
        #chromosomes
        for n, ch in enumerate(chromosomes):
            # Adding pluged_history and mero_history
            chp = array(ch.pluged_history)
            rp = chp[:, 0]
            lp = chp[:, 1]
            chm = array(ch.mero_history)
            rm = chm[:, 0]
            lm = chm[:, 1]

            rch = SubElement(experiment, "trajectory", name="rightkineto",
                             index = str(n),
                             column=str(col_num), units='mu m')
            SubElement(rch, "description").text="chromosome %i right kinetochore trajectory" %n
            wavelist.append(array(ch.righttraj)); col_num += 1
            
            rchp = SubElement(experiment, "numberpluged",
                              name="rightkineto", index = str(n),
                              column=str(col_num))
            wavelist.append(rp); col_num += 1
            rchp = SubElement(experiment, "numbermero",
                              name="rightkineto", index = str(n),
                              column=str(col_num))
            wavelist.append(rm); col_num += 1
            
            lch = SubElement(experiment, "trajectory",
                             name="leftkineto", index=str(n), 
                             column=str(col_num), units='mu m')
            SubElement(lch, "description").text="chromosome %s left kinetochore trajectory" %n
            wavelist.append(array(ch.lefttraj)); col_num += 1
            lchp = SubElement(experiment, "numberpluged", name="leftkineto",
                              index = str(n),
                              column=str(col_num))
            wavelist.append(lp); col_num += 1
            lchm = SubElement(experiment, "numbermero", name="leftkineto",
                              index = str(n),
                              column=str(col_num))
            wavelist.append(lm); col_num += 1
            #Plug Sites
            for m, rplug in enumerate(ch.rplugs.values()):
                rpt_traj = array(rplug.traj)
                rpt = SubElement(experiment, "trajectory",
                                 name="rightplugsite", index = str((n, m)),
                                 column=str(col_num), units='mu m')
                wavelist.append(rpt_traj); col_num += 1
                rpt_plug = array(rplug.state_hist)
                rpp = SubElement(experiment, "state",
                                 name="rightplugsite", index = str((n, m)),
                                 column=str(col_num), units='')
                wavelist.append(rpt_plug); col_num += 1
            
            for m, lplug in enumerate(ch.lplugs.values()):
                lpt_traj = array(lplug.traj)
                lpt = SubElement(experiment, "trajectory",
                                 name="leftplugsite", index = str((n, m)),
                                 column=str(col_num), units='mu m')
                wavelist.append(lpt_traj); col_num += 1
                lpt_plug = array(lplug.state_hist)
                lpp = SubElement(experiment, "state",
                                 name="leftplugsite", index = str((n, m)),
                                 column=str(col_num), units='')
                wavelist.append(lpt_plug); col_num += 1

        #Observations
        obs_elem = SubElement(experiment, "observations")
        if not hasattr(self, 'observations'):
            self.evaluate()
        for key, val in self.observations.items():
            obs = SubElement(obs_elem, key).text = str(val)

        #Now we write down the whole experiment XML element
        indent(experiment)
        out.write(tostring(experiment))
        out.close()
        #And the numbers in the file datafname
        dataout = file(datafname, 'w+')
        data = column_stack(wavelist)
        if datafname.endswith('.npy'):
            save(dataout, data)
        else:
            dataout.write("# desciptor: "+xmlfname+"\n")
            savetxt(dataout, data, delimiter=' ')

    def _make_movie(self, t, imsize):

        N = self.KD.params['N']
        Mk = self.KD.params['Mk']

        try:
            os.mkdir("movie")
        except OSError:
            pass
        
        imname = "movie/spindle%03i.tif" % t
        if os.path.isfile(imname):
            os.remove(imname)

        imarray = zeros(imsize, uint8)

        xspbR = scale(self.KD.spbR.pos, imsize[1])
        yspbR = scale(0, imsize[0])
        xspbL = scale(self.KD.spbL.pos, imsize[1])
        yspbL = scale(0, imsize[0])
        imarray[yspbR, xspbR, 0] += 128
        imarray[yspbL, xspbL, 0] += 128

        chromosomes = self.KD.chromosomes.values()
        for n, ch in enumerate(chromosomes):
            xchR = scale(ch.rightpos, imsize[1])
            ychR = scale(0.1 * n - 0.1, imsize[0])
            xchL = scale(ch.leftpos, imsize[1])
            ychL = scale(0.1 * n - 0.1, imsize[0])
            imarray[ychR, xchR, 1] += 42
            imarray[ychL, xchL, 1] += 42

        if 80< t <100:
            self.testim =  imarray
            print xchR, ychR, t
        im = image_fromarray(imarray, 'RGB')
        im.save(imname)

        del im



    def show_one(self, n = 0, fig = None):

        """ Shows one chromosome's trajectory and plug state """

        dt = self.KD.params['dt']
        if fig == None:
            fig = figure()

        ch = self.KD.chromosomes[n]
        fig.clear()
        fig.add_subplot(312)
        fig.gca().plot(self.timelapse, ch.righttraj, 'b', lw = 1.5)
        fig.gca().plot(self.timelapse, ch.lefttraj, 'k', lw = 1.5)
        rps = ch.rplugs
        for ps in rps.values():
            fig.gca().plot(self.timelapse, ps.traj, 'g')
        lps = ch.lplugs
        for ps in lps.values():
            fig.gca().plot(self.timelapse, ps.traj, 'purple')
        mero_hist = asarray(ch.mero_history)
        pluged_hist = asarray(ch.pluged_history)

        
        fig.add_subplot(311)
        fig.gca().plot(self.timelapse, mero_hist[:,0], 'r',
                       label= 'number of merotellic MTs')
        fig.gca().plot(self.timelapse, pluged_hist[:,0], 'g',
                       label= 'number of pluged MTs')
        fig.gca().axis((0, self.nb_steps*dt, -0.5, 4.5))

        fig.add_subplot(313)
        fig.gca().plot(self.timelapse, mero_hist[:,1], 'r',
                       label= 'number of merotelic MTs')
        fig.gca().plot(self.timelapse, pluged_hist[:,1], 'g',
                       label= 'number of pluged MTs')
        fig.gca().axis((0, self.nb_steps*dt, -0.5, 4.5))
        show()

    def show_one_article(self, n = 0, fig = None):

        """ Shows one chromosome's trajectory and plug state """        

        dt = self.KD.params['dt']
        if fig == None:
            fig = figure()

        ch = self.KD.chromosomes[n]
        fig.clear()
        ax_traj = fig.add_subplot(312)
        spbRtraj = array(self.KD.spbR.traj)
        spbLtraj = array(self.KD.spbL.traj)

        ax_traj.plot(self.timelapse, spbRtraj, 'r-')
        ax_traj.plot(self.timelapse, spbLtraj, 'r-')
        
        ax_traj.plot(self.timelapse, ch.righttraj, 'g-', lw = 1., alpha = 1.)
        ax_traj.plot(self.timelapse, ch.lefttraj, 'g-', lw = 1., alpha = 1.)
        rps = ch.rplugs
        for ps in rps.values():
            ax_traj.plot(self.timelapse, ps.traj, 'purple', lw = 0.5,
                         alpha = 0.5)
        lps = ch.lplugs
        for ps in lps.values():
            ax_traj.plot(self.timelapse, ps.traj, 'purple', lw = 0.5,
                         alpha = 0.5)

        axis((0, self.nb_steps*dt, -2, 2))
        yticks(range(-2,4,2))
        xticks(range(0, 850 , 240), ['0', '4', '8', '12'])
        ylabel(u'Position (µm)', fontsize = 'small')
        xlabel('Time (min)', fontsize = 'small')

        
        mero_hist = array(ch.mero_history)
        pluged_hist = array(ch.pluged_history)

        ax_plug = fig.add_subplot(311)
        ax_plug.plot(self.timelapse, mero_hist[:,0], 'r-',
                     label= 'Number of merotellic kMTs')
        ax_plug.plot(self.timelapse, pluged_hist[:,0], 'g-',
                     label= 'Number of synthelic kMTs')
        ax_plug.axis((0, self.nb_steps*dt, -0.5, 4.5))
        xticks(range(0, 850 , 240), (''))
        ylabel('kMTs state', fontsize = 'small')
        yticks(arange(0,5,2))

        ax_plug2 = fig.add_subplot(313)
        ax_plug2.plot(self.timelapse, mero_hist[:,1], 'r-',
                      label= 'number of merotellic MTs')
        ax_plug2.plot(self.timelapse, pluged_hist[:,1], 'g-',
                      label= 'number of synthelic MTs')
        ax_plug2.axis((0, self.nb_steps*dt, -0.5, 4.5))
        xticks(range(0, 850 , 240), (''))
        ylabel('kMTs state', fontsize = 'small')
        ax_plug.set_aspect(10, anchor = 'S')
        show()

        return fig
        
    def get_ch(self, n = 0):

        return self.KD.chromosomes[n]

    def get_2Dtraj_list(self, ang_noise = 5e-3, pos_noise = 1e-2, cen = None):
        
        """returns a list of 2D trajectories, adding a random mouvement
        to the whole spindle

        keyword arguments:
        ang_noise : angular variation in radians per seconds
        pos_noise : center of mass displacement in µm/s
        cen: 0, 1 or 2, returns the trajectory of a simulated cen marker
        attached to the chromosome number cen. Adds a random coiled coil
        movement around the kinetochore position with sigma=0.13 µm
        if None (default), returns all 6 trajectories plus the SPBs
        """

        dt = self.KD.params['dt']
        trajectories = []
        ang_noise *= dt # rad.s^-1
        pos_noise *=  dt # µm.s^-1
        xspbR = array(self.KD.spbR.traj)
        xspbL = array(self.KD.spbL.traj)
        yspbR = zeros(xspbR.shape)
        yspbL = zeros(xspbL.shape)
        trajectories.append(column_stack((yspbR, xspbR)))
        trajectories.append(column_stack((yspbL, xspbL)))

        traj_length = xspbR.shape[0]

        if cen is not None:
            ch = self.get_ch(cen)
            xrt = array(ch.righttraj) + normal(0, scale = 0.13,
                                               size = traj_length)  
            yrt = zeros(xrt.shape) + normal(0, scale = 0.13,
                                            size = traj_length)
            trajectories.append(column_stack((yrt, xrt)))
            xlt = array(ch.lefttraj) + normal(0, scale = 0.13,
                                               size = traj_length)
            ylt = zeros(xlt.shape)+ normal(0, scale = 0.13,
                                           size = traj_length)
            trajectories.append(column_stack((ylt, xlt)))
            
        else:
            for n in range(3):
                ch = self.get_ch(n)
                xrt = array(ch.righttraj)
                yrt = zeros(xrt.shape) + (1 - n) * 0.3
                trajectories.append(column_stack((yrt, xrt)))
                xlt = array(ch.lefttraj)
                ylt = zeros(xlt.shape) + (1 - n) * 0.3
                trajectories.append(column_stack((ylt, xlt)))

        
        xcs = []
        ycs = []
        rots = []
        xd, yd, thetad = (0.,)*3

        for n in range(xspbR.size):
            xd += normal(0, scale = pos_noise)
            yd += normal(0, scale = pos_noise)
            thetad += normal(0, scale = ang_noise)
            xcs.append(xd)
            ycs.append(yd)
            rots.append([[cos(thetad), sin(thetad)],
                         [cos(thetad), -sin(thetad)]])
        
        for traj in trajectories:
            traj += column_stack((array(ycs), array(xcs)))
            n = 0
            for pos, rot in zip(traj, rots):
                new_pos = dot(pos, rot)
                traj[n] = new_pos
                n += 1

        return trajectories

def get_fromfile(xmlfname = "results.xml"):

    '''re-creates the Kinetochores Dynamics datas structure
    i.e a dictionnary of "waves" for each element
    returns a Metaphase instance
    '''
    restree = ResultTree(xmlfname)
    param_root = restree.root.find('parameters')
    paramtree = ParamTree(root = param_root)
    paramtree.create_dic()
    params = paramtree.dic

    measure_root = restree.root.find('measures')
    measuretree = ParamTree(root = measure_root)
    measuretree.create_dic(adimentionalized = False)
    metaphase = Metaphase(paramtree = paramtree,
                          measuretree = measuretree)
    
    traj_matrix = restree.get_all_trajs()
    pluged_matrix = restree.get_all_pluged()
    mero_matrix = restree.get_all_mero()
    state_hist_matrix = restree.get_all_plug_state()
    KD = KinetoDynamics(params)
    KD.spbR.traj = traj_matrix[:,0]
    KD.spbL.traj = traj_matrix[:,1]
    N = params['N']
    Mk = params['Mk']
    col_num = 2
    state_num = 0
    for n in range(N):
        KD.chromosomes[n].righttraj = traj_matrix[:, col_num]
        col_num += 1
        KD.chromosomes[n].lefttraj = traj_matrix[:, col_num]
        col_num += 1
        KD.chromosomes[n].pluged_history = (pluged_matrix[:, n*2: n*2 + 2])
        KD.chromosomes[n].mero_history = (mero_matrix[:, n*2 : n*2 + 2])
        for m in range(Mk):
            KD.chromosomes[n].rplugs[m].traj = traj_matrix[:, col_num]
            col_num += 1
            KD.chromosomes[n].rplugs[m].state_hist = state_hist_matrix[:, state_num]
            state_num += 1
        for m in range(Mk):
            KD.chromosomes[n].lplugs[m].traj = traj_matrix[:, col_num]
            col_num += 1
            KD.chromosomes[n].lplugs[m].state_hist = state_hist_matrix[:, state_num]
            state_num += 1
    metaphase.KD = KD

    return metaphase





def scale(x, size, pix_size = 0.0645):
    '''
    Scale the position x on a line of size "size" from microns to pixels
    origin is put in the center of the line
    '''
    return int(x / pix_size + size / 2)
    
def circ_box(xy, rad):
    return (xy[0] - rad, xy[1] - rad , xy[0] + rad , xy[1] + rad)
    
    
