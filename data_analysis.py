#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:uillaume Gay<elagachado AT  gmail DOT com>
## Commentary:


from pylab import *
from numpy import *
import os

# pcs1s = linspace(0., 1., 21)[::-1] # 14/02/2011
# auroras = logspace(-1.5, .2, 41)

pcs1s = linspace(0., 1., 41)[::-1] # 16/02/2011
auroras = logspace(-1.5, .2, 41)


# num_steps = 1001 # 14/02/2011
# num_ech = 10 #14/02/2011

# num_steps = 501 # 14/02/2011
# num_ech = 15 #14/02/2011

#num_steps = 401 # 14/02/2011
#num_ech = 30 #14/02/2011

# num_steps = 401 # 14/02/2011
# num_ech = 50 #14/02/2011

num_steps = 801 # 18/10/2011
num_ech = 50 #18/10/2011


time_step = 1.

Mk = 4
N = 3
num_proc = 4


class Explored(object):

    def __init__(self, num_steps, num_ech, base_name, default):
        
        self.pcs1s = load('%s_pcs1s.npy' %base_name)
        extra_pt = self.pcs1s[0] - self.pcs1s[1]
        self.pcs1s_L = hstack((extra_pt, self.pcs1s))
        
        self.auroras = load('%s_auroras.npy' %base_name)
        extra_pt = self.auroras[0] - self.auroras[1]
        self.auroras_L = hstack((extra_pt, self.auroras))

        
        self.num_steps = num_steps
        self.num_ech = num_ech
        self.timelapse = arange(num_steps) * time_step
        self.base_name = base_name
        self.default = default
        
        self.defect_keys = ('amphitelic', 'merotelic','monotelic',
                            'syntelic','unattached')

        self.mero_keys = ('corrected', 'cut',
                          'monotelic', 'syntelic')
 
        self.linestyles = {'amphitelic':'k--',
                           'merotelic':'r-',
                           'monotelic':'g-',
                           'syntelic':'purple',
                           'unattached':'b-'}
        self.def_tags = {'merotelic':3, 'syntelic':4, 'amphitelic':2,
                    'monotelic':1, 'unattached':0}

        self.space_shape = (self.pcs1s.size, self.auroras.size, num_ech, num_steps)
        self.sub_shape =  (self.pcs1s.size, self.auroras.size, num_steps)
        self.balance_shape = (self.pcs1s.size, self.auroras.size, 2 * Mk - 3,
                              num_ech,  num_steps)
        self.sub_bshape =  (self.pcs1s.size, self.auroras.size, 2 * Mk - 3, num_steps)

        self.trans_shape = (self.pcs1s.size, self.auroras.size, Mk + 1, num_ech, Mk + 1) 
        self.sub_tshape =  (self.pcs1s.size, self.auroras.size, Mk + 1, Mk + 1) 

    def __add__(self, other):

        if self.space_shape != other.space_shape :
            raise ValueError("Data shapes don't match. First shape: %s, second shape: %s " %(str(self.space_shape),
                                                                                             str(other.space_shape)))

        
        summed = Explored(self.num_steps, self.num_ech + other.num_ech, self.base_name, self.default)
        summed.total = (self.total + other.total) / 2.

        summed.defects = {}
        summed.mero_types = {}

        for defect in self.defects.keys():
            summed.defects[defect] = ( self.defects[defect] + other.defects[defect] ) / 2.
            
        for m_type in self.mero_types.keys():
            summed.mero_types[m_type] = ( self.mero_types[m_type] + other.mero_types[m_type] ) / 2.

        return summed


    def get_all(self):

        defects = {}
        balance = zeros(self.sub_bshape)
        transition = zeros(self.sub_tshape)
        
        total = zeros(self.sub_shape)
        for key in self.defect_keys:
            value = zeros(self.sub_shape)
            for i in r_[1:num_proc+1]:
                npy_name = '%s%idefect_%s_%s.npy' %(self.base_name, i, key, self.default)
                def_i = load(npy_name)
                def_i = def_i.reshape(self.space_shape).mean(axis = 2)
                value += def_i
                
            defects[key] = value/num_proc
            if 'amphi' not in key:
                total += value

        total_mero = zeros(self.sub_shape)
        mero_types={}
        for key in self.mero_keys:
            value = zeros(self.sub_shape)
            for i in r_[1:num_proc+1]:
                npy_name = '%s%imero_types_%s_%s.npy' %(self.base_name, i, key, self.default)
                def_i = load(npy_name)
                def_i = def_i.reshape(self.space_shape).mean(axis = 2)
                value += def_i
                
            mero_types[key] = value/num_proc
        
        for i in r_[1:num_proc+1]:
            # balance_name = '%s%ibalance_%s.npy' %(self.base_name, i, self.default)
            # bal = load(balance_name).reshape(self.balance_shape).mean(axis = 3)
            # balance += bal

            trans_name = '%s%itransition_%s.npy' %(self.base_name, i, self.default)
            tran = load(trans_name).reshape(self.trans_shape).mean(axis = 3)
            transition += tran


        self.balance = balance/num_proc
        self.total = total/num_proc
        self.defects = defects
        self.transition = transition/num_proc
        self.mero_types = mero_types

    def get_from_filelist(self, xml_file_list):
        
        from kt_simul import simul_spindle as sim
        from kt_simul import attachment_state as atst
        
        self.defects = {}
        self.mero_types = {}
        self.xml_file_list = xml_file_list


        
        for i in range(len(xml_file_list[::num_ech])):

            for j in range(num_ech):
                mp = sim.get_fromfile( xml_file_list[i * num_ech + j])
                if j == 0:
                    tmp_defects, were_defects = atst.defect_histories(mp.KD)
                    balance, tmp_mero_types = atst.balance_histories(mp.KD)
                else:
                    defects, were_defects = atst.defect_histories(mp.KD)
                    balance, mero_types = atst.balance_histories(mp.KD)
                    for defect in defects.keys():
                        tmp_defects[defect] += defects[defect] 
                    for mero_type in mero_types.keys():
                        tmp_mero_types[mero_type] += mero_types[mero_type] 
                    del defects, were_defects, balance, mero_types
                del mp

            if i == 0:
                for defect in tmp_defects.keys():
                    self.defects[defect] = tmp_defects[defect] / float(num_ech)
                for mero_type in tmp_mero_types.keys():
                    self.mero_types[mero_type] = tmp_mero_types[mero_type] / float(num_ech)
                
            else:
                for defect in tmp_defects.keys():
                    self.defects[defect] = hstack(( self.defects[defect],
                                                    tmp_defects[defect] / float(num_ech) ))
                for mero_type in tmp_mero_types.keys():
                    self.mero_types[mero_type] = hstack(( self.mero_types[mero_type],
                                                          tmp_mero_types[mero_type] / float(num_ech) ))
                
            
            print '%i / %i'%( ( i + 1 ) * num_ech , len(xml_file_list))

            
        self.total = zeros(self.sub_shape)
        for defect in tmp_defects.keys():
            self.defects[defect] = self.defects[defect].reshape(self.sub_shape)
            if 'amphitelic' not in defect:
                self.total += self.defects[defect]

        for mero_type in tmp_mero_types.keys():
            self.mero_types[mero_type] = self.mero_types[mero_type].reshape(self.sub_shape)

            
          

    def plot_defects(self, t_start, t_stop):

        sub_total = self.total[:,:,t_start:t_stop].mean(2)

        figure(1)

        pcolor(self.auroras_L, self.pcs1s_L, sub_total,
               vmin = 0., vmax = 3.)
        colorbar()
        xscale('log')
        axis((self.auroras_L.min(), self.auroras_L.max(), self.pcs1s_L.min(), self.pcs1s_L.max()))
        ylabel(r'Orientation effect')
        xlabel(r'$d_\alpha$')
        title('Total')
        figure(2)
        i = 1
        for key, defect in self.defects.items():
            if 'amphi' not in key:
                subplot(2,2,i)
                sub_defect = defect[:,:,t_start:t_stop].mean(2)
                pcolor(self.auroras_L, self.pcs1s_L, sub_defect,
                       vmin = 0., vmax = 3.)
                colorbar()
                xscale('log')
                ylabel(r'Orientation effect')
                xlabel(r'$d_\alpha$')
                axis((self.auroras_L.min(), self.auroras_L.max(), self.pcs1s_L.min(), self.pcs1s_L.max()))
                title(key)
                i+=1

    def plot_dynamics(self, pcs1, aurora):
    
        sub_total = self.total[pcs1,aurora,:]
        figure()
        plot(self.timelapse/60, sub_total, 'k-', label = 'Total', lw = 3)
        for key, defect in self.defects.items():
            sub_defect = defect[pcs1,aurora,:]
            plot(self.timelapse/60, sub_defect, self.linestyles[key],
                 label = key, lw = 1.5)

def renormalize(in_array, reference):

    out_array = in_array.copy().astype(float)
    out_array[reference != 0] /=  reference[reference != 0]
    #out_array[reference == 0] 
    return out_array

def plot_mero_types(explo, t_start, t_stop, optimum, p = 1, a_sample=5, b_sample=5):

    optimum = argmin(explo.total[-1,:,700:740].mean(-1))
    meros = explo.defects['merotelic']
    cor = explo.mero_types['corrected']
    cut = explo.mero_types['cut']
    synt = explo.mero_types['syntelic']
    total = cor + cut + synt
    
    cor = renormalize(explo.mero_types['corrected'], total)
    cut = renormalize(explo.mero_types['cut'], total )
    synt = renormalize(explo.mero_types['syntelic'], total)
    mono = renormalize(explo.mero_types['monotelic'], total)

    # figure()
    # for i, pcs1 in enumerate(explo.pcs1s[1:]):
        
    #     plot(explo.auroras, cor[i,:,t_start:t_stop].mean(-1),
    #          'go', alpha = 0.5)
    #     plot(explo.auroras, cut[i,:,t_start:t_stop].mean(-1),
    #          'ro', alpha = 0.5)
    #     plot(explo.auroras, synt[i,:,t_start:t_stop].mean(-1),
    #          'bo', alpha = 0.5)
    #     plot(explo.auroras, mono[i,:,t_start:t_stop].mean(-1),
    #          'o', c='grey', alpha = 0.5)
    figure(20)
    a_width = (log10(explo.auroras[::a_sample][1]/explo.auroras[::a_sample][0])) / 3.
    b_width = (explo.pcs1s[::b_sample][0] - explo.pcs1s[::b_sample][1]) / 3.
    for s, (defect, c) in enumerate( zip((cor, cut, synt), ('g','r','b'))):

        ys = defect[p,::a_sample,t_start:t_stop].mean(-1)*100
        ys_err = defect[p,::a_sample,t_start:t_stop].std(-1)*100
        bar(log10(explo.auroras[::a_sample])+(s-3/2.)*a_width, ys, yerr=ys_err, width=a_width,
                 color=c, lw=1, ecolor='k')

    plot(log10(explo.auroras[::a_sample]), explo.total[p,::a_sample,t_start:t_stop].mean(-1)*100/3.,
         'k-', lw=2)
    plot(log10(explo.auroras[::a_sample]), meros[p,::a_sample,t_start:t_stop].mean(-1)*100/3.,
         'ko-', lw=2)
    xticks(log10(explo.auroras[::a_sample]), ['%.2f' %a for a in explo.auroras[::a_sample]])
    xlabel(u'dα (µm)', fontsize = 12)
    ylabel('Fraction of merotelic chromosomes', fontsize = 12)

    axis((-3, 1., 0, 110))

    figure(21)
    for s, (defect, c) in enumerate( zip((cor, cut, synt), ('g','r','b'))):

        ys = defect[::b_sample,optimum,t_start:t_stop].mean(-1)*100
        ys_err = defect[::b_sample,optimum,t_start:t_stop].std(-1)*100

        bar(explo.pcs1s[::b_sample]+(s-3/2.)*b_width, ys, yerr=ys_err,width=b_width,
            color=c, lw=1, ecolor='k')

        # errorbar(explo.pcs1s, ys, yerr=ys_err,
        #          color='k', lw=1, ls = '.')
        
    plot(explo.pcs1s[::b_sample], explo.total[::b_sample,optimum,t_start:t_stop].mean(-1)*100/3.,
         'k-', lw=2)
    plot(explo.pcs1s[::b_sample], meros[::b_sample,optimum,t_start:t_stop].mean(-1)*100/3.,
         'ko-', lw=2)
    xlabel(u'β', fontsize = 12)
    ylabel('Fraction of merotelic chromosomes (%)', fontsize = 12)
    axis((-0.1, 1.5, 0, 110))


    # xscale('log')
    # xlabel(u'dα (µm)', fontsize = 12)
    # ylabel('Fraction of merotelic chromosomes', fontsize = 12)

    #    figure(4)

    


    # xs = explo.auroras[::5]
    # zs = explo.pcs1s[::5]
    # x_width = (log10(xs[1]/xs[0])) / 3.
    # z_width = (zs[0]- zs[1]) / 3.
    # for s, (defect, c) in enumerate( zip((cor, cut, synt), ('g','r','b'))):
    #     ys = defect[::5,::10, 350:370].mean(-1)*100
    #     ys_err = defect[::5,::10, 350:370].std(-1)*100
        
        
    #     #for n, z in enumerate(zs):
    #     figure(10)
        
    #     bar(log10(xs) + s * x_width, ys[0,:], width = x_width, color=c, alpha=0.8)
    #     errorbar(log10(xs)+ (s + 0.5)*x_width, ys[0,:],
    #              yerr=ys_err[0,:], color='k', ls='.')


    #     figure(11)
    #     bar(zs + s * z_width, ys[:,3], width = z_width, color=c, alpha=0.8)
    #     errorbar(zs+ (s + 0.5) * z_width, ys[:,3],
    #              yerr=ys_err[:,3], color='k', ls='.')


    # figure(10)
    # xlabel(u'dα (µm)', fontsize = 12)
    # ylabel('Fraction of merotelic chromosomes (%)', fontsize = 12)
    # xticks(log10(xs)+x_width*(s+1)/2., ['%0.2f'%x for x in xs])

    # figure(11)
    # xlabel(u'β', fontsize = 12)
    # ylabel('Fraction of merotelic chromosomes (%)', fontsize = 12)
    
    # xticks(zs+z_width*s/2., ['%0.2f'%z for z in zs])
        

    return cor, cut, synt

        
def plot_defects_slice(explo, pcs1, t_start, t_stop):

    sub_total = explo.total[pcs1,:,t_start:t_stop].mean(-1)
    #sub_balance = explo.balance[pcs1,:,:,t_start:t_stop].mean(-1)

    figure()
    plot(explo.auroras[:], sub_total*100/3., 'k-', label = 'Total', lw = 1)

    for key, defect in explo.defects.items():
        if 'amphi' not in key:
            sub_defect = defect[pcs1,:,t_start:t_stop].mean(-1)
            plot(explo.auroras[:], sub_defect*100/3., explo.linestyles[key],
                 label = key, lw = 1)

    ax = gca()

    for line in ax.lines:
        line.set_marker('o')
    

    # figure()
    # #cumul = zeros(sub_balance.shape[0])
    # cumul = explo.defects['unattached'][pcs1,:,t_start:t_stop].mean(-1)*100/3.
    # plot(explo.auroras, cumul, 'k-')
    # cumul += explo.defects['monotelic'][pcs1,:,t_start:t_stop].mean(-1)*100/3.
    # plot(explo.auroras, cumul, 'k-')
    
    # cumul += explo.defects['syntelic'][pcs1,:,t_start:t_stop].mean(-1)*100/3.
    # plot(explo.auroras, cumul, 'k-')

    # for i in range(sub_balance.shape[1]):
    #     cumul += sub_balance[:,i]*100/3.
    #     plot(explo.auroras, cumul, 'k-')

    # #xscale('log')
    # axis((explo.auroras[0], explo.auroras[-1], 0, 100))


    # figure()
    # total_aneu = explo.defects['syntelic'][pcs1,:,t_start:t_stop].mean(-1)*100/3.
    # total_aneu += sub_balance[:,0]*100/3. #balance = -2
    # total_aneu += sub_balance[:,1]*100/3. #balance = -1

    # plot(explo.auroras, total_aneu, color='pink', ls='-', marker='o', label='Aneuploidy')

    # cut = sub_balance[:,2]*100/3. #balance = 0
    # plot(explo.auroras, cut, color='red', ls='-', marker='o', label = 'Cut')
    
    # lag = sub_balance[:,3]*100/3. #balance = 1
    # lag += sub_balance[:,4]*100/3. #balance = 2
    # plot(explo.auroras, lag, color='g', ls='-', marker='o', label = 'Corrected')
    
    # #xscale('log')
    # axis((explo.auroras[0], explo.auroras[-1], 0, 100))
