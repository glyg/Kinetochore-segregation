# -*- coding: utf-8 -*-
from kt_simul.cluster import Launcher

if __name__=='__main__':

    PARAMFILE = "params.xml"
    MEASUREFILE = "measures.xml"
    result_path = "/media/thor/data/ktsimu"
    number_simu = 100

    l = Launcher(result_path,
                 number_simu,
                 paramfile = PARAMFILE,
                 measurefile = MEASUREFILE)

    l.run()
