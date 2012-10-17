# -*- coding: utf-8 -*-

import sys

from kt_simul.cluster import Launcher

if __name__ == '__main__':

    PARAMFILE = "params.xml"
    MEASUREFILE = "measures.xml"
    result_path = sys.argv[1]
    number_simu = 1000
    ncore = 4

    l = Launcher(result_path,
                 number_simu,
                 ncore=ncore,
                 paramfile=PARAMFILE,
                 measurefile=MEASUREFILE)

    l.run()
