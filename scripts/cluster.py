# -*- coding: utf-8 -*-

import sys
import argparse

from kt_simul.cluster import Launcher

if __name__ == '__main__':

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Launcher')
    parser.add_argument('--nsimu', "-n", type=int, default=10, help='Number of simulations to launch (default = 10)')
    parser.add_argument("--path", "-p", type=str, required=True, help='Directory to store results')
    args = parser.parse_args()

    PARAMFILE = "params.xml"
    MEASUREFILE = "measures.xml"
    result_path = args.path
    number_simu = args.nsimu
    ncore = 4

    l = Launcher(result_path,
                 number_simu,
                 ncore=ncore,
                 paramfile=PARAMFILE,
                 measurefile=MEASUREFILE)

    l.run()
