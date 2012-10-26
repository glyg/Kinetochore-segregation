# -*- coding: utf-8 -*-

import argparse

from kt_simul.cluster import Launcher
from kt_simul.cluster.process import Process

if __name__ == '__main__':

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Launcher')
    parser.add_argument('--nsimu', "-n", type=int, default=10,
                        help='Number of simulations to launch (default = 10)')
    parser.add_argument("--path", "-p", type=str, required=True,
                        help='Directory to store results')
    parser.add_argument("--name", "-a", type=str, default="",
                        help='Name of the simulations')
    parser.add_argument("--eval", default=False,
                        action="store_true",
                        help='Launch pool evaluations')
    args = parser.parse_args()

    PARAMFILE = "params.xml"
    MEASUREFILE = "measures.xml"
    result_path = args.path
    number_simu = args.nsimu
    name = args.name
    launch_eval = args.eval
    ncore = 4

    l = Launcher(result_path,
                 number_simu,
                 name=name,
                 ncore=ncore,
                 paramfile=PARAMFILE,
                 measurefile=MEASUREFILE)

    l.run()

    if launch_eval:
        p = Process(results_path=l.results_path)
        resu = p.evaluate(groups=['attachment_state'], debug=True)
