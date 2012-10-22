# -*- coding: utf-8 -*-

import argparse

from kt_simul.cluster.process import Process

if __name__ == '__main__':

    # Arguments parser
    parser = argparse.ArgumentParser(description='KtSimu Pool Evaluator')
    parser.add_argument("--path", "-p", type=str, required=True,
                        help='Directory where simulations are stored')
    args = parser.parse_args()

    results_path = args.path

    p = Process(results_path=results_path)
    resu = p.evaluate(groups=['attachment_state'], debug=True)

