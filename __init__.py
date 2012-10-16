#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This packages contains the kinetochore dynamics simulation

modules provided:
    - simul_spindle : main simulation
    - gui: a subpackage with a GUI of the simulation
    - explore_params : parameter exploration routines
"""
__all__ = ["core", "data_analysis", "gui"]


# from core import simul_spindle# Metaphase, get_fromfile

import logging

# Setup logging
logging.basicConfig( level = logging.DEBUG,
                     format = '%(asctime)s:%(levelname)s: %(message)s',
                     datefmt = '%Y/%m/%d-%H:%M:%S')


