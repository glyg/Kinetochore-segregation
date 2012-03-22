#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This packages contains the kinetochore dynamics simulation

modules provided:
    - simul_spindle : main simulation
    - gui: a subpackage with a GUI of the simulation
    - eval_simul: evaluation routines
    - explore_params : parameter exploration routines
"""
__all__ = ["simul_spindle", "explore_params"]


import core.simul_spindle 
import analysis.eval_simul
#import gui


