#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This packages contains the kinetochore dynamics simulation

modules provided:
    - core : main simulation
    - gui: a subpackage with a GUI of the simulation
    - analysis : Data processing, analysis, batch scripts
"""
__all__ = ["core", "analysis", "gui"]

import logging
import os
import sys
import numpy

# Setup logging
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s:%(levelname)s: %(message)s',
                     datefmt='%Y/%m/%d-%H:%M:%S')
