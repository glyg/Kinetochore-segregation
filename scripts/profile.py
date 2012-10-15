#!/usr/bin/env python
# encoding: utf-8
# filename: profile.py

import pstats, cProfile

import pyximport
pyximport.install()

from kt_simul import simul_spindle

m = simul_spindle.Metaphase()
cProfile.runctx("m.simul()", globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()
