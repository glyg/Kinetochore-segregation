#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("spindle_dynamics", ["spindle_dynamics.pyx"])]

setup(name='kt_simul',
      version='0.2.a',
      description='Kinetochore dynamics simulation',
      author='Guillaume Gay',
      author_email='gay@cict.fr',
      url='http://www-lbcmcp.ups-tlse.fr',
      cmdclass = {'build_ext': build_ext},
      packages = ['kt_simul', 'kt_simul.gui'],
      package_dir = {'kt_simul':''},
      package_data = {'kt_simul': ['params.xml', 'measures.xml']},
      data_files = [('gui/images', ['gui/images/exec.svg',
                                    'gui/images/open.png',
                                    'gui/images/save.png',
                                    'gui/images/player_pause.svg',
                                    'gui/images/player_play.svg'])]
      )
