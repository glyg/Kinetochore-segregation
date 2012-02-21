#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("spindle_dynamics", ["spindle_dynamics.pyx"])]

setup(name='kt_simul',
      version='1.0',
      description='Kinetochore dynamics simulation',
      author='Guillaume Gay',
      author_email='gllm.gay@gmail.com',
      url='http://www-lbcmcp.ups-tlse.fr',
      cmdclass = {'build_ext': build_ext},
      packages = ['kt_simul', 'kt_simul.gui'],
      package_dir = {'kt_simul':''},
      package_data = {'kt_simul': ['LICENCE.txt',
                                   'default/params.xml',
                                   'default/measures.xml']},
      data_files = [('gui/images', ['gui/images/exec.svg',
                                    'gui/images/open.png',
                                    'gui/images/save.png',
                                    'gui/images/player_pause.svg',
                                    'gui/images/player_play.svg'])]
      )
