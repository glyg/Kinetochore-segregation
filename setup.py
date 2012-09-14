#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("spindle_dynamics", ["spindle_dynamics.pyx"])]

setup(name='kt_simul',
      version='0.6',
      description='Kinetochore dynamics simulation',
      author='Guillaume Gay',
      author_email='gllm.gay@gmail.com',
      url='https://github.com/glyg/Kinetochore-segregation',
      cmdclass = {'build_ext': build_ext},
      packages = ['kt_simul', 'kt_simul.gui', 'kt_simul.core', 'kt_simul.analysis'],
      package_dir = {'kt_simul':''},
      package_data = {'kt_simul': ['LICENCE.txt',
                                   'default/params.xml',
                                   'default/measures.xml',
                                   'core/components.pxd',
                                   'core/components.pyx',
                                   'core/spindle_components.pyx',
                                   'core/spindle_dynamics.pyx',]},
      data_files = [('gui/images', ['gui/images/exec.svg',
                                    'gui/images/open.png',
                                    'gui/images/save.png',
                                    'gui/images/player_pause.svg',
                                    'gui/images/player_play.svg'])]
      )
