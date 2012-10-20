# -*- coding: utf-8 -*-
"""
Pool evaluations are small plugin that can be automatically
launch after severals simulations.

PoolEvaluation take a path where .zip files are stored. To open simulations
uses kt_simul.simuio.SimuIO().read()
"""

import sys
import os

from kt_simul.utils.filesystem import tree

EVAL_PATH = os.path.abspath(os.path.dirname(__file__))

def find_pool_evaluations(groups=[]):
    """
    This function return a list of PoolEvaluation classes that are enabled. In the future, it will have the capability to filter evaluations.

    .. note:: Source from http://www.luckydonkey.com/2008/01/02/python-style-plugins-made-easy/

    .. todo:: add selector and filter

    :param groups: Select plugins from these groups
    :type groups: str

    :rtype: list
    :return: a list if classes that are subclasses of cls
    """

    path = EVAL_PATH
    cls = PoolEvaluation

    subclasses=[]

    def look_for_subclass(modulename):
        module=__import__(modulename)

        # walk the dictionaries to get to the last one
        d=module.__dict__
        for m in modulename.split('.')[1:]:
            try:
                d=d[m].__dict__
            except:
                pass

        #look through this dictionary for things
        #that are subclass of Job
        #but are not Job itself
        for key, entry in d.items():
            if key == cls.__name__:
                continue
            try:
                if issubclass(entry, cls):
                    subclasses.append(entry)
            except TypeError:
                #this happens when a non-type is passed in to issubclass. We
                #don't care as it can't be a subclass of Job if it isn't a
                #type
                continue

    sys.path.insert(0, path)

    plugin_files = []
    for x in tree(path):
        if x.endswith(".py"):
            mod = x[:-3].replace(os.sep, ".")
            plugin_files.append(mod)

    sys.path.insert(0, path)

    for plugin in plugin_files:
        look_for_subclass("kt_simul.analysis.pool_evaluations." + plugin)

    # Filter plugin here
    evaluations = []
    for cls in subclasses:
        if cls.enable:
            if not groups and cls.group == None:
                evaluations.append(cls)
            else:
                if cls.group in groups:
                    evaluations.append(cls)

    return evaluations

class RunFunctionNotImplemented(Exception):
    pass

class PoolEvaluation(object):
    """
    Abstract class. All pool evaluation plugins need to write a
    class that inherit from PoolEvaluation.

    .. warning:: run() method need to be implemented.
    """

    def __init__(self,):
        pass

    def run(self, *args):
        raise RunFunctionNotImplemented("run() method need to be implemented in your Evaluation plugin.")

__all__ = []
for ev in find_pool_evaluations():
    __all__.append(ev)

__all__ += ['find_pool_evaluations', 'PoolEvaluation', 'RunFunctionNotImplemented']
