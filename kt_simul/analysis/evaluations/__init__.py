# -*- coding: utf-8 -*-
"""
Data processing, analysis, batch scripts
"""

import sys
import os

def find_evaluations(path, cls):
    """
    Source: http://www.luckydonkey.com/2008/01/02/python-style-plugins-made-easy/

    Find all subclass of cls in py files located below path
    (does look in sub directories)

    :param path: the path to the top level folder to walk
    :type path: str
    :param cls: the base class that all subclasses should inherit from
    :type cls: class
    :rtype: list
    :return: a list if classes that are subclasses of cls
    """

    subclasses=[]

    def look_for_subclass(modulename):
        module=__import__(modulename)
        print module
        # walk the dictionaries to get to the last one
        d=module.__dict__
        for m in modulename.split('.')[1:]:
            print m
            d=d[m].__dict__

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

    plugin_files = [x[:-3] for x in os.listdir(path) if x.endswith(".py")]
    sys.path.insert(0, path)

    for plugin in plugin_files:
        look_for_subclass("kt_simul.analysis.evaluations." + plugin)

    return subclasses

class Evaluation(object):
    def __init__(self,):
        pass

    def hello(self):
        print "salut grand"