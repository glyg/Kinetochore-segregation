from __future__ import with_statement
from fabric.api import *
from fabric.colors import *

import os
import socket

# Local paths
project = os.path.abspath(os.path.dirname(__file__))
project = os.path.dirname(os.path.dirname(project))
if socket.gethostname() == "aragorn":
    results = "/media/thor/data/ktsimu/"
else:
    results = "/home/hadim/local/data/"

# Loki paths
host = 'hadim@130.120.107.234'
rproject = "/home/hadim/dev/kt_simul/"
rresults = "/home/hadim/local/data/ktsimu/"
rpython = "/home/hadim/local/virtualenvs/ktsimu/bin/python "
rpythonpath = "/home/hadim/dev/kt_simul"

env.hosts = [host, ]

NSIMU = 10000

@task
def push():
    """
    Push code to host
    """
    print red("Push kt_simul to loki")
    with lcd(project):
        local("rsync --progress -a --delete ../kt_simul/ %s:%s" % (host, rproject))


@task
def launch(simu = None):
    """
    Launch simulations
    """
    if not simu:
        simu = NSIMU
    push()
    with cd(os.path.join(rproject, "scripts")):
        run("workon ktsimu")
        cmd = rpython + "cluster.py --path %s --nsimu %s" % (rresults, str(simu))
        # Allow to run in background
        run("screen -dmS ktsimu " + cmd)


@task
def kill():
    """
    Ugly way : kill all python process
    """
    run("killall -9 python")

@task
def status(keep=False):
    """
    Display the status of all running simulations
    """
    with cd(rresults):
        for simu in run("ls").split(" "):
            if simu:
                simu_path = os.path.join(rresults, simu)
                files = run("ls %s" % simu_path).split(" ")
                # If no simu.log then simu is running
                if not "simu.log" in files:
                    _show_status(simu_path, keep)


def _show_status(path, keep):
    print red("Simulation %s" % path)
    if keep:
        last = " -f "
    else:
        last = " -1 "
    run("tail %s %s" % (last, os.path.join(path, "run.log")))
