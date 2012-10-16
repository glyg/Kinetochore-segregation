from kt_simul.core import simul_spindle as sim
from kt_simul.analysis import processing

meta = sim.Metaphase(verbose = True)
meta.simul()

proc = processing.Processor(meta)
proc.show_trajs(fname = "trajs.pdf")
proc.show_one(fname = "one.png")

proc.write_results()
