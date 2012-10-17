from kt_simul.core import simul_spindle as sim
from kt_simul.draw import Animator

PARAMFILE = "params.xml"
meta = sim.Metaphase(verbose=True)

anim = Animator(meta)
anim.play()
