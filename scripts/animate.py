import sys

from kt_simul.draw.animate import Animator
from kt_simul.io.simuio import SimuIO
from kt_simul.core.simul_spindle import Metaphase

if len(sys.argv) > 1 and sys.argv[1] == "--new":

    meta = Metaphase(paramfile="params.xml", verbose=True)
    meta.simul()
    io = SimuIO(meta)
    io.save("results.zip")

    anim = Animator(meta)
    anim.play()

else:

    meta = SimuIO().read("results.zip")

    anim = Animator(meta)
    anim.play()
