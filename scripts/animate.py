from kt_simul.draw.animate import Animator
from kt_simul.io import SimuIO
from kt_simul.core.simul_spindle import Metaphase

simu = False

if simu:
    meta = Metaphase(verbose=True)
    meta.simul()
    io = SimuIO(meta)
    io.save("results.xml", "data.npy")
else:
    meta = SimuIO().read("results.xml")

anim = Animator(meta)
anim.play()
