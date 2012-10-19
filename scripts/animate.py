from kt_simul.draw.animate import Animator
from kt_simul.io import SimuIO
from kt_simul.core.simul_spindle import Metaphase

simu = True

meta = Metaphase(verbose=True)
meta.simul()
io = SimuIO(meta)
io.save("results.xml", "data.npy")

anim = Animator(meta)
anim.play()

########################

meta = SimuIO().read("results.xml")

anim = Animator(meta)
anim.play()
