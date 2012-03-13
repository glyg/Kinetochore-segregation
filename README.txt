ABOUT:
------

This code provides the simulation of mitotic spindle elements (for now, the
kinetochore and the spindle pole bodies), during cell division in
fission yeast. 

The underlying model is fully described in: 

Gay et al. 'A stochastic model of kinetochore–microtubule attachment
accurately describes fission yeast chromosome segregation' 
J. Cell Biol 2012 doi: 10.1083/jcb.201107124

LICENCE:
--------

This code is provided under the GPL compatible CeCILL licence (see
LICENCE.txt for full details).


DEPENDENCIES:
-------------
I generally use my package manager versions of the python libraries. 

-Python >= 2.5
-Numpy >= 1.4 and Scipy >= 0.9 
-Cython >= 0.14
-Qt4 and PyQt4

INSTALLATION:
-------------

You should first clone the github version of this code, then
use the setup script, whether via:

# python setup.py build 

for a local install, or 

# python setup.py install 

for a system install.
You will need a C compiler for the cython part.


USAGE:
------

>>> import kt_simul.simul_spindle as sim
>>> metaph = sim.Metaphase()
>>> metaph.simul()
>>> metaph.show_one()
>>> help(metaph)

Should provide a first view of the simulation. 

Using the GUI:

# python [path_to_package]/kt_smul/gui/kineto_simulation.py 

Further details can be found by looking at the code 
Please send e-mail to gllm.gay-at-gmail.com for further details



