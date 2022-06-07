# Uncomment the line below and run the cell to reload the empty exercise
#%load ./answers/blank2.py

from steps.geom import *
from steps.rng import *

# Create a well-mixed geometry container
wmgeom = Geometry()

# Create a compartment of  0.1um^3 and associate it to volume system 'vsys'
with wmgeom:
    cyt = Compartment.Create(vsys, 0.1e-18)

# Create a 'r123' random number generator
rng = RNG('r123', 256, 1)
