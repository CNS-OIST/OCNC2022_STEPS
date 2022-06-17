# Uncomment the line below and run the cell to load the answer
#%load ./answers/ans1.py

import steps.interface

# Import biochemical model module
from steps.model import *

# Create model container
mdl = Model()

with mdl:
    # Create chemical species
    MEKp, ERK, MEKpERK, ERKp = Species.Create()

    # Create reaction set container (volume system)
    vsys = VolumeSystem.Create()

    with vsys:
        r = ReactionManager()

        # MEKp + ERK <-> MEKpERK, rate constants 16.2*10e6 (forward) and 0.6 (backward)

        # MEKpERK -> MEKp + ERKp, rate constant 0.15

# Check that everything is in the model
print(*mdl.ALL())
