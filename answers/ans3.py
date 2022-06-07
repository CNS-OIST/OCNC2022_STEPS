# Uncomment the line below and run the cell to reload the empty exercise
#%load ./answers/blank3.py

from steps.sim import *
from steps.saving import *

# Create a simulation using the 'Wmdirect' solver
sim_wm = Simulation('Wmdirect', mdl, wmgeom, rng)

# Set-up automatic recording of the concentration of each molecule every 0.01 seconds
rs = ResultSelector(sim_wm)

counts_wm = rs.cyt.ALL(Species).Count

sim_wm.toSave(counts_wm, dt=0.01)

# Start a new run and set the initial conditions
sim_wm.newRun()

sim_wm.cyt.MEKp.Conc = 1e-6
sim_wm.cyt.ERK.Conc = 1.5e-6

# Run the simulation for 30 seconds
sim_wm.run(30)
