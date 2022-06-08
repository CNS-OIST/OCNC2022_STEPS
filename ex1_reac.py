# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Example 1: Second-order reaction, well-mixed simulation 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import the new steps interface
import steps.interface

# Import biochemical model module
from steps.model import *

# Create model container
mdl = Model()

# Create reaction manager
r = ReactionManager()

# Using this model, declare species and reactions
with mdl:
    # Create chemical species
    SA, SB, SC  = Species.Create()

    # Create reaction set container
    vsys = VolumeSystem.Create()

    with vsys:
        # Create reaction
        # SA + SB - > SC with rate 200 /uM.s
        SA + SB >r[1]> SC
        r[1].K = 200e6

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import geometry module
from steps.geom import *

# Create well-mixed geometry container
wmgeom = Geometry()

# Using this geometry, declare compartments
with wmgeom:
    # Create cytosol compartment, associate it with vsys, give volume in m^3
    cyt = Compartment.Create(vsys, vol=1e-18)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import random number generator module
from steps.rng import *

# Create random number generator, with buffer size as 256 and starting seed 899
rng = RNG('mt19937', 256, 899)

# Could use time to get random seed
# import time
# rng = RNG('mt19937', 256, int(time.time()))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import simulation and data saving modules
from steps.sim import *
from steps.saving import *

# Create well-mixed stochastic simulation object
sim_wm = Simulation('Wmdirect', mdl, wmgeom, rng)

# Create a result selector associated to sim_wm
rs = ResultSelector(sim_wm)

# Specify which values should be saved, here all species counts in compartment "cyt"
counts_wm = rs.cyt.ALL(Species).Count

# Save these values every 0.001s
sim_wm.toSave(counts_wm, dt=0.001)

# Signalize the start of a new run
sim_wm.newRun()

# Initialize the number of 'SA' molecules to 10
sim_wm.cyt.SA.Count = 10

# Or you can set the concentration (M), as for 'SB'
sim_wm.cyt.SB.Conc = 3.32e-08

# Run simulation until 0.5s, data is recorded automatically
sim_wm.run(0.5)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create well-mixed deterministic solver
sim_rk4 = Simulation('Wmrk4', mdl, wmgeom, rng)

# Set the integration time-step (s)
sim_rk4.setRk4DT(0.00001)

# Declare values to be saved for the deterministic solver
rs = ResultSelector(sim_rk4)
counts_rk4 = rs.cyt.ALL(Species).Count
sim_rk4.toSave(counts_rk4, dt=0.001)

# Repeat the simulation process for the deterministic solver
sim_rk4.newRun()
sim_rk4.cyt.SA.Count = 10
sim_rk4.cyt.SB.Conc = 3.32e-08
sim_rk4.run(0.5)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Import matplotlib
from matplotlib import pyplot as plt

# Use Matlpotlib functions to plot data from both simulations
# Note that we access data from run 0
plt.plot(counts_wm.time[0], counts_wm.data[0])
plt.plot(counts_rk4.time[0], counts_rk4.data[0], color='black')

plt.ylabel('Number of molecules')
plt.xlabel('Time (sec)')
plt.legend(counts_wm.labels)
plt.show()
