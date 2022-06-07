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

        # Create diffusion rules for species SA, SB, and SC
        diff_SA = Diffusion.Create(SA, 0.02e-9)
        diff_SB = Diffusion.Create(SB, 0.02e-9)
        # First argument gives the affected species, second argument gives the diffusion coefficient in m^2/s
        diff_SC = Diffusion.Create(SC, 0.02e-9)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import geometry module
from steps.geom import *

# Import the mesh and scale it so that dimensions are in meters
mesh = TetMesh.LoadAbaqus('meshes/1x1x1_cube.inp', 1e-6)

with mesh:
    # Create the compartment by providing a list of all tetrahedrons
    cyt = Compartment.Create(mesh.tets, vsys)

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
sim_tet = Simulation('Tetexact', mdl, mesh, rng)

# Create a result selector associated to sim_tet
rs = ResultSelector(sim_tet)

# Specify which values should be saved, here all species counts in compartment "cyt"
counts_tet = rs.cyt.ALL(Species).Count

# Save these values every 0.001s
sim_tet.toSave(counts_tet, dt=0.001)

# Signalize the start of a new run
sim_tet.newRun()

# Initialize the number of 'SA' molecules to 10
sim_tet.cyt.SA.Count = 10

# Or you can set the concentration (M), as for 'SB'
sim_tet.cyt.SB.Conc = 3.32e-08

# Run simulation until 0.5s, data is recorded automatically
sim_tet.run(0.5)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import matplotlib via pylab
from matplotlib import pyplot as plt

# Use Matlpotlib functions to plot data from both simulations
# Note that we access data from run 0
plt.plot(counts_tet.time[0], counts_tet.data[0])

plt.ylabel('Number of molecules')
plt.xlabel('Time (sec)')
plt.legend(counts_tet.labels)
plt.show()

