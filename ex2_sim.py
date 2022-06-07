# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Example 2: Surface reactions: IP3 receptor model 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import steps.interface

from steps.rng import *
from steps.sim import *
from steps.saving import *

from matplotlib import pyplot as plt
import numpy  as np

import ex2_ip3model as ip3r_model 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Simulation control variables
NITER = 100
T_END = 0.201
DT = 0.001

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import model
mdl = ip3r_model.getModel()

# Import geometry 
geom = ip3r_model.getGeom(mdl)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create random number generator
r = RNG('mt19937', 512, 654)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create well-mixed solver object
sim = Simulation('Wmdirect', mdl, geom, r)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

rs = ResultSelector(sim)

results = rs.cyt.Ca.Conc * 1e6 << rs.ERmemb.Ropen.Count

sim.toSave(results, dt=DT)

# Run the simulation
for i in range(NITER):
    
    # Start a new run
    sim.newRun()
    
    # Set initial conditions
    sim.cyt.Ca.Conc = 3.30657e-08
    sim.cyt.IP3.Conc = 2.5e-06
    sim.ER.Ca.Conc = 0.00015
    sim.ERmemb.R.Count = 16
    
    # Clamp ER calcium to fixed concentration
    sim.ER.Ca.Clamped = True

    # Run the simulation
    sim.run(T_END)

print(f'Ran {NITER} sim iterations')

# Numpy array manipulation
mean_res = np.mean(results.data, 0)
std_res = np.std(results.data, 0)
res_std1 = mean_res + std_res
res_std2 = mean_res - std_res

tpnts = results.time[0]

# Plot Ca concentration
plt.plot(tpnts, mean_res[:, 0], color = 'black', linewidth = 1.0, label = 'mean')
plt.plot(tpnts, res_std1[:, 0], color = 'gray', linewidth = 0.5, label = 'standard deviation')
plt.plot(tpnts, res_std2[:, 0], color = 'gray', linewidth = 0.5)
plt.xlabel('Time (sec)')
plt.ylabel('Cytosolic Ca concentration ($\mu$M)')
plt.title(f'{NITER} iterations')
plt.ylim(0)
plt.legend()
plt.show()

# Plot open IP3 receptors
plt.plot(tpnts, mean_res[:, 1], color = 'black', linewidth = 1.0, label = 'mean')
plt.plot(tpnts, res_std1[:, 1], color = 'gray', linewidth = 0.5, label = 'standard deviation')
plt.plot(tpnts, res_std2[:, 1], color = 'gray', linewidth = 0.5)
plt.xlabel('Time (sec)')
plt.ylabel('Number of open IP3 receptors')
plt.title(f'{NITER} iterations')
plt.ylim(0)
plt.legend()
plt.show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
