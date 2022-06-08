# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Example 3: Unbounded diffusion 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

from matplotlib import pyplot as plt
import numpy as np

########################################################################

# The Abaqus file to import
MESHFILE = 'meshes/sp_10r_3875'
# The number of iterations to run 
NITER = 10
# The data collection time increment (s)
DT = 0.001
# The simulation endtime (s)
INT = 0.1

# The number of molecules to be injected into the centre
NINJECT = 1000

# The number of tetrahedral elements to sample data from. 
SAMPLE = 500

# The diffusion constant for our diffusing species (m^2/s)
DCST= 20.0e-12

########################################################################

def gen_model():         
    mdl = Model()
    with mdl:
        SA = Species.Create()
        vsys = VolumeSystem.Create()
        with vsys:
            Diffusion(SA, DCST)
    return mdl

########################################################################

def gen_geom(meshPath, mdl):

    print("Loading mesh...")
    mesh = TetMesh.LoadAbaqus(meshPath, 1e-6)
    print("Mesh Loaded")

    with mesh:
        # Create a compartment containing all tetrahedron
        cyto = Compartment.Create(mesh.tets, mdl.vsys)

        print("Finding tetrahedron samples...")
        # List to hold tetrahedrons
        tets = TetList()

        # Fetch the central tetrahedron
        ctet = mesh.tets[0, 0, 0]
        tets.append(ctet)
        tets += ctet.neighbs

        # Find the maximum and minimum coordinates of the mesh
        bmin = mesh.bbox.min
        bmax = mesh.bbox.max

        # Run a loop until we have stored all tet indices we require
        while len(tets) < SAMPLE:
            # Pick a random position within the bounding box
            pos = bmin + (bmax - bmin) * np.random.random(3)
            if pos in mesh.tets:
                tets.append(mesh.tets[pos])

        # Find the radial distance of the tetrahedrons to mesh center:
        tetrads = np.array([np.linalg.norm(tet.center - ctet.center)*1e6 for tet in tets])

        print("Tetrahedron samples found")

    return mesh, tets, tetrads

########################################################################

model = gen_model()
mesh, tets, tetrads = gen_geom(MESHFILE, model)

rng = RNG('mt19937', 512, 2903)

sim = Simulation('Tetexact', model, mesh, rng)

rs = ResultSelector(sim)

AConc = rs.TETS(tets).SA.Conc

sim.toSave(AConc, dt=DT)

for i in range(NITER):
    sim.newRun()
    print("Running iteration", i)
    # Inject all molecules into the central tet:
    sim.TET(0, 0, 0).SA.Count = NINJECT
    sim.run(INT)

########################################################################

def plotres(tidx, tpnts, res_mean, tetrads):
    if tidx >= len(tpnts):
        raise ValueError('Time index is out of range.')
    t = tpnts[tidx]

    plt.scatter(tetrads, res_mean[tidx,:], s=2)
    plt.xlabel('Radial distance of tetrahedron ($\mu$m)')            
    plt.ylabel('Concentration in tetrahedron ($\mu$M)')
    plt.title(f'Unbounded diffusion. Time: {t}s')

    plotAnalytical(t)

    plt.xlim(0.0, max(tetrads))
    plt.ylim(0.0)
    plt.show()

########################################################################

def plotAnalytical(t):     
    segs = 100     
    anlytconc = np.zeros((segs))     
    radialds = np.zeros((segs))     
    maxrad = 0.0     
    for i in tetrads:         
        if (i > maxrad): maxrad = i     
    maxrad *= 1e-6     
    intervals = maxrad/segs     
    rad = 0.0     
    for i in range((segs)):         
        # Find the conc from analytical solution, and convert to mol/L         
        anlytconc[i]=1.0e3*(1/6.022e23)* \
                ((NINJECT/((4*np.pi*DCST*t) ** 1.5))* \
                (np.exp((-1.0*(rad*rad))/(4*DCST*t))))         
        radialds[i] = rad*1e6         
        rad += intervals     
    plt.plot(radialds, anlytconc, color = 'red')

########################################################################

tpnts = AConc.time[0]
res_mean = np.mean(AConc.data, axis=0)*1e6

plotres(10, tpnts, res_mean, tetrads)
plotres(100, tpnts, res_mean, tetrads)
