# Uncomment the line below and run the cell to reload the empty exercise
#%load ./answers/blank5.py

# Create a new Model
mdl_tet = Model()

with mdl_tet:
    MEKp, ERK, MEKpERK, ERKp = Species.Create()
    vsys = VolumeSystem.Create()
    with vsys:
        r = ReactionManager()
        # Previously declared reactions
        MEKp + ERK <r[1]> MEKpERK
        r[1].K = 16.2*10e6, 0.6
        MEKpERK >r[1]> MEKp + ERKp
        r[1].K = 0.15

        # Add diffusion rules
        diff_MEKp = Diffusion.Create(MEKp, 30e-12)
        diff_ERK = Diffusion.Create(ERK, 30e-12)
        diff_MEKpERK = Diffusion.Create(MEKpERK, 10e-12)

# Load the 'meshes/spine.msh' mesh and create a cyt compartment
# The mesh is in micrometers and should be rescaled to meters
mesh = TetMesh.LoadGmsh('meshes/spine.msh', 1e-6)

with mesh:
    cyt = Compartment.Create(mesh.tets, vsys)

# Create a new simulation that uses the 'Tetexact' solver 
sim_tet = Simulation('Tetexact', mdl_tet, mesh, rng)

# Set-up automatic recording of the concentration of each molecule every 0.01 seconds
rs = ResultSelector(sim_tet)
counts_tet = rs.cyt.ALL(Species).Count
sim_tet.toSave(counts_tet, dt=0.01)

# Start a new run and set the initial conditions
sim_tet.newRun()

sim_tet.cyt.MEKp.Conc = 1e-6
sim_tet.cyt.ERK.Conc = 1.5e-6

# Run the simulation for 30 seconds
sim_tet.run(30)

# Plot results
plt.plot(counts_tet.time[0], counts_tet.data[0])
plt.ylabel('Number of molecules')
plt.xlabel('Time [s]')
plt.legend(counts_tet.labels)
plt.show()
