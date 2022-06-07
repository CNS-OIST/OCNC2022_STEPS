# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# Example 4: Hodgkin-Huxley Action Potential propagation model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # IMPORTS # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

import numpy as np
import os
import math

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # PARAMETERS  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # #

# Potassium conductance = 0.036 S/cm2
# Sodium conductance = 0.120 S/cm2

# Potassium single-channel conductance
K_G = 20.0e-12 # Siemens

# Potassium channel density
K_ro = 18.0e12 # per square meter

# Potassium reversal potential
K_rev = -77e-3 # volts

# Sodium single-channel conductance
Na_G = 20.0e-12 # Siemens

# Sodium channel density
Na_ro = 60.0e12 # per square meter

# Sodium reversal potential
Na_rev = 50e-3 # volts

# Leak single-channel conductance
L_G = 0.3e-12 # Siemens

# Leak density
L_ro = 10.0e12 # per square meter

# Leak reveral potential
leak_rev = -54.4e-3 # volts


# A table of potassium channel population factors: 
# n0, n1, n2, n3, n4
K_facs = [ 0.21768, 0.40513, 0.28093, 0.08647, 0.00979 ]

# A table of sodium channel population factors
# m0h0, m1h0, m2h0, m3h0, m0h1, m1h1, m2h1, m3h1:
Na_facs = [[0.34412, 0.05733, 0.00327, 6.0e-05],
           [0.50558, 0.08504, 0.00449, 0.00010]]

# # # # # # # # # # # # # # # # # # RATE FUNCTION # # # # # # # # # # # # # # # #

def HHRateFunction(A, B, C, D, F, H, V, abs_tol=1e-13):
    num = A + B * V
    denom = C + H * math.exp((V + D) / F)
    if math.isclose(num, 0, abs_tol=abs_tol) and math.isclose(denom, 0, abs_tol=abs_tol):
        return F * B / (H * math.exp((V + D) / F))
    else:
        return num / denom

# # # # # # # # # # # # # # # # # # MESH  # # # # # # # # # # # # # # # # # # # # 

meshfile_ab = 'meshes/axon_cube_L1000um_D443nm_equiv0.5_19087tets.inp'

# # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # # # #

# Temperature for gating kinetics
celsius = 20.0

# Current injection
Iclamp = 50.0e-12 #	amps

# Voltage range for gating kinetics in Volts
Vrange = [-100.0e-3, 50e-3, 1e-4]

# The simulation dt
DT_sim = 1.0e-4 # seconds

# The time until which the simulation should be run
ENDT = 4.0e-3

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # BIOCHEMICAL MODEL # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

model = Model()

r = ReactionManager()

with model:
    ssys = SurfaceSystem.Create()

    #  Potassium channel
    Ko, Kc = SubUnitState.Create()
    KSU = SubUnit.Create([Ko, Kc])
    VGKC = Channel.Create([KSU]*4)

    # Sodium channel
    Na_mo, Na_mc, Na_hi, Na_ha = SubUnitState.Create()
    NamSU, NahSU = SubUnit.Create(
        [Na_mo, Na_mc],
        [Na_hi, Na_ha]
    )
    VGNaC = Channel.Create([NamSU, NamSU, NamSU, NahSU])

    # Leak channel
    lsus = SubUnitState.Create()
    Leak = Channel.Create([lsus])

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    # Hodgkin-Huxley gating kinetics

    # Temperature dependence
    thi = math.pow(3.0, ((celsius-6.3)/10.0))

    _a_n = VDepRate(lambda V: thi * 1e3 * HHRateFunction(-0.55, -0.01, -1, 55, -10, 1, V*1e3), vrange=Vrange)
    _b_n = VDepRate(lambda V: thi * 1e3 * HHRateFunction(1, 0, 0, 65, 80, 8, V*1e3), vrange=Vrange)
    
    _a_m = VDepRate(lambda V: thi * 1e3 * HHRateFunction(-4, -0.1, -1, 40, -10, 1, V*1e3), vrange=Vrange)
    _b_m = VDepRate(lambda V: thi * 1e3 * HHRateFunction(1, 0, 0, 65, 18, 0.25, V*1e3), vrange=Vrange)

    _a_h = VDepRate(lambda V: thi * 1e3 * HHRateFunction(1, 0, 0, 65, 20, 1 / 0.07, V*1e3), vrange=Vrange)
    _b_h = VDepRate(lambda V: thi * 1e3 * HHRateFunction(1, 0, 1, 35, -10, 1, V*1e3), vrange=Vrange)

    with ssys:

        with VGKC[...]:
            Kc.s <r[1]> Ko.s
            r[1].K = _a_n, _b_n

        with VGNaC[...]:
            Na_hi.s <r[1]> Na_ha.s
            r[1].K = _a_h, _b_h
            
            Na_mc.s <r[1]> Na_mo.s
            r[1].K = _a_m, _b_m

        # Create ohmic current objects
        VGKC_I = OhmicCurr.Create(VGKC[Ko, Ko, Ko, Ko], K_G, K_rev)
        VGNaC_I = OhmicCurr.Create(VGNaC[Na_mo, Na_mo, Na_mo, Na_ha], Na_G, Na_rev)
        Leak_I = OhmicCurr.Create(Leak[lsus], L_G, leak_rev)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # TETRAHEDRAL MESH  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mesh = TetMesh.LoadAbaqus(meshfile_ab, scale=1e-6)

# # # # # # # # # # # # # # # MESH MANIPULATION # # # # # # # # # # # # # # # # #

facetris = TriList([tri for tri in mesh.tris if tri.center.z == mesh.bbox.min.z])
injverts = facetris.verts

print("Found ", len(injverts), "I_inject vertices")
print("Found ", len(facetris), "triangles on bottom face")

memb_tris = mesh.surface - facetris

# The points along (z) axis at which to record potential
pot_pos = np.arange(mesh.bbox.min.z, mesh.bbox.max.z, 10e-6)
pot_tet = [mesh.tets[(0, 0, z)] for z in pot_pos]

# # # # # # # # # # # # # # # GEOMETRY OBJECTS  # # # # # # # # # # # # # # # # #

with mesh:
    # Create cytosol compartment
    cyto = Compartment.Create(mesh.tets)

    # Create the patch and associate with surface system ssys
    patch = Patch.Create(memb_tris, cyto, None, ssys)

    # Create the membrane across which the potential will be solved
    membrane = Membrane.Create([patch], opt_method = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, 1234)

sim = Simulation('Tetexact', model, mesh, rng, True)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rs = ResultSelector(sim)

NaCurrs = rs.TRIS(memb_tris).VGNaC_I.I
KCurrs = rs.TRIS(memb_tris).VGKC_I.I
CellPot = rs.TETS(pot_tet).V

NaCurrs.metaData['trizpos'] = [tri.center.z for tri in memb_tris]
KCurrs.metaData['trizpos'] = [tri.center.z for tri in memb_tris]
NaCurrs.metaData['triarea'] = [tri.Area for tri in memb_tris]
KCurrs.metaData['triarea'] = [tri.Area for tri in memb_tris]
CellPot.metaData['tetzpos'] = pot_pos

sim.toSave(NaCurrs, KCurrs, CellPot, dt=DT_sim)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sim.newRun()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Inject channels
surfarea = sim.patch.Area

for state in VGNaC:
    prop = Na_facs[state.Count(Na_ha)][state.Count(Na_mo)]
    sim.patch.VGNaC[state].Count = Na_ro * surfarea * prop

for state in VGKC:
    prop = K_facs[state.Count(Ko)]
    sim.patch.VGKC[state].Count = K_ro * surfarea * prop

sim.patch.Leak[lsus].Count = L_ro * surfarea

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set some simulation variables:

# Set dt for membrane potential calculation to 0.01ms
sim.EfieldDT = 1.0e-5

# Initialize potential to -65mV
sim.membrane.Potential = -65e-3

# Set capacitance of the membrane to 1 uF/cm^2 = 0.01 F/m^2
sim.membrane.Capac = 1.0e-2

# Set resistivity of the conduction volume to 100 ohm.cm = 1 ohm.meter
sim.membrane.VolRes = 1.0

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Set the current clamp

sim.VERTS(injverts).IClamp = Iclamp/len(injverts)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Run the simulation
sim.run(ENDT)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from matplotlib import pyplot as plt

NBINS = 100

def getTIdx(saver, t):
    return min(enumerate(saver.time[0]), key=lambda v: abs(v[1] - t))[0]

def plotPotential(t):
    tidx = getTIdx(CellPot, t)
    plt.plot(
        CellPot.metaData['tetzpos'] * 1e6, 
        CellPot.data[0, tidx, :] * 1e3, 
        label=f'{CellPot.time[0, tidx]*1e3} ms'
    )

def plotCurrents(t):
    tidx = getTIdx(NaCurrs, t)
    for results, currName in zip([NaCurrs, KCurrs], ['Na', 'K']):
        data = results.data[0, tidx, :] * 1e12
        pos = results.metaData['trizpos'] * 1e6
        areas = results.metaData['triarea'] * 1e12
        bins = np.histogram_bin_edges(pos, NBINS)
        dig = np.digitize(pos, bins)
        # Ignore empty bins
        with np.errstate(invalid='ignore'):
            meanData = np.bincount(dig, weights=data) / np.bincount(dig, weights=areas)
            meanPos  = np.bincount(dig, weights=pos) / np.bincount(dig)
        plt.plot(meanPos, meanData, label=f'{currName} {results.time[0, tidx]*1e3} ms')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
plotPotential(1e-3)
plotPotential(2e-3)
plotPotential(3e-3)
plt.xlabel('Z-axis (um)')
plt.ylabel('Membrane potential (mV)')
plt.legend()
plt.show()
    
plotCurrents(1e-3)
plotCurrents(2e-3)
plotCurrents(3e-3)
plt.xlabel('Z-axis (um)')
plt.ylabel('Current  (pA/um^2)')
plt.legend()
plt.show()












# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # # IMPORTS # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# from __future__ import print_function
# import steps.model as smodel
# import steps.geom as sgeom
# import steps.rng as srng
# import steps.solver as ssolver
# import steps.utilities.meshio as meshio

# import numpy
# import math
# import time
# from random import *
# from pylab import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # # PARAMETERS  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # #

# # Potassium conductance = 0.036 S/cm2
# # Sodium conductance = 0.120 S/cm2

# # Potassium single-channel conductance
# K_G = 20.0e-12 # Siemens

# # Potassium channel density
# K_ro = 18.0e12 # per square meter

# # Potassium reversal potential
# K_rev = -77e-3 # volts

# # Sodium single-channel conductance
# Na_G = 20.0e-12 # Siemens

# # Sodium channel density
# Na_ro = 60.0e12 # per square meter

# # Sodium reversal potential
# Na_rev = 50e-3 # volts

# # Leak single-channel conductance
# L_G = 0.3e-12 # Siemens

# # Leak density
# L_ro = 10.0e12 # per square meter

# # Leak reveral potential
# leak_rev = -54.4e-3 # volts


# # A table of potassium channel population factors: 
# # n0, n1, n2, n3, n4
# K_facs = [ 0.21768, 0.40513, 0.28093, 0.08647, 0.00979 ]

# # A table of sodium channel population factors
# # m0h0, m1h0, m2h0, m3h0, m0h1, m1h1, m2h1, m3h1:
# Na_facs = [ 0.34412, 0.05733, 0.00327, 6.0e-05, \
                # 0.50558, 0.08504, 0.00449, 0.00010 ]

# # # # # # # # # # # # # # # # # # # MESH  # # # # # # # # # # # # # # # # # # # # 

# meshfile_ab = 'axon_cube_L1000um_D443nm_equiv0.5_19087tets.inp'

# # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # # # #

# # Temperature for gating kinetics
# celsius = 20.0		

# # Current injection
# Iclamp = 50.0e-12 #	amps

# # Voltage range for gating kinetics in Volts
# Vrange = [-100.0e-3, 50e-3, 1e-4]

# # The number of simulation time-points
# N_timepoints = 41

# # The simulation dt
# DT_sim = 1.0e-4 # seconds

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # BIOCHEMICAL MODEL # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# mdl = smodel.Model()
# ssys = smodel.Surfsys('ssys', mdl)

# # Potassium channel
# K = smodel.Chan('K', mdl)
# K_n0 = smodel.ChanState('K_n0', mdl, K)		
# K_n1 = smodel.ChanState('K_n1', mdl, K)
# K_n2 = smodel.ChanState('K_n2', mdl, K)
# K_n3 = smodel.ChanState('K_n3', mdl, K)
# K_n4 = smodel.ChanState('K_n4', mdl, K)

# # Sodium channel
# Na = smodel.Chan('Na', mdl)
# Na_m0h0 = smodel.ChanState('Na_m0h0', mdl, Na)
# Na_m1h0 = smodel.ChanState('Na_m1h0', mdl, Na)
# Na_m2h0 = smodel.ChanState('Na_m2h0', mdl, Na)
# Na_m3h0 = smodel.ChanState('Na_m3h0', mdl, Na)
# Na_m0h1 = smodel.ChanState('Na_m0h1', mdl, Na)
# Na_m1h1 = smodel.ChanState('Na_m1h1', mdl, Na)
# Na_m2h1 = smodel.ChanState('Na_m2h1', mdl, Na)
# Na_m3h1 = smodel.ChanState('Na_m3h1', mdl, Na)

# # Leak channel
# L = smodel.Chan('L', mdl)
# Leak = smodel.ChanState('Leak', mdl, L)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # Hodgkin-Huxley gating kinetics

# # Temperature dependence
# thi = math.pow(3.0, ((celsius-6.3)/10.0))

# _a_n = lambda mV: thi*((0.01*(10-(mV+65.))/(math.exp((10-(mV+65.))/10.)-1)))

# _b_n = lambda mV: thi*((0.125*math.exp(-(mV+65.)/80.)))

# _a_m = lambda mV: thi*((0.1*(25-(mV+65.))/(math.exp((25-(mV+65.))/10.)-1)))

# _b_m = lambda mV: thi*((4.*math.exp(-(mV+65.)/18.)))


# _a_h = lambda mV: thi*((0.07*math.exp(-(mV+65.)/20.)))

# _b_h = lambda mV: thi*((1./(math.exp((30-(mV+65.))/10.)+1)))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # Potassium channel voltage-dependent gating dynamics

# Kn0n1 = smodel.VDepSReac('Kn0n1', ssys, slhs = [K_n0], srhs = [K_n1], \
                            # k=lambda V: 1.0e3 *4.*_a_n(V*1.0e3), vrange = Vrange)
# Kn1n2 = smodel.VDepSReac('Kn1n2', ssys, slhs = [K_n1], srhs = [K_n2], \
                            # k=lambda V: 1.0e3 *3.*_a_n(V*1.0e3), vrange = Vrange)
# Kn2n3 = smodel.VDepSReac('Kn2n3', ssys, slhs = [K_n2], srhs = [K_n3], \
                            # k=lambda V: 1.0e3 *2.*_a_n(V*1.0e3), vrange = Vrange)
# Kn3n4 = smodel.VDepSReac('Kn3n4', ssys, slhs = [K_n3], srhs = [K_n4], \
                            # k=lambda V: 1.0e3 *1.*_a_n(V*1.0e3), vrange = Vrange)

# Kn4n3 = smodel.VDepSReac('Kn4n3', ssys, slhs = [K_n4], srhs = [K_n3], \
                            # k=lambda V: 1.0e3 *4.*_b_n(V*1.0e3), vrange = Vrange)
# Kn3n2 = smodel.VDepSReac('Kn3n2', ssys, slhs = [K_n3], srhs = [K_n2], \
                            # k=lambda V: 1.0e3 *3.*_b_n(V*1.0e3), vrange = Vrange)
# Kn2n1 = smodel.VDepSReac('Kn2n1', ssys, slhs = [K_n2], srhs = [K_n1], \
                            # k=lambda V: 1.0e3 *2.*_b_n(V*1.0e3), vrange = Vrange)
# Kn1n0 = smodel.VDepSReac('Kn1n0', ssys, slhs = [K_n1], srhs = [K_n0], \
                            # k=lambda V: 1.0e3 *1.*_b_n(V*1.0e3), vrange = Vrange)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # Sodium channel voltage-dependent gating dynamics

# Na_m0h1_m1h1 = smodel.VDepSReac('Na_m0h1_m1h1', ssys, \
                                # slhs=[Na_m0h1], srhs=[Na_m1h1], \
                                # k=lambda V:1.0e3*3.*_a_m(V*1.0e3), vrange=Vrange)
# Na_m1h1_m2h1 = smodel.VDepSReac('Na_m1h1_m2h1', ssys, \
                                # slhs=[Na_m1h1], srhs=[Na_m2h1], \
                                # k=lambda V:1.0e3*2.*_a_m(V*1.0e3), vrange=Vrange)
# Na_m2h1_m3h1 = smodel.VDepSReac('Na_m2h1_m3h1', ssys, \
                                # slhs=[Na_m2h1], srhs=[Na_m3h1], \
                                # k=lambda V:1.0e3*1.*_a_m(V*1.0e3), vrange=Vrange)

# Na_m3h1_m2h1 = smodel.VDepSReac('Na_m3h1_m2h1', ssys, \
                                # slhs=[Na_m3h1], srhs=[Na_m2h1], \
                                # k=lambda V:1.0e3*3.*_b_m(V*1.0e3), vrange=Vrange)
# Na_m2h1_m1h1 = smodel.VDepSReac('Na_m2h1_m1h1', ssys, \
                                # slhs=[Na_m2h1], srhs=[Na_m1h1], \
                                # k=lambda V:1.0e3*2.*_b_m(V*1.0e3), vrange=Vrange)
# Na_m1h1_m0h1 = smodel.VDepSReac('Na_m1h1_m0h1', ssys, \
                                # slhs=[Na_m1h1], srhs=[Na_m0h1], \
                                # k=lambda V:1.0e3*1.*_b_m(V*1.0e3), vrange=Vrange)

# Na_m0h0_m1h0 = smodel.VDepSReac('Na_m0h0_m1h0', ssys, \
                                # slhs=[Na_m0h0], srhs=[Na_m1h0], \
                                # k=lambda V:1.0e3*3.*_a_m(V*1.0e3), vrange=Vrange)
# Na_m1h0_m2h0 = smodel.VDepSReac('Na_m1h0_m2h0', ssys, \
                                # slhs=[Na_m1h0], srhs=[Na_m2h0], \
                                # k=lambda V:1.0e3*2.*_a_m(V*1.0e3), vrange=Vrange)
# Na_m2h0_m3h0 = smodel.VDepSReac('Na_m2h0_m3h0', ssys, \
                                # slhs=[Na_m2h0], srhs=[Na_m3h0], \
                                # k=lambda V:1.0e3*1.*_a_m(V*1.0e3), vrange=Vrange)

# Na_m3h0_m2h0 = smodel.VDepSReac('Na_m3h0_m2h0', ssys, \
                                # slhs=[Na_m3h0], srhs=[Na_m2h0], \
                                # k=lambda V:1.0e3*3.*_b_m(V*1.0e3), vrange=Vrange) 
# Na_m2h0_m1h0 = smodel.VDepSReac('Na_m2h0_m1h0', ssys, \
                                # slhs=[Na_m2h0], srhs=[Na_m1h0], \
                                # k=lambda V:1.0e3*2.*_b_m(V*1.0e3), vrange=Vrange)
# Na_m1h0_m0h0 = smodel.VDepSReac('Na_m1h0_m0h0', ssys, \
                                # slhs=[Na_m1h0], srhs=[Na_m0h0], \
                                # k=lambda V:1.0e3*1.*_b_m(V*1.0e3), vrange=Vrange)

# Na_m0h0_m0h1 = smodel.VDepSReac('Na_m0h0_m0h1', ssys, \
                                # slhs=[Na_m0h0], srhs=[Na_m0h1], \
                                # k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)
# Na_m1h0_m1h1 = smodel.VDepSReac('Na_m1h0_m1h1', ssys, \
                                # slhs=[Na_m1h0], srhs=[Na_m1h1], \
                                # k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)
# Na_m2h0_m2h1 = smodel.VDepSReac('Na_m2h0_m2h1', ssys, \
                                # slhs=[Na_m2h0], srhs=[Na_m2h1], \
                                # k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)
# Na_m3h0_m3h1 = smodel.VDepSReac('Na_m3h0_m3h1', ssys, \
                                # slhs=[Na_m3h0], srhs=[Na_m3h1], \
                                # k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)

# Na_m0h1_m0h0 = smodel.VDepSReac('Na_m0h1_m0h0', ssys, \
                                # slhs=[Na_m0h1], srhs=[Na_m0h0], \
                                # k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)
# Na_m1h1_m1h0 = smodel.VDepSReac('Na_m1h1_m1h0', ssys, \
                                # slhs=[Na_m1h1], srhs=[Na_m1h0], \
                                # k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)
# Na_m2h1_m2h0 = smodel.VDepSReac('Na_m2h1_m2h0', ssys, \
                                # slhs=[Na_m2h1], srhs=[Na_m2h0], \
                                # k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)
# Na_m3h1_m3h0 = smodel.VDepSReac('Na_m3h1_m3h0', ssys, \
                                # slhs=[Na_m3h1], srhs=[Na_m3h0], \
                                # k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # Create ohmic current objects

# OC_K = smodel.OhmicCurr('OC_K', ssys, chanstate=K_n4, g=K_G, erev=K_rev)	
# OC_Na = smodel.OhmicCurr('OC_Na', ssys, chanstate=Na_m3h1, g=Na_G, erev=Na_rev)
# OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate=Leak, g=L_G, erev=leak_rev) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # TETRAHEDRAL MESH  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# mesh = meshio.importAbaqus('meshes/'+meshfile_ab, 1e-6)[0]

# # # # # # # # # # # # # # # # MESH MANIPULATION # # # # # # # # # # # # # # # # #

# # Find the vertices for the current clamp and store in a list
# injverts = []
# for i in range(mesh.nverts):
	# if ((mesh.getVertex(i)[2] < (mesh.getBoundMin()[2]+0.1e-6))):
		# injverts.append(i)
# print("Found ", injverts.__len__(), "I_inject vertices")

# facetris = []
# for i in range(mesh.ntris):
	# tri = mesh.getTri(i) 
	# if ((tri[0] in injverts) and (tri[1] in injverts) and (tri[2] in injverts)):
		# facetris.append(i)
# print("Found ", facetris.__len__(), "triangles on bottom face")

# memb_tris = list(mesh.getSurfTris()) 

# # Remove triangles on bottom face from membrane triangles
# for t in facetris: memb_tris.remove(t)


# # Bin the surface triangles for recording current
# bins_n = 100
# memb_tris_binned = [None]*bins_n
# mtb_area = numpy.zeros(bins_n)

# # In m
# bin_dz = 1000.0e-6/bins_n

# # The centre positions of the bins
# bin_pos = numpy.arange((bin_dz/2.0), 1000e-6, bin_dz)

# for m in range(bins_n): memb_tris_binned[m]=[]

# # Bin the triangles 
# for t in memb_tris:
    # barycz = mesh.getTriBarycenter(t)[2]
    # idx = 0
    # for p in bin_pos:
        # if (barycz >= p-(bin_dz/2.0) and barycz < p+(bin_dz/2.0)): 
            # memb_tris_binned[idx].append(t)
            # mtb_area[idx]+=(mesh.getTriArea(t)*1.0e12)
            # break
        # idx +=1

# # The points along (z) axis at which to record potential
# pot_pos = numpy.arange(mesh.getBoundMin()[2], mesh.getBoundMax()[2], 10e-6)
# pot_n = len(pot_pos)

# pot_tet = numpy.zeros(pot_n, dtype = 'uint')

# i=0
# for p in pot_pos:
    # # Axis is aligned with z-axis
    # pot_tet[i] = mesh.findTetByPoint([0.0, 0.0, pot_pos[i]])
    # i=i+1

# # # # # # # # # # # # # # # # GEOMETRY OBJECTS  # # # # # # # # # # # # # # # # #

# # Create cytosol compartment
# cyto = sgeom.TmComp('cyto', mesh, range(mesh.ntets))

# # Create the patch and associate with surface system 'ssys'
# patch = sgeom.TmPatch('patch', mesh, memb_tris, cyto)
# patch.addSurfsys('ssys')

# # Create the membrane across which the potential will be solved
# membrane = sgeom.Memb('membrane', mesh, [patch], opt_method = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # Create the random number generator
# r = srng.create('mt19937',512)
# r.initialize(int(time.time()%10000))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Create solver object: 'True' argument tells solver to activate EField object
# sim = ssolver.Tetexact(mdl, mesh, r, True)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Inject channels

# surfarea = sim.getPatchArea('patch')

# sim.setPatchCount('patch', 'Na_m0h0', Na_ro*surfarea*Na_facs[0])	 
# sim.setPatchCount('patch', 'Na_m1h0', Na_ro*surfarea*Na_facs[1])	 
# sim.setPatchCount('patch', 'Na_m2h0', Na_ro*surfarea*Na_facs[2])	 
# sim.setPatchCount('patch', 'Na_m3h0', Na_ro*surfarea*Na_facs[3])	 
# sim.setPatchCount('patch', 'Na_m0h1', Na_ro*surfarea*Na_facs[4])	 
# sim.setPatchCount('patch', 'Na_m1h1', Na_ro*surfarea*Na_facs[5])	 
# sim.setPatchCount('patch', 'Na_m2h1', Na_ro*surfarea*Na_facs[6])	 
# sim.setPatchCount('patch', 'Na_m3h1', Na_ro*surfarea*Na_facs[7])

# sim.setPatchCount('patch', 'K_n0', K_ro*surfarea*K_facs[0])
# sim.setPatchCount('patch', 'K_n1', K_ro*surfarea*K_facs[1])			
# sim.setPatchCount('patch', 'K_n2', K_ro*surfarea*K_facs[2])			
# sim.setPatchCount('patch', 'K_n3', K_ro*surfarea*K_facs[3])
# sim.setPatchCount('patch', 'K_n4', K_ro*surfarea*K_facs[4])

# sim.setPatchCount('patch', 'Leak', L_ro * surfarea)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Set some simulation variables:

# # Set dt for membrane potential calculation to 0.01ms
# sim.setEfieldDT(1.0e-5)

# # Initialize potential to -65mV
# sim.setMembPotential('membrane', -65e-3)

# # Set capacitance of the membrane to 1 uF/cm^2 = 0.01 F/m^2
# sim.setMembCapac('membrane', 1.0e-2)

# # Set resistivity of the conduction volume to 100 ohm.cm = 1 ohm.meter
# sim.setMembVolRes('membrane', 1.0)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # Set the current clamp
# niverts = injverts.__len__()
# for t in injverts:
    # sim.setVertIClamp(t, Iclamp/niverts) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Create result structures
# res = numpy.zeros((N_timepoints, pot_n))
# res_I_Na = numpy.zeros((N_timepoints, bins_n))
# res_I_K = numpy.zeros((N_timepoints, bins_n))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Run the simulation
# for l in range(N_timepoints):
    # print("\nTpnt: ", l,)
    
    # sim.run(DT_sim*l)
    
    # for b in range(bins_n):
        # # Record currents across membrane triangles
        # for mt in memb_tris_binned[b]:
            # res_I_Na[l,b]+= sim.getTriOhmicI(mt, 'OC_Na')*1.0e12
            # res_I_K[l,b]+= sim.getTriOhmicI(mt, 'OC_K')*1.0e12

        # res_I_Na[l,b]/=mtb_area[b]
        # res_I_K[l,b]/=mtb_area[b]
    
    # # Record voltage along axis
    # for p in range(pot_n):
        # res[l,p] = sim.getTetV(int(pot_tet[p]))*1.0e3

# results = (res,  pot_pos, res_I_Na, res_I_K, bin_pos)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# tpnt = arange(0.0, N_timepoints*DT_sim, DT_sim)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Function to plot voltage along the axis
# def plotVz(tidx):
    # if (tidx >= tpnt.size): 
        # print('Time index out of range')
        # return
    # plot(results[1]*1e6, results[0][tidx], \
         # label=str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    # legend(numpoints=1)
    # xlim(0, 1000)
    # ylim(-80,40)
    # xlabel('Z-axis (um)')
    # ylabel('Membrane potential (mV)')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# # Function to plot membrane currents along the axis
# def plotIz(tidx, plotstyles = ['-', '--']):
    # if (tidx >= tpnt.size): 
        # print('Time index out of range')
        # return
    # plot(results[4]*1e6, results[2][tidx], plotstyles[0],\
         # label = 'Na: '+str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    # plot(results[4]*1e6, results[3][tidx], plotstyles[1],\
         # label = 'K: '+str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    # legend(loc='best')
    # xlim(0, 1000)
    # ylim(-10, 15)
    # xlabel('Z-axis (um)')
    # ylabel('Current  (pA/um^2)')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# # END
