# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Example 2: Surface reactions: IP3 receptor model 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import steps.interface

from steps.model import *
from steps.geom import *

###############################################################################

def getModel():

    # Create model container object
    mdl = Model()
    r = ReactionManager()
    with mdl:
        # Chemical species
        Ca, IP3 = Species.Create() # Calcium and IP3
        
        # Receptor state objects
        R     = Species.Create() # IP3 receptor in 'naive' state
        RIP3  = Species.Create() # bound IP3 
        Ropen = Species.Create() # bound IP3 and Ca (open)
        RCa   = Species.Create() # 1 bound Ca to inactivation site
        R2Ca  = Species.Create() # 2 bound Ca to inactivation sites
        R3Ca  = Species.Create() # 3 bound Ca to inactivation sites
        R4Ca  = Species.Create() # 4 bound Ca to inactivation sites
        
        surfsys = SurfaceSystem.Create()

        with surfsys:
            # The binding reactions:
            R.s + IP3.o   <r[1]> RIP3.s
            RIP3.s + Ca.o <r[2]> Ropen.s
            R.s + Ca.o    <r[3]> RCa.s
            RCa.s + Ca.o  <r[4]> R2Ca.s
            R2Ca.s + Ca.o <r[5]> R3Ca.s
            R3Ca.s + Ca.o <r[6]> R4Ca.s

            # The reaction constants
            r[1].K = 1000000000.0, 25800.0
            r[2].K = 8000000000.0, 2000.0
            r[3].K = 8889000.0, 5.0
            r[4].K = 20000000.0, 10.0
            r[5].K = 40000000.0, 15.0
            r[6].K = 60000000.0, 20.0

            # Ca ions passing through open IP3R channel
            Ca.i + Ropen.s <r[1]> Ropen.s + Ca.o
            # Corresponds to Ca input ~ 20000/ms for open receptor
            r[1].K = 8000000.0, 8000000.0

    # return model container object
    return mdl

###############################################################################

def getGeom(mdl):

    geom = Geometry()
    with geom:
        # Create cytosol compartment
        cyt = Compartment.Create()
        # Assign volume to cytosol 
        cyt.Vol = 1e-19
        
        # Create ER compartment
        ER = Compartment.Create()
        # Assign volume to ER 
        ER.Vol = 2e-20
        
        # Create the ER membrane and associate it with a surface system
        ERmemb = Patch.Create(ER, cyt, mdl.surfsys)

    # return geometry container object
    return geom

###############################################################################
