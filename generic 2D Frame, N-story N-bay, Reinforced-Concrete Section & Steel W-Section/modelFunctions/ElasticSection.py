# import packages
# opensees packages
import openseespy.opensees as ops
import opsvis as ovs
import opstool as otl
#other packages
import numpy as np
import matplotlib.pyplot as plt

# user def modules
from modelUnits import *

# Sign Convention
# All column related tags will proceed with 1
# All beam related tags will proceed with 2

# Global Variables for model
# Section tags for Elastic Beam and Column Section
ColSecTag = 1
BeamSecTag = 2

# Additional tags for Inelastic Beam Column Section
ColMatTagAxial = 101                              # Represent column axial behaviour
ColMatTagFlex = 102                               # Represent column flexural behaviour
BeamMatTagAxial = 201                             # Represent beam axial behaviour
BeamMatTagFlex = 202                              # Represent beam flexural behaviour 

# Geometric Transformation Tags
ColSecTransf = 1
BeamSecTransf = 2
ColSecTransfType = 'Linear'

# model building function
def getModel(NBay, NStory, LBeam, LCol, sectionType = 'Elastic'):
    ops.wipe()

    ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)

    for j in range(NStory + 1):
        for i in range(NBay + 1):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            nodeTag = int(f"{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            nodeCoor = (i * LBeam, j * LCol)
            
            print(nodeTag, *nodeCoor)
            ops.node(nodeTag, *(nodeCoor))

    # fixity to the nodes
    ops.fixY(0.0, 1, 1, 1)

    # determine support nodes where ground motions are input, for multiple-support excitation
    level = 1
    iSupportNode = [int(f"{level}{i+1}") for i in range(0, NBay + 1 )]

    # Set up parameters that are particular to the model for displacement control
    IDctrlNode =  (NStory+1)*10+1		# node where displacement is read for displacement control
    IDctrlDOF =  1	                    # DoF of displacement read for displacement control
    LBuilding =  NStory*LCol	        # total building height

    # building up the elements
    if sectionType == 'Elastic':
        getElasticSection()

    for nod in iSupportNode:
        print(ops.nodeCoord(nod))

    # Building the elements
    numItgrPts = 5
    colNameId = 1                               #elements starting with this for col
    beamNameId = 2

    for j in range(NStory):
        for i in range(NBay + 1):
            # nodeTag = j * (NBay+1) + i + 1    # conventional node tag
            eleColTag = int(f"{colNameId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            inodeTag = int(f"{j+1}{i+1}")
            jnodeTag = int(f"{j+2}{i+1}")
            print(eleColTag, inodeTag, jnodeTag, ops.nodeCoord(inodeTag), ops.nodeCoord(jnodeTag))

            # element('nonlinearBeamColumn', eleTag, *eleNodes, numIntgrPts, secTag, transfTag)
            ops.element('nonlinearBeamColumn', eleColTag, *[inodeTag, jnodeTag], numItgrPts, ColSecTag, ColSecTransf)

    for j in range(1, NStory + 1): #because we dont have beam at ground level
        for i in range(NBay):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            eleBeamTag = int(f"{beamNameId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            inodeTag = int(f"{j+1}{i+1}")
            jnodeTag = int(f"{j+1}{i+2}")
            print(eleBeamTag, inodeTag, jnodeTag, ops.nodeCoord(inodeTag), ops.nodeCoord(jnodeTag))

            ops.element('nonlinearBeamColumn', eleBeamTag, *[inodeTag, jnodeTag], numItgrPts, BeamSecTag, BeamSecTransf)


def getElasticSection():
    # section geometry
    # column sections: W27x114
    AgCol = 33.5*sqinch	                # cross-sectional area
    IzCol = 4090.*inch4	                # moment of Inertia
    # beam sections: W24x94
    AgBeam = 27.7*sqinch		        # cross-sectional area
    IzBeam = 2700.*inch4	            # moment of Inertia

    # material properties
    Es = 29000*ksi		# Steel Young's Modulus

    ops.geomTransf(ColSecTransfType, ColSecTransf)
    ops.geomTransf('Linear', BeamSecTransf)

    ops.section('Elastic', BeamSecTag, Es, AgBeam, IzBeam)
    ops.section('Elastic', ColSecTag, Es, AgCol, IzCol)


def getInelasticSection():
    # MATERIAL properties 
    Fy =  6.0*ksi       # Yield Stress 
    Es =  29000*ksi		# Steel Young's Modulus
    nu =  0.3           # Poisson's Ratio
    Gs =  Es/(2.*(1+nu)) 	# Torsional stiffness Modulus
    
    # SECTION PROPERTIES
    # COLUMN section W27x114
    AgCol = 33.5*pow(inch,2)		# cross-sectional area
    IzCol = 4090.*pow(inch,4)		# moment of Inertia
    EICol = Es*IzCol				# EI, for moment-curvature relationship
    EACol = Es*AgCol				# EA, for axial-force-strain relationship
    MyCol = 2e4*kip*inch	   		# yield moment
    PhiYCol = 0.25e-3/inch	   		# yield curvature
    PhiYCol = MyCol/EICol			# yield curvature
    EIColCrack = MyCol/PhiYCol		# cracked section inertia
    b = 0.01  					    # strain-hardening ratio 

    ops.uniaxialMaterial('Steel01', ColMatTagFlex, MyCol, EIColCrack, b)
    ops.uniaxialMaterial('Elastic', ColMatTagAxial, EACol)
    ops.section('Aggregator', ColSecTag, *[ColMatTagAxial, 'P', ColMatTagFlex, 'Mz'])

    # BEAM SECTION W24x94
    AgBeam = 27.7*pow(inch,2)		# cross-sectional area
    IzBeam = 2700.*pow(inch,4)		# moment of Inertia
    EIBeam = Es*IzBeam				# EI, for moment-curvature relationship
    EABeam = Es*AgBeam				# EA, for axial-force-strain relationship
    MyBeam = 1.5e4*kip*inch	   		# yield moment
    PhiYBeam = 0.25e-3/inch	   		# yield curvature
    PhiYBeam = MyBeam/EIBeam			# yield curvature
    EIBeamCrack = MyBeam/PhiYBeam		# cracked section inertia
    b = 0.01  					    # strain-hardening ratio 

    ops.uniaxialMaterial('Steel01', BeamMatTagFlex, MyBeam, EIBeamCrack, b)  # bilinear behavior for flexure
    ops.uniaxialMaterial('Elastic', BeamMatTagAxial, EABeam)
    ops.section('Aggregator', BeamSecTag, *[BeamMatTagAxial, 'P', BeamMatTagFlex, 'Mz'])

getModel(NBay=3, NStory=3, LBeam=3, LCol=3)
