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
    #element('nonlinearBeamColumn', eleTag, *eleNodes, numIntgrPts, secTag, transfTag)
    numItgrPts = 5
    colId = 1 #elements starting with this for col
    beamId = 2
    # section tags
    ColSecTag = 1
    BeamSecTag = 2
    for j in range(NStory):
        for i in range(NBay + 1):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            eleColTag = int(f"{colId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            inodeTag = int(f"{j+1}{i+1}")
            jnodeTag = int(f"{j+2}{i+1}")
            print(eleColTag, inodeTag, jnodeTag)

            ops.element('nonlinearBeamColumn', eleColTag, *[inodeTag, jnodeTag], numItgrPts, ColSecTag, )

    for j in range(NStory + 1):
        for i in range(NBay):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            eleColTag = int(f"{beamId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            inodeTag = int(f"{j+1}{i+1}")
            jnodeTag = int(f"{j+1}{i+2}")
            print(eleColTag, inodeTag, jnodeTag)

            #ops.element('nonlinearBeamColumn', eleColTag, *[inodeTag, jnodeTag], numItgrPts, ColSecTag)
            #print(nodeTag, *nodeCoor)
            #ops.node(nodeTag, *(nodeCoor))

def getElasticSection():
    # section tags
    ColSecTag = 1
    BeamSecTag = 2

    # section geometry
    # column sections: W27x114
    AgCol = 33.5*sqinch	                # cross-sectional area
    IzCol = 4090.*inch4	                # moment of Inertia

    # beam sections: W24x94
    AgBeam = 27.7*sqinch		        # cross-sectional area
    IzBeam = 2700.*inch4	            # moment of Inertia

    # material properties
    Es = 29000*ksi		# Steel Young's Modulus
    nu = 0.3
    Gs = Es/2./(1+nu)  # Torsional stiffness Modulus

    # define geometric transformation tag
    ColSecTransf = 1
    BeamSecTransf = 2
    ops.geomTransf('Linear', ColSecTransf)
    ops.geomTransf('Linear', BeamSecTransf)

    ops.section('Elastic', BeamSecTag, Es, AgBeam, IzBeam)
    ops.section('Elastic', ColSecTag, Es, AgCol, IzCol)

getModel(NBay=3, NStory=3, LBeam=2, LCol=2)