# import packages
# opensees packages
import openseespy.opensees as ops
import opsvis as ovs
import opstool as otl
#other packages
import numpy as np
import matplotlib.pyplot as plt
from modelUnits import *

gap = 0*cm


Building10RNodesDisp = np.loadtxt('Building10RightNodes_Disp.txt', delimiter=" ")
Building20LNodesDisp = np.loadtxt('Building20LeftNodes_Disp.txt', delimiter=" ")
collison = np.zeros_like(Building10RNodesDisp)

poundingForceKelvinVoigt = np.loadtxt('testForceKelvinVoigt.txt', delimiter=" ")
poundingForceEPP = np.loadtxt('testForceSpringEPPGAP.txt', delimiter=" ")

NStorey = min(len(Building10RNodesDisp[0]) -1, len(Building20LNodesDisp[0])-1)
timeSeries = Building10RNodesDisp[:, 0]

fig, ax = plt.subplots(NStorey,1)

for storey in range(1, NStorey+1):
    collisonSeries = gap + Building10RNodesDisp[:, storey] - (Building20LNodesDisp[:, storey])
    #print(collisonSeries * 1000)
    isCollided = collisonSeries >= 0
    collisonSeries[~isCollided] = 0
    poundingForceKelvinVoigt[~isCollided] = 0
    poundingForceEPP[~isCollided] = 0

    ax[storey - 1].plot(timeSeries, poundingForceKelvinVoigt[:, storey])
    ax[storey - 1].plot(timeSeries, isCollided/1e13)
    #print(collisonSeries * 1000)
    #ax[storey - 1][1].plot(timeSeries, poundingForce[:, storey])
    #ax[storey - 1].plot(timeSeries, Building10RNodesDisp[:, storey], label='left Building')
    #ax[storey - 1].plot(timeSeries, Building20LNodesDisp[:, storey], label='right Building')


plt.show()


