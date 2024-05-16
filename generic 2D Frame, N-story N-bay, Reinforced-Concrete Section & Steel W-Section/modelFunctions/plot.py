# import packages
# opensees packages
import openseespy.opensees as ops
import opsvis as ovs
import opstool as otl
#other packages
import numpy as np
import matplotlib.pyplot as plt


Building10RNodesDisp = np.loadtxt('Building10RightNodes_Disp.txt', delimiter=" ")
Building20LNodesDisp = np.loadtxt('Building20LeftNodes_Disp.txt', delimiter=" ")
collison = np.zeros_like(Building10RNodesDisp)
print(collison)
poundingForce = np.loadtxt('testForce.txt', delimiter=" ")


timeSeries = Building10RNodesDisp[:, 0]
print(len(Building10RNodesDisp))

fig, ax = plt.subplots(5, 1)

for storey in range(1, 6):
    collisonSeries = Building10RNodesDisp[:, storey] - Building20LNodesDisp[:, storey]
    isCollided = collisonSeries < 0
    collisonSeries[~isCollided] = 0

    ax[storey - 1].plot(timeSeries, poundingForce[:, storey]*1e28)
    ax[storey - 1].plot(timeSeries, collisonSeries)
    
    #ax[storey - 1][1].plot(timeSeries, poundingForce[:, storey])
    #ax[storey - 1].plot(timeSeries, Building10RNodesDisp[:, storey], label='left Building')
    #ax[storey - 1].plot(timeSeries, Building20LNodesDisp[:, storey], label='right Building')

plt.legend()
plt.show()


