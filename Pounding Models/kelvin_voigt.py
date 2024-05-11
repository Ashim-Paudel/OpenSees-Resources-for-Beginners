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

ops.wipe()
ops.model("BasicBuilder", '-ndm', 1, '-ndf', 1)

gap = 2*cm
ops.node(1, *[0])
ops.node(2, *[0])
ops.node(3, *[gap])
#ops.node(4, *[0, 0])

ops.mass(2, 100)
ops.mass(3, 100)

ops.fix(1, 1)
#ops.fix(4, 1, 1, 1)


#hertzDampID = 1
#K1 = 0.0
#K2 = 0.0
#delY = 0.0
#ops.uniaxialMaterial('ImpactMaterial', hertzDampID, K1, K2, delY, gap)

# viscous material
viscousID = 2
C = 683 
alpha = 1
ops.uniaxialMaterial('Viscous', viscousID, C, alpha)

# spring materil
springID = 3
Fy = 1000
E0 = 200
b = 0.1
ops.uniaxialMaterial('Steel01', springID, Fy, E0, b)

# epp GAP
eppGAPMatID = 4
E = 2* E0
Fy = 250*Mpa
ops.uniaxialMaterial('ElasticPPGap', eppGAPMatID, -1*E, -1*Fy, -1*gap, 0.1)


### kelvin voigt construction
parallelTag = 100
ops.uniaxialMaterial('Parallel', parallelTag, *[viscousID, springID])



#element('zeroLength', eleTag, *eleNodes, '-mat', *matTags, '-dir', *dirs)
ops.element('zeroLength', 1, *[1, 2], '-mat', springID, '-dir', *[1])
#ops.element('twoNodeLink', 2, *[2,3], '-mat', eppGAPMatID, '-dir', *[1])
#ops.element('zeroLength', 3, *[2,4], '-mat', springID, '-dir', *[1, 2, 6])

# load assignment
# load assignment
linTS = 1
#ops.timeSeries('Path', 1, '-dt', 0.01, '-filePath', "Pounding Models/TakY.th", '-factor', g)
#sine wave
#timeSeries('Trig', tag, tStart, tEnd, period, '-factor', factor=1.0, '-shift', shift=0.0, '-zeroShift', zeroShift=0.0)
ops.timeSeries('Trig', 1, 0, 5000, 2, '-factor', 1000.0)
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)

# recorders
ops.recorder('Node', '-file', 'kelvin_voigt_Disp.txt', '-time', '-closeOnWrite','-node',2 , '-dof', 1, 'disp')
ops.recorder('Node', '-file', 'kelvin_voigt_Reactions.txt', '-time', '-closeOnWrite','-node',1 , '-dof', 1, 'reaction')
ops.record()

# analysis

ops.constraints('Transformation')
ops.numberer('RCM')
ops.test('EnergyIncr', 1.0e-10, 100)
ops.algorithm('KrylovNewton')
ops.system('UmfPack')
ops.integrator('Newmark', .5, .25)
#ops.integrator('DisplacementControl', 2, 1, 0.001)
ops.analysis('Transient')
ops.analyze(50000, 0.001)

# ovs.plot_model()
# plt.show()
#ops.analyze(10)



disp = np.loadtxt("kelvin_voigt_Disp.txt")
rxn = np.loadtxt("kelvin_voigt_Reactions.txt")
plt.plot(disp[:, 1], -rxn[:, 1])
plt.show()

