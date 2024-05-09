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
#ops.node(3, *[gap])
#ops.node(4, *[0, 0])

ops.fix(1, 1)
#ops.fix(4, 1, 1, 1)


#hertzDampID = 1
#K1 = 0.0
#K2 = 0.0
#delY = 0.0
#ops.uniaxialMaterial('ImpactMaterial', hertzDampID, K1, K2, delY, gap)

# viscous material
viscousID = 2
C = 683 * kN*sec/m
alpha = 1
ops.uniaxialMaterial('Viscous', viscousID, C, alpha)

# spring materil
springID = 3
Fy = 250*Mpa
E0 = 93500 * kN/m**2
b = -0.1
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
linTS = 1
ops.timeSeries('Linear', linTS)
ops.pattern('Plain', 1, linTS)

#0.000106952 displacement occurs for 10000 N load
px = 1000*kN
ops.load(2, *[-px])
#ops.sp(2, 1, 0.000106952)


# recorders
ops.recorder('Node', '-file', 'Pounding Models/linear_springDisp.txt', '-time', '-closeOnWrite','-node',2 , '-dof', 1, 'disp')
ops.recorder('Node', '-file', 'Pounding Models/linear_springReactions.txt', '-time', '-closeOnWrite','-node',1 , '-dof', 1, 'reaction')
ops.record()

# analysis


ops.constraints('Plain')
ops.numberer('RCM')
ops.test('NormDispIncr', 1.0e-6, 6)
ops.algorithm('Newton')
ops.system('ProfileSPD')
ops.integrator('LoadControl', 1)
#ops.integrator('DisplacementControl', 2, 1, 0.026738)
ops.analysis('Static')

# ovs.plot_model()
# plt.show()
#ops.analyze(10)

for i in range(0, 100):
    status = ops.analyze(1)
    
    if status != 0:
        print("Trying other analysis params.")
        ops.algorithm('ModifiedNewton')
        ops.test('NormDispIncr', 1.0e-10, 10, 0)
        status = ops.analyze(1)
    
    if status != 0:
        print("Breaking analysis")
        break
    #print(ops.getTime(), ops.nodeReaction(1), ops.nodeReaction(2))


disp = np.loadtxt("Pounding Models/linear_springDisp.txt")
rxn = np.loadtxt("Pounding Models/linear_springReactions.txt")
plt.plot(disp[:, 1], -rxn[:, 1])
plt.show()

