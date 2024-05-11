# elastic perfectly plastic gap material

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

ops.fix(1, 1)

# epp GAP
eppGAPMatID = 4

E = 1425 * ksi
Fy = -15.46 * ksi
gap = -0.001
eta = -0.01
ops.uniaxialMaterial('ElasticPPGap', eppGAPMatID, E, Fy, gap, eta)

#element('zeroLength', eleTag, *eleNodes, '-mat', *matTags, '-dir', *dirs)
ops.element('twoNodeLink', 1, *[1, 2], '-mat', eppGAPMatID, '-dir', *[1])


# load assignment
linTS = 1
ops.timeSeries('Linear', linTS)
ops.pattern('Plain', 1, linTS)

#0.000106952 displacement occurs for 10000 N load
px = 1*kN
ops.load(2, *[-1*px])
#ops.sp(2, 1, 0.000106952)

# recorders
ops.recorder('Node', '-file', 'Pounding Models/eppDisp.txt', '-time', '-closeOnWrite','-node',2 , '-dof', 1, 'disp')
ops.recorder('Node', '-file', 'Pounding Models/eppReactions.txt', '-time', '-closeOnWrite','-node',1 , '-dof', 1, 'reaction')
ops.record()

# analysis


ops.constraints('Plain')
ops.numberer('RCM')
ops.test('NormDispIncr', 1.0e-6, 6, 0)
ops.algorithm('Newton')
ops.system('BandGeneral')  #system profileSPD can't solve for negative strain hardening ratio
ops.integrator('LoadControl', .1)  #load control can't push our model beyond yield point for negative strain hardening
#ops.integrator('DisplacementControl', 2, 1, 0.01)
ops.analysis('Static')

for i in range(0, 10):
    status = ops.analyze(1)
    
    if status != 0:
        print("Trying other analysis params.")
        ops.algorithm('ModifiedNewton')
        ops.test('NormDispIncr', 1.0e-3, 10, 0)
        status = ops.analyze(1)
    
    if status != 0:
        print("Breaking analysis")
        break
    #print(ops.getTime(), ops.nodeReaction(1), ops.nodeReaction(2))


disp = np.loadtxt("Pounding Models/eppDisp.txt")
rxn = np.loadtxt("Pounding Models/eppReactions.txt")
plt.plot(disp[:, 1], -rxn[:, 1], '.')
plt.show()

