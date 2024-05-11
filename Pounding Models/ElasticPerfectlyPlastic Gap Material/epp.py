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
ops.node(2, *[gap])

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


# sine series
sineWave = 1
#timeSeries('Trig', tag, tStart, tEnd, period, '-factor', factor=1.0, '-shift', shift=0.0, '-zeroShift', zeroShift=0.0)
ops.timeSeries('Trig', sineWave, 0, 40, 4, '-factor', 1.0) #apply sine series for 40s
ops.pattern('UniformExcitation', sineWave, 1, '-accel', 1)

# ground motion
eqLoad = 2
# ops.timeSeries('Path', eqLoad, '-dt', 0.01, '-filePath', "TakY.th", '-factor', g)
# ops.pattern('UniformExcitation', eqLoad, 1, '-accel', eqLoad)


# recorders
ops.recorder('Node', '-file', 'eppDisplacement.txt', '-time', '-closeOnWrite','-node',2 , '-dof', 1, 'disp')
ops.recorder('Node', '-file', 'eppcousReactions.txt', '-time', '-closeOnWrite','-node',1 , '-dof', 1, 'reaction')

# analysis


ops.constraints('Transformation')
ops.numberer('RCM')
ops.test('EnergyIncr', 1.0e-10, 100)
ops.algorithm('ModifiedNewton')
ops.system('UmfPack')
ops.integrator('Newmark', .5, .25)
#ops.integrator('DisplacementControl', 2, 1, 0.001)
ops.analysis('Transient')
ops.analyze(4000, 0.01)
# run analysis up to 4000 step, with time step = 0.01s, so total time = 40s

disp = np.loadtxt("eppDisplacement.txt")
rxn = np.loadtxt("eppReactions.txt")
plt.plot(disp[:, 1], -rxn[:, 1])
plt.show()

