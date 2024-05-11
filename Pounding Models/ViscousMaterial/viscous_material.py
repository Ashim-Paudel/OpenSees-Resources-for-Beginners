# A simple implementation of Viscous Material to record force deformation properties
# the material is subjected to uniform sinewave 


import openseespy.opensees as ops

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

ops.mass(2, 100) # 100 kg mass

ops.fix(1, *[1])

# viscous material
viscousID = 2
C = 683 * kN*sec/m
alpha = 1
#ops.uniaxialMaterial('Viscous', viscousID, C, alpha)
ops.uniaxialMaterial('Viscous', viscousID, C, alpha)

#element('zeroLength', eleTag, *eleNodes, '-mat', *matTags, '-dir', *dirs)
ops.element('zeroLength', 1, *[1, 2], '-mat', viscousID, '-dir', *[1, 2, 3])

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
ops.recorder('Node', '-file', 'viscousDisplacement.txt', '-time', '-closeOnWrite','-node',2 , '-dof', 1, 'disp')
ops.recorder('Node', '-file', 'viscousReactions.txt', '-time', '-closeOnWrite','-node',1 , '-dof', 1, 'reaction')

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

disp = np.loadtxt("viscousDisplacement.txt")
rxn = np.loadtxt("viscousReactions.txt")
plt.plot(disp[:, 1], -rxn[:, 1])
plt.show()
