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

gap = 5*cm
ops.node(1, *[0])
ops.node(2, *[gap/2])

ops.fix(1, 1)

# epp GAP
eppGAPMatID = 4
E = 2e9
Fy = -250*Mpa
eta = 0.0

ops.uniaxialMaterial('ElasticPPGap', eppGAPMatID, E, Fy, -1*gap, eta, 'damage')
ops.element('twoNodeLink', 1, *[1, 2], '-mat', eppGAPMatID, '-dir', *[1])


# load assignment
linTS = 1
ops.timeSeries('Linear', linTS)
ops.pattern('Plain', 1, linTS)

px = -1000*kN
ops.load(2, *[px])


# recorders
ops.recorder('Element', '-file', 'eppDisplacement.txt', '-time', '-closeOnWrite','-ele',1 , '-dof', 1, 'deformation')
ops.recorder('Element', '-file', 'eppReactions.txt', '-time', '-closeOnWrite','-ele',1 , '-dof', 1, 'force')

ops.constraints('Plain')
ops.numberer('RCM')
ops.test('NormDispIncr', 1.0e-6, 6, 0)
ops.algorithm('ModifiedNewton')
ops.system('BandGeneral') 
ops.analysis('Static')

NumSteps = 40
ops.integrator('DisplacementControl', 2, 1, -0.01)
ops.analyze(NumSteps)
ops.integrator('DisplacementControl', 2, 1, 0.01)
ops.analyze(NumSteps)

NumSteps = 60
ops.integrator('DisplacementControl', 2, 1, -0.01)
ops.analyze(NumSteps)
ops.integrator('DisplacementControl', 2, 1, 0.01)
ops.analyze(NumSteps)

NumSteps = 80
ops.integrator('DisplacementControl', 2, 1, -0.01)
ops.analyze(NumSteps)
ops.integrator('DisplacementControl', 2, 1, 0.01)
ops.analyze(NumSteps)

#print(ops.getTime(), ops.nodeReaction(1), ops.nodeReaction(2))
disp = np.loadtxt("eppDisplacement.txt")
rxn = np.loadtxt("eppReactions.txt")
plt.plot(disp[:, 1]*1000, -rxn[:, 1]/100)
plt.xlabel("Strain (mm)")
plt.ylabel("Stress (MPa)")
plt.title("Stress-Strain Relationship of Elastic Perfectly Plastic Gap Material in Compression")
plt.grid()
plt.savefig("EPP_Gap_compression.png")
plt.show()
