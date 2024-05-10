import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# Translational spring #Static analysis
#
# Written: MHS
# Date: Feb 2000
# Units: kip, in
# Define the model builder
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)
# Set some variables
# Define nodes
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 0.0)
ops.node(3, 120, 0.0)

# Define single point constraints
ops.fix(1, 1, 1, 1)
ops.fix(3, 1, 1, 1)


spring_orientation = 0

# ndR ndC dofs
if spring_orientation == 0:
    ops.equalDOF(1, 2, *[1, 3])
elif spring_orientation == 45:
    ops.equalDOF(1, 2, *[3])


# Define force-deformation relationship for spring
ops.uniaxialMaterial('ElasticPP', 2, 1050, 0.02)
ops.uniaxialMaterial('Elastic', 3, -50)
ops.uniaxialMaterial('Parallel', 1, 2, 3)

#element('zeroLength', eleTag, *eleNodes, '-mat', *matTags, '-dir', *dirs)
# id ndI ndJ mat dir
if spring_orientation == 0:
    ops.element('zeroLength', 1, *[1, 2], '-mat', 1, '-dir', 2)
elif spring_orientation == 45:
    ops.element('zeroLength', 1, *[1, 2], '-mat', 1, '-dir', 1, '-orient', 1, 1, 0, -1, 1, 0)
#


# Geometric transformation
ops.geomTransf('Linear', 1)

# id ndI ndJ A E I transf
ops.element('elasticBeamColumn', 2, *[2, 3], 20, 30000, 1400, 1)


ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 0.0, 10.0, 0.0)

ops.recorder('Node', '-file', 'translationalSpring_disp.out', '-time','-closeOnWrite',  '-node', 2, '-dof', 2, 'disp')
ops.recorder('Node', '-file', 'translationalSpring_reaction.out', '-time','-closeOnWrite',  '-node', 1,'-dof',1, 2,3, 'reaction')

ops.numberer('Plain')
ops.constraints('Transformation', 1.0)
ops.system('SparseGeneral', '-piv')
ops.integrator('LoadControl', 1)
ops.test('EnergyIncr', 1e-06, 10, 1)
ops.algorithm('Newton')
ops.analysis('Static')

with open("translationalSpring_appliedForce.txt", 'w') as file:
    for i in range(10):
        ops.analyze(1)
        #print(ops.nodeUnbalance(2))
        file.write(str(ops.nodeUnbalance(2)[1])+'\n')

disp22 = np.loadtxt('translationalSpring_disp.out')[:, 1]
rxn12 = np.loadtxt('translationalSpring_reaction.out')[:, 2]
force_applied = np.loadtxt('translationalSpring_appliedForce.txt')

fig, ax = plt.subplots()
ax.plot(disp22, -rxn12, label="Reaction on Node 1 dof 2")
ax.plot(disp22, force_applied, label="Force applied on Node 2 dof 2")
ax.set_title('Plot of Displacement, Reaction and Force Applied')
ax.set_xlabel("displacement (inch)")
ax.set_ylabel("force (kip)")
ax.set_yticks(force_applied, force_applied)
#ax[0].plot(rxn12, force_applied, label="rxn vs applied force")
#ax[0].set_title('Plot of Reaction and Force Applied')
plt.legend()
plt.show()