# Translational spring # Static analysis
#

# Written: MHS

# Date: Feb 2000

# Units: kip, in

# Define the model builder

model BasicBuilder -ndm 2 -ndf 3

# Set some variables

set L 120

set A 20

set E 30000

set I 1400

# Define nodes

node 1 0.0 0.0

node 2 0.0 0.0

node 3 $L 0.0

# Define single point constraints

fix 1 1 1 1

fix 3 1 1 1

# ndR ndC dofs

equalDOF 1 2 1 3

# Define force-deformation relationship for spring

uniaxialMaterial ElasticPP 2 1050 0.02

uniaxialMaterial Elastic 3 -50

uniaxialMaterial Parallel 1 2 3

# id ndI ndJ mat dir

element zeroLength 1 1 2 -mat 1 -dir 2

# Geometric transformation

geomTransf Linear 1

# id ndI ndJ A E I transf

element elasticBeamColumn 2 2 3 $A $E $I 1

pattern Plain 1 Linear {

load 2 0.0 10.0 0.0

}

recorder Node ZeroLength1.out disp -time -node 2 -dof 2

integrator LoadControl 1 1 1 1

test EnergyIncr 1.0e-6 10 1

algorithm Newton

numberer Plain

constraints Transformation 1.0

system SparseGeneral -piv

analysis Static

analyze 5