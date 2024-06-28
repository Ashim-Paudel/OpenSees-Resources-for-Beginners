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
import buildFiberSection

# Sign Convention
# All column related tags will proceed with 1
# All beam related tags will proceed with 2

# Global Variables for model
# Section tags for Elastic Beam and Column Section
ColSecTag = 1
BeamSecTag = 2

# Additional tags for Inelastic Beam Column Section
ColMatTagAxial = 101                              # Represent column axial behaviour
ColMatTagFlex = 102                               # Represent column flexural behaviour
BeamMatTagAxial = 201                             # Represent beam axial behaviour
BeamMatTagFlex = 202                              # Represent beam flexural behaviour 

# Geometric Transformation Tags
ColSecTransf = 1
BeamSecTransf = 2
ColSecTransfType = 'PDelta'


# model building function
def getModel(buildingID, NBay, NStory, LBeam, LCol, startCoor, **kwargs):
    """
    # Build a building model with parameters given
    ## Arguments
    buildingID = start building id with 100 [100, 200, 300, 400] etc
    NBay = no. of bays
    NStory = No. of storeys
    LBeam = Beam Length
    LCol = Column Length
    sectionType = any one of ['Elastic', 'InElastic', 'RCFiber', 'SteelFiber']
    startCoor = Coordinate of bottom left node of building default (0,0)
    """

    ### GLOBAL VARIALBLES ###
    #print(kwargs['HCol'])
    HCol =  kwargs.get('HCol')
    BCol =  kwargs.get('BCol')
    HBeam = kwargs.get('HBeam')
    BBeam = kwargs.get('BBeam')

    #ops.wipe()
    #ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)

    #we will create list to store all the extreme left and right nodes of two buildings
    LNodes = []
    RNodes = []

    for j in range(NStory + 1):
        for i in range(NBay + 1):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            nodeTag = int(f"{buildingID}{j+1}{i+1}") #example: nodeTag 10011 means 1st building(100) 1st node of ground floor (1st floor) 
            nodeCoor = (startCoor[0] + i * LBeam, startCoor[1] + j * LCol)
            if i == 0:
                LNodes.append(nodeTag)

            if (i+1 == NBay+1):
                RNodes.append(nodeTag)

            print(nodeTag, *nodeCoor)
            ops.node(nodeTag, *(nodeCoor))
        #ops.equalDOF(LNodes[j], RNodes[j], 2,3)


    exec(f"LNodes{buildingID} = {LNodes}", globals())
    exec(f"RNodes{buildingID} = {RNodes}", globals())

    # fixity to the nodes
    ops.fixY(0.0, 1, 1, 1)

    # determine support nodes where ground motions are input, for multiple-support excitation
    level = 1
    iSupportNode = [int(f"{buildingID}{level}{i+1}") for i in range(0, NBay + 1 )]

    # up parameters that are particular to the model for displacement control
    IDctrlNode =  (NStory+1)*10+1		# node where displacement is read for displacement control
    IDctrlDOF =  1	                    # DoF of displacement read for displacement control
    LBuilding =  NStory*LCol	        # total building height


    getRCFiberSection(buildingID, **kwargs)

    #### GRAVITY LOADS, MASSES AND WEIGHT ###
    ## unit weights ##
    GammaConcrete = 25*kN/m**3   		    # Reinforced-Concrete floor slabs
    GammaMasonry = 20*kN/m**3 
    FloorFin = 1*kN/m**2
    LiveLoad = 2*kN/m**2

    Twall = 230*mm                      # masonry wall thickness
    Qwall = GammaMasonry*Twall*LCol     # masonry wall weight per unit length of beam
    Tslab = 150*mm			            # 150 mm slab
    Lslab = 2*LBeam/2 			        # assume slab extends a distance of LBeam1/2 in/out of plane
    Qslab = (GammaConcrete + FloorFin + LiveLoad) * Tslab*Lslab   # slab dead weight
    QBeam = GammaConcrete*HBeam*BBeam

    QdlCol = GammaConcrete*HCol*BCol	# self weight of Column, weight per length
    QdlBeam = Qwall + Qslab + QBeam # dead load distributed along beam.
    WeightCol = QdlCol*LCol  		    # total single Column weight
    WeightBeam = QdlBeam*LBeam 	        # total single Beam weight

    ### Assigning masses to node and calculating building weight ###
    iFloorWeight = []                   # to store weight of each floor
    WeightTotal = 0.0
    sumWiHi = 0.0                       # sum of storey weight and its height for lateral load distribution
    ## uppermost node has weight from beam in the storey and 1/2column connecting it only
    ## outermost nodes has weight from half beam connecting it only
    for j in range(1, NStory + 1):
        FloorWeight = 0.0

        if (j+1) == (NStory + 1):      # Uppermost storey
            ColWeightFact = 1          # load shared from only one column
        if (j+1) < (NStory +1 ):
            ColWeightFact = 2          # load shared from 2 columns connecting it

        for i in range(NBay + 1):
            if (i+1) == 1 or (i+1) == (NBay+1):   # Outermost nodes 
                BeamWeightFact = 1                # Load shared from one beam only
            else:
                BeamWeightFact = 2

            nodeTag = int(f"{buildingID}{j+1}{i+1}") #example: nodeTag 10011 means 1st building(100) 1st node of ground floor (1st floor) 
            nodeCoor = ops.nodeCoord(nodeTag)
            
            WeightNode = ColWeightFact * WeightCol/2 + BeamWeightFact * WeightBeam/2
            MassNode = WeightNode / g
            ops.mass(nodeTag, *[MassNode, 0, 0])

            FloorWeight += WeightNode                 # FloorWeight =  Sum of all weights concn in nodes

            print(nodeTag, nodeCoor, MassNode)

        iFloorWeight.append(FloorWeight)
        WeightTotal += FloorWeight
        sumWiHi += FloorWeight * j * LCol
        
    MassTotal = WeightTotal / g                      # Total Building Mass

    for nod in iSupportNode:
        print(ops.nodeCoord(nod))
    
    # Building the Beam and Column Elements
    numItgrPts = 5
    colNameId = 1                               #elements starting with this for col
    beamNameId = 2
    # Column elements
    for j in range(NStory):
        for i in range(NBay + 1):
            # nodeTag = j * (NBay+1) + i + 1    # conventional node tag
            eleColTag = int(f"{buildingID}{colNameId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            inodeTag = int(f"{buildingID}{j+1}{i+1}")
            jnodeTag = int(f"{buildingID}{j+2}{i+1}")
            print(eleColTag, inodeTag, jnodeTag, ops.nodeCoord(inodeTag), ops.nodeCoord(jnodeTag))

            # element('nonlinearBeamColumn', eleTag, *eleNodes, numIntgrPts, secTag, transfTag)
            ops.element('nonlinearBeamColumn', eleColTag, *[inodeTag, jnodeTag], numItgrPts, ColSecTag, ColSecTransf)
    
    # Beam elements
    for j in range(1, NStory + 1): #because we dont have beam at ground level
        for i in range(NBay):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            eleBeamTag = int(f"{buildingID}{beamNameId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            inodeTag = int(f"{buildingID}{j+1}{i+1}")
            jnodeTag = int(f"{buildingID}{j+1}{i+2}")
            print(eleBeamTag, inodeTag, jnodeTag, ops.nodeCoord(inodeTag), ops.nodeCoord(jnodeTag))

            ops.element('nonlinearBeamColumn', eleBeamTag, *[inodeTag, jnodeTag], numItgrPts, BeamSecTag, BeamSecTransf)
    
    # define GRAVITY -------------------------------------------------------------
    # GRAVITY LOADS # define gravity load applied to beams and columns -- eleLoad applies loads in local coordinate axis
    linearTS = int(f"{buildingID}{1}")
    ops.timeSeries('Linear', linearTS)
    ops.pattern('Plain', int(f"{buildingID}{100}"), linearTS)

    # Column elements
    for j in range(NStory):
        for i in range(NBay + 1):
            # nodeTag = j * (NBay+1) + i + 1    # conventional node tag
            eleColTag = int(f"{buildingID}{colNameId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            ops.eleLoad('-ele', eleColTag, '-type', '-beamUniform', 0, -QdlCol)

    # Beam elements
    for j in range(1, NStory + 1): #because we dont have beam at ground level
        for i in range(NBay):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            eleBeamTag = int(f"{buildingID}{beamNameId}{j+1}{i+1}") #example: nodeTag 11 means 1st node of ground floor (1st floor)
            ops.eleLoad('-ele', eleBeamTag, '-type', '-beamUniform', -QdlBeam)    

def runGravityAnalysis(NStepGravity):
    Tol = 1.0e-8
    DGravity = 1/NStepGravity

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', Tol, 6)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', DGravity)
    ops.analysis('Static')

    status = ops.analyze(NStepGravity)

    ops.loadConst('-time', 0.)

    return status
    

def getRCFiberSection(buildingID, plotSection = False, **kwargs):
    ### GLOBAL VARIALBLES ###
    HCol = kwargs.get('HCol')
    BCol = kwargs.get('BCol')
    HBeam = kwargs.get('HBeam')
    BBeam = kwargs.get('BBeam')

    ### CHANGING GLOBAL PARAMETERS FOR DIFFERNET BUILDING IDS ###
    global ColSecTag, BeamSecTag, ColSecTransf, BeamSecTransf  # assigning global

    ### SECTION TAGS ###
    ColSecTag = int(f"{buildingID}{1}")
    BeamSecTag = int(f"{buildingID}{2}")
    ### TRANSFORMATION TAGS ###
    ColSecTransf = int(f"{buildingID}{1}")
    BeamSecTransf = int(f"{buildingID}{2}")

    ### MATERIAL PROPERTIES  ###
    # General Material parameters
    G = 1.e10		# make stiff shear modulus
    J = 1.0			# torsional section stiffness (G makes GJ large)
    GJ =  G*J

    # confined and unconfined CONCRETE
    # nominal concrete compressive strength
    fc = -25.0*N/mm**2	            # CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
    Ec = 5000*pow(-fc, .5)	# Concrete Elastic Modulus
    nu = 0.2
    Gc = Ec/(2.*(1+nu))	        # Torsional stiffness Modulus

    # confined concrete
    Kfc = 1.3			        # ratio of confined to unconfined concrete strength
    Kres = 0.2			        # ratio of residual/ultimate to maximum stress
    fc1C = Kfc*fc		        # CONFINED concrete (mander model), maximum stress
    eps1C = 2.*fc1C/Ec	        # strain at maximum stress 
    fc2C = Kres*fc1C		    # ultimate stress
    eps2C =  20*eps1C		    # strain at ultimate stress 
    lambda_ =  0.1			    # ratio between unloading slope at eps2 and initial slope Ec
    # unconfined concrete
    fc1U = fc			        # UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.003			    # strain at maximum strength of unconfined concrete
    fc2U = Kres*fc1U		    # ultimate stress
    eps2U = -0.01			    # strain at ultimate stress
    # tensile-strength properties
    ftC = -0.14*fc1C	        # tensile strength +tension
    ftU = -0.14*fc1U		    # tensile strength +tension
    Ets = ftU/0.002		        # tension softening stiffness

    # Core and Cover Concrete define
    IDconcCore = int(f"{buildingID}{1}")  # material id followed by building id
    IDconcCover = int(f"{buildingID}{2}") # material id followed by building id
    ops.uniaxialMaterial('Concrete02', IDconcCore, fc1C, eps1C, fc2C, eps2C, lambda_, ftC, Ets)	# Core concrete (confined)
    ops.uniaxialMaterial('Concrete02', IDconcCover, fc1U, eps1U, fc2U, eps2U, lambda_, ftU, Ets)	# Cover concrete (unconfined)

    # REINFORCING STEEL parameters
    Fy = 500*N/mm**2	        # STEEL yield stress
    Es = 2e11*N/m**2        	# modulus of steel
    Bs = 0.01			        # strain-hardening ratio 
    R0 = 18			            # control the transition from elastic to plastic branches
    cR1 = 0.925			        # control the transition from elastic to plastic branches
    cR2 = 0.15			        # control the transition from elastic to plastic branches

    IDSteel = int(f"{buildingID}{3}") # material id followed by building id
    ops.uniaxialMaterial('Steel02', IDSteel, Fy, Es, Bs, R0, cR1, cR2)

    ### RC FIBER SECTION PARAMETERS ###
    # Column section geometry:
    cover = 40*mm	            # rectangular-RC-Column cover
    numBarsTopCol = 4	        # number of longitudinal-reinforcement bars on top layer
    numBarsBotCol = 4		        # number of longitudinal-reinforcement bars on bottom layer
    numBarsIntCol = 8		        # TOTAL number of reinforcing bars on the intermediate layers
    barAreaTopCol = np.pi/4*(20*mm)**2	    # longitudinal-reinforcement bar area
    barAreaBotCol = np.pi/4*(20*mm)**2	    # longitudinal-reinforcement bar area
    barAreaIntCol = np.pi/4*(20*mm)**2	    # longitudinal-reinforcement bar area

    #Beam Section Geometry
    numBarsTopBeam = 4		        # number of longitudinal-reinforcement bars on top layer
    numBarsBotBeam = 4		        # number of longitudinal-reinforcement bars on bottom layer
    numBarsIntBeam = 8		        # TOTAL number of reinforcing bars on the intermediate layers
    barAreaTopBeam =  np.pi/4*(20*mm)**2	    # longitudinal-reinforcement bar area
    barAreaBotBeam =  np.pi/4*(20*mm)**2	    # longitudinal-reinforcement bar area
    barAreaIntBeam =  np.pi/4*(20*mm)**2	    # longitudinal-reinforcement bar area

    nfCoreY = 10		            # number of fibers in the core patch in the y direction
    nfCoreZ = 1	                    # number of fibers in the core patch in the z direction
    nfCoverY = 5		            # number of fibers in the cover patches with long sides in the y direction
    nfCoverZ = 1		            # number of fibers in the cover patches with long sides in the z direction
    # rectangular section with one layer of steel evenly distributed around the perimeter and a confined core.
    buildFiberSection.BuildRCrectSection(ColSecTag, HCol, BCol, cover, cover, IDconcCore, IDconcCover, 
                       IDSteel, numBarsTopCol, barAreaTopCol, numBarsBotCol, barAreaBotCol, numBarsIntCol,
                       barAreaIntCol, nfCoreY, nfCoreZ, nfCoverY, nfCoverZ, plotSection)
    
    buildFiberSection.BuildRCrectSection(BeamSecTag, HBeam, BBeam, cover, cover, IDconcCore, IDconcCover, 
                       IDSteel, numBarsTopBeam, barAreaTopBeam, numBarsBotBeam, barAreaBotBeam, numBarsIntBeam,
                       barAreaIntBeam, nfCoreY, nfCoreZ, nfCoverY, nfCoverZ, plotSection)
    
    ops.geomTransf(ColSecTransfType, ColSecTransf)
    ops.geomTransf('Linear', BeamSecTransf)

def setRayleighDamping():
    # define DAMPING--------------------------------------------------------------------------------------
    # apply Rayleigh DAMPING from $xDamp
    # D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial
    xDamp= 0.02;                        # 2% damping ratio
    lambda_ =  ops.eigen(1)[0]			# eigenvalue mode 1
    omega = pow(lambda_,0.5)
    alphaM =  0.				        # M-prop. damping; D = alphaM*M
    betaKcurr =  0.         			# K-proportional damping;      +beatKcurr*KCurrent
    betaKcomm =  2.*xDamp/omega  	    # K-prop. damping parameter;   +betaKcomm*KlastCommitt
    betaKinit =  0.		                # initial-stiffness proportional damping      +beatKinit*Kini
    # define damping
    ops.rayleigh(alphaM, betaKcurr, betaKinit, betaKcomm) # RAYLEIGH damping


def runGroundMotionAnalysis(gmData, GM_fact, dt):
    ops.loadConst('-time', 0.0)

    DtAnalysis = 0.01 #for analysis
    TmaxAnalysis = 10 #for analysis
    Nstep = int(TmaxAnalysis/DtAnalysis)
    
    GM_dirn = 1
    gmTS = 2
    ops.timeSeries("Path", gmTS, '-dt', dt, '-values', *gmData)
    #pattern('UniformExcitation', patternTag, dir, '-disp', dispSeriesTag, '-vel', velSeriesTag, '-accel', accelSeriesTag, '-vel0', vel0, '-fact', fact)
    ops.pattern('UniformExcitation', 300, GM_dirn, '-accel', gmTS, '-fact', GM_fact)

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("ProfileSPD")
    tol = 1.e-8
    maxNumIter = 50
    ops.test("EnergyIncr", tol, maxNumIter)
    ops.algorithm("ModifiedNewton")
    ops.integrator("Newmark", .5, .25)
    ops.analysis("Transient")

    for i in range(Nstep):
        status = ops.analyze(1, DtAnalysis)
        # ovs.plot_defo()
        # plt.savefig(f'test_images/{(i+1)*DtAnalysis}s_{TmaxAnalysis}.png')
        # plt.close()
        print(f"Ground Motion - {i+1}/{Nstep}")
        if status != 0:
            print("Analysis failed trying Krylov Newton...")
            ops.algorithm('KrylovNewton')
            status = ops.analyze(1, DtAnalysis)
        if status != 0:
            print("Analysis failed trying with more iterations and less tolerance...")
            tol = 1.e-6
            maxNumIter = 100
            ops.system("BandSPD")
            ops.test("RelativeNormDispIncr", tol, maxNumIter)
            status = ops.analyze(1, DtAnalysis)
        if status!=0:
            print("Analysis failed trying Bryoden Algorithm")
            tol = 1.e-8
            maxNumIter = 50
            ops.test("EnergyIncr", tol, maxNumIter)
            ops.algorithm('Broyden')
            status = ops.integrator("HHT", 0.85, .5, .25)
            status = ops.analyze(1, DtAnalysis)
            #
        if status != 0:
            print("Breaking analysis")
            break
        
        ops.algorithm('ModifiedNewton')
        ops.system("ProfileSPD")
        tol = 1.e-10
        maxNumIter = 50    
        ops.test("EnergyIncr", tol, maxNumIter)
        

def kelvinVoigtMaterials(idKelvin, LBuildingRNodes, RBuildingLNodes, gap):
    # viscous material
    viscousID = int(f"{idKelvin}{1}")
    C = 683 * kN*sec/m
    alpha = 1
    ops.uniaxialMaterial('Viscous', viscousID, C, alpha)

    # spring materil
    springID = int(f"{idKelvin}{2}")
    E0 = 93500 * kN/m**2
    ops.uniaxialMaterial('Elastic', springID, E0)

    # epp GAP
    eppGAPMatID = int(f"{idKelvin}{3}")
    E = 100* E0
    Fy = 250*Mpa
    eta = 0.1
    ops.uniaxialMaterial('ElasticPPGap', eppGAPMatID, 1*E, -1*Fy, -1*gap, eta, 'damage')

    ### kelvin voigt construction
    parallelTag = int(f"{idKelvin}{4}")
    ops.uniaxialMaterial('Parallel', parallelTag, *[viscousID, springID])

    # creating new nodes for kelvin voigt plus the epp material
    adjacent_nodes = []   
    for rNode, lNode in zip(LBuildingRNodes, RBuildingLNodes): #zipping it will help to create less no. of nodes
        node_tag_kv = int(f"{900}{rNode}")
        ops.node(node_tag_kv, *ops.nodeCoord(rNode))
        adjacent_nodes.append(node_tag_kv)
        ops.equalDOF(rNode, node_tag_kv, *[2,3]) #to remove problems
    
    kvEleID = []
    eppEleID = []
    for lNode, midNode, rNode in zip(LBuildingRNodes, adjacent_nodes, RBuildingLNodes):
        #element('zeroLength', eleTag, *eleNodes, '-mat', *matTags, '-dir', *dirs, <'-doRayleigh', rFlag=0>, <'-orient', *vecx, *vecyp>)
        kvEleTag = int(f"{9000}{lNode}")
        ops.element('zeroLength', kvEleTag, *[lNode, midNode], '-mat', parallelTag, '-dir', *[1])
        #element('twoNodeLink', eleTag, *eleNodes, '-mat', *matTags, '-dir', *dir, <'-orient', *vecx, *vecyp>, <'-pDelta', *pDeltaVals>, <'-shearDist', *shearDist>, <'-doRayleigh'>, <'-mass', m>)
        eppEleTag = int(f"{9000}{rNode}")
        ops.element('twoNodeLink', eppEleTag, *[midNode, rNode], '-mat', eppGAPMatID, '-dir', *[1])
        kvEleID.append(kvEleTag)
        eppEleID.append(eppEleTag)

    ops.recorder('Element', '-file', 'testForceKelvinVoigt.txt', '-time', '-closeOnWrite', '-ele', *kvEleID,'-dof',1, 'force')
    ops.recorder('Element', '-file', 'testForceSpringEPPGAP.txt', '-time', '-closeOnWrite', '-ele', *eppEleID,'-dof',1, 'force')
    print(kvEleID)







ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)
#LNodes10 = []
#RNodes10 = []

gap = 30*mm
building3storey = {'HCol':350*mm, 'BCol':330*mm, 'HBeam':350*mm, 'BBeam':250*mm}
building5storey = {'HCol':450*mm, 'BCol':450*mm, 'HBeam':450*mm, 'BBeam':300*mm}
building7storey = {'HCol':550*mm, 'BCol':550*mm, 'HBeam':550*mm, 'BBeam':350*mm}

getModel(buildingID=10, NBay=1, NStory=7, LBeam=4, LCol=3, startCoor=(0,0), **building7storey)
getModel(buildingID=20, NBay=1, NStory=3, LBeam=4, LCol=3, startCoor=(4+gap,0), **building3storey)
ovs.plot_model()
plt.show()

print(LNodes10, RNodes10)
print(LNodes20, RNodes20)

for a,b in zip(RNodes10, LNodes20):
    print(a,b)

lomaPrietaEq2 = "data/0493a.smc"
lomaPrietaEq = "data/A10000.dat"
casi68Eq = "data/BM68elc.dat"
elcentro = "data/elcentro.txt"

setRayleighDamping()

runGravityAnalysis(100)
ovs.plot_defo()
plt.show()

kelvinVoigtMaterials(100, RNodes10, LNodes20, gap)
print("Added Kelvin Voigt Impact Element")

ops.recorder('Node', '-file', "Building10RightNodes_Disp.txt", '-time', '-closeOnWrite', '-node', *RNodes10,'-dof', 1, 'disp')
ops.recorder('Node', '-file', "Building20LeftNodes_Disp.txt", '-time', '-closeOnWrite', '-node', *LNodes20,'-dof',1, 'disp')

gmDatacasi68Eq = np.loadtxt(casi68Eq).ravel()
gmDataLomaPrietaEq = np.loadtxt(lomaPrietaEq).ravel()
gmDataLomaPrietaEq2 = np.loadtxt(lomaPrietaEq2, comments='#').ravel()


runGroundMotionAnalysis(gmDataLomaPrietaEq, GM_fact=9.81, dt=0.01)