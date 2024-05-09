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


# Section Geometry:
HCol = 24*inch	# square-Column width
BCol = HCol
HBeam = 42*inch	# Beam depth -- perpendicular to bending axis
BBeam = 24*inch	# Beam width -- parallel to bending axis


# Additional tags for Inelastic Beam Column Section
ColMatTagAxial = 101                              # Represent column axial behaviour
ColMatTagFlex = 102                               # Represent column flexural behaviour
BeamMatTagAxial = 201                             # Represent beam axial behaviour
BeamMatTagFlex = 202                              # Represent beam flexural behaviour 

# Geometric Transformation Tags
ColSecTransf = 1
BeamSecTransf = 2
ColSecTransfType = 'Linear'

# model building function
def getModel(buildingID, NBay, NStory, LBeam, LCol, sectionType = 'Elastic', startCoor = (0.0 ,0.0)):
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
    global HCol, BCol, HBeam, BBeam

    #ops.wipe()

    #ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)

    for j in range(NStory + 1):
        for i in range(NBay + 1):
            # nodeTag = j * (NBay+1) + i + 1 # conventional node tag
            nodeTag = int(f"{buildingID}{j+1}{i+1}") #example: nodeTag 10011 means 1st building(100) 1st node of ground floor (1st floor) 
            nodeCoor = (startCoor[0] + i * LBeam, startCoor[1] + j * LCol)
            
            print(nodeTag, *nodeCoor)
            ops.node(nodeTag, *(nodeCoor))

    # fixity to the nodes
    ops.fixY(0.0, 1, 1, 1)

    # determine support nodes where ground motions are input, for multiple-support excitation
    level = 1
    iSupportNode = [int(f"{buildingID}{level}{i+1}") for i in range(0, NBay + 1 )]

    # up parameters that are particular to the model for displacement control
    IDctrlNode =  (NStory+1)*10+1		# node where displacement is read for displacement control
    IDctrlDOF =  1	                    # DoF of displacement read for displacement control
    LBuilding =  NStory*LCol	        # total building height



    # Define Section and weight of section per unit length of element
    # beam and column section weight for all materials except rc concrete
    QBeam = 94*lbf/ft		            # W-section weight per length
    QdlCol = 114*lbf/ft	                # W-section weight per length

    GammaConcrete = 150*pcf   		    # Reinforced-Concrete floor slabs
    if sectionType == 'Elastic':
        getElasticSection()

    elif sectionType == 'InElastic':
        getInelasticSection()

    elif sectionType == 'RCFiber':
        getRCFiberSection(buildingID)
        QBeam = GammaConcrete*HBeam*BBeam	# self weight of Beam, weight per length
        QdlCol = GammaConcrete*HCol*BCol	# self weight of Column, weight per length
        
    elif sectionType == 'SteelFiber':
        getSteelFiberSection(buildingID)


    #### GRAVITY LOADS, MASSES AND WEIGHT ###
    Tslab = 6*inch			            # 6-inch slab
    Lslab = 2*LBeam/2 			        # assume slab extends a distance of LBeam1/2 in/out of plane
    Qslab = GammaConcrete*Tslab*Lslab   # slab dead weight 
    QdlBeam = Qslab + QBeam 	        # dead load distributed along beam.
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
    
    #### GRAVITY LOADS, MASSES AND WEIGHT ###
    #GammaConcrete = 150*pcf   		# Reinforced-Concrete floor slabs
    #Tslab = 6*inch			# 6-inch slab
    #Lslab = 2*LBeam/2 			# assume slab extends a distance of LBeam1/2 in/out of plane
    #Qslab = GammaConcrete*Tslab*Lslab 
    #Qslab
    #QdlCol = GammaConcrete*HCol*BCol	# self weight of Column, weight per length
    #QBeam = GammaConcrete*HBeam*BBeam	# self weight of Beam, weight per length
    #QdlBeam = Qslab + QBeam 	# dead load distributed along beam.
    #WeightCol = QdlCol*LCol  		# total Column weight
    #WeightBeam = QdlBeam*LBeam 	# total Beam weight

    #Qslab = GammaConcrete*Tslab*Lslab
    #QdlCol = GammaConcrete*HCol*BCol 
    #QBeam = 94*lbf/ft		            # W-section weight per length
    #QdlBeam = Qslab + QBeam	            # dead load distributed along beam.
    #QdlCol = 114*lbf/ft	                # W-section weight per length
    #WeightCol = QdlCol*LCol     		# total Column weight
    #WeightBeam = QdlBeam*LBeam      	# total Beam weight
                
def getElasticSection():
    # section geometry
    # column sections: W27x114
    AgCol = 33.5*sqinch	                # cross-sectional area
    IzCol = 4090.*inch4	                # moment of Inertia
    # beam sections: W24x94
    AgBeam = 27.7*sqinch		        # cross-sectional area
    IzBeam = 2700.*inch4	            # moment of Inertia

    # material properties
    Es = 29000*ksi		# Steel Young's Modulus

    ops.geomTransf(ColSecTransfType, ColSecTransf)
    ops.geomTransf('Linear', BeamSecTransf)

    ops.section('Elastic', BeamSecTag, Es, AgBeam, IzBeam)
    ops.section('Elastic', ColSecTag, Es, AgCol, IzCol)


def getInelasticSection():
    # MATERIAL properties 
    Fy =  6.0*ksi       # Yield Stress 
    Es =  29000*ksi		# Steel Young's Modulus
    nu =  0.3           # Poisson's Ratio
    Gs =  Es/(2.*(1+nu)) 	# Torsional stiffness Modulus
    
    # SECTION PROPERTIES
    # COLUMN section W27x114
    AgCol = 33.5*pow(inch,2)		# cross-sectional area
    IzCol = 4090.*pow(inch,4)		# moment of Inertia
    EICol = Es*IzCol				# EI, for moment-curvature relationship
    EACol = Es*AgCol				# EA, for axial-force-strain relationship
    MyCol = 2e4*kip*inch	   		# yield moment
    PhiYCol = 0.25e-3/inch	   		# yield curvature
    PhiYCol = MyCol/EICol			# yield curvature
    EIColCrack = MyCol/PhiYCol		# cracked section inertia
    b = 0.01  					    # strain-hardening ratio 

    ops.uniaxialMaterial('Steel01', ColMatTagFlex, MyCol, EIColCrack, b)
    ops.uniaxialMaterial('Elastic', ColMatTagAxial, EACol)
    ops.section('Aggregator', ColSecTag, *[ColMatTagAxial, 'P', ColMatTagFlex, 'Mz'])

    # BEAM SECTION W24x94
    AgBeam = 27.7*pow(inch,2)		# cross-sectional area
    IzBeam = 2700.*pow(inch,4)		# moment of Inertia
    EIBeam = Es*IzBeam				# EI, for moment-curvature relationship
    EABeam = Es*AgBeam				# EA, for axial-force-strain relationship
    MyBeam = 1.5e4*kip*inch	   		# yield moment
    PhiYBeam = 0.25e-3/inch	   		# yield curvature
    PhiYBeam = MyBeam/EIBeam			# yield curvature
    EIBeamCrack = MyBeam/PhiYBeam		# cracked section inertia
    b = 0.01  					    # strain-hardening ratio 

    ops.uniaxialMaterial('Steel01', BeamMatTagFlex, MyBeam, EIBeamCrack, b)  # bilinear behavior for flexure
    ops.uniaxialMaterial('Elastic', BeamMatTagAxial, EABeam)
    ops.section('Aggregator', BeamSecTag, *[BeamMatTagAxial, 'P', BeamMatTagFlex, 'Mz'])



def getRCFiberSection(buildingID, plotSection = False):
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
    fc = -4.0*ksi		        # CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
    Ec = 57*ksi*pow((-fc/psi), .5)	# Concrete Elastic Modulus
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
    Fy = 66.8*ksi		        # STEEL yield stress
    Es = 29000.*ksi		        # modulus of steel
    Bs = 0.01			        # strain-hardening ratio 
    R0 = 18			            # control the transition from elastic to plastic branches
    cR1 = 0.925			        # control the transition from elastic to plastic branches
    cR2 = 0.15			        # control the transition from elastic to plastic branches

    IDSteel = int(f"{buildingID}{3}") # material id followed by building id
    ops.uniaxialMaterial('Steel02', IDSteel, Fy, Es, Bs, R0, cR1, cR2)

    ### RC FIBER SECTION PARAMETERS ###
    # Column section geometry:
    cover = 2.5*inch	            # rectangular-RC-Column cover
    numBarsTopCol = 8		        # number of longitudinal-reinforcement bars on top layer
    numBarsBotCol = 8		        # number of longitudinal-reinforcement bars on bottom layer
    numBarsIntCol = 6		        # TOTAL number of reinforcing bars on the intermediate layers
    barAreaTopCol = 1.*sqinch	    # longitudinal-reinforcement bar area
    barAreaBotCol = 1.*sqinch	    # longitudinal-reinforcement bar area
    barAreaIntCol = 1.*sqinch	    # longitudinal-reinforcement bar area

    #Beam Section Geometry
    numBarsTopBeam = 6		        # number of longitudinal-reinforcement bars on top layer
    numBarsBotBeam = 6		        # number of longitudinal-reinforcement bars on bottom layer
    numBarsIntBeam = 2		        # TOTAL number of reinforcing bars on the intermediate layers
    barAreaTopBeam =  1.*sqinch	    # longitudinal-reinforcement bar area
    barAreaBotBeam =  1.*sqinch	    # longitudinal-reinforcement bar area
    barAreaIntBeam =  1.*sqinch	    # longitudinal-reinforcement bar area

    nfCoreY = 20		            # number of fibers in the core patch in the y direction
    nfCoreZ = 20	            # number of fibers in the core patch in the z direction
    nfCoverY = 20		            # number of fibers in the cover patches with long sides in the y direction
    nfCoverZ = 20		            # number of fibers in the cover patches with long sides in the z direction
    # rectangular section with one layer of steel evenly distributed around the perimeter and a confined core.
    buildFiberSection.BuildRCrectSection(ColSecTag, HCol, BCol, cover, cover, IDconcCore, IDconcCover, 
                       IDSteel, numBarsTopCol, barAreaTopCol, numBarsBotCol, barAreaBotCol, numBarsIntCol,
                       barAreaIntCol, nfCoreY, nfCoreZ, nfCoverY, nfCoverZ, plotSection)
    
    buildFiberSection.BuildRCrectSection(BeamSecTag, HBeam, BBeam, cover, cover, IDconcCore, IDconcCover, 
                       IDSteel, numBarsTopBeam, barAreaTopBeam, numBarsBotBeam, barAreaBotBeam, numBarsIntBeam,
                       barAreaIntBeam, nfCoreY, nfCoreZ, nfCoverY, nfCoverZ, plotSection)
    
    ops.geomTransf(ColSecTransfType, ColSecTransf)
    ops.geomTransf('Linear', BeamSecTransf)


def getSteelFiberSection(buildingID, plotSection = False):
    ### CHANGING GLOBAL PARAMETERS FOR DIFFERNET BUILDING IDS ###
    global ColSecTag, BeamSecTag, ColSecTransf, BeamSecTransf  # assigning global

    ### SECTION TAGS ###
    ColSecTag = int(f"{buildingID}{1}")
    BeamSecTag = int(f"{buildingID}{2}")
    ### TRANSFORMATION TAGS ###
    ColSecTransf = int(f"{buildingID}{1}")
    BeamSecTransf = int(f"{buildingID}{2}")

    
    #### MATERIAL PROPERTIES ####
    Fy = 60.0*ksi
    Es = 29000*ksi		# Steel Young's Modulus
    nu = 0.3
    Gs = Es/(2.*(1+nu)) # Torsional stiffness Modulus
    Hiso = 0
    Hkin = 1000
    matIDhard = int(f"{buildingID}1")  # material id followed by building id
    # uniaxialMaterial('Hardening', matTag, E, sigmaY, H_iso, H_kin, eta=0.0)
    ops.uniaxialMaterial('Hardening', matIDhard, Es, Fy, Hiso, Hkin)

    #### SECTION PROPERTIES ####
    # Structural-Steel W-section properties
    # COLUMN SECTION: W27x114
    dCol = 27.29*inch	        # depth
    bfCol = 10.07*inch	        # flange width
    tfCol = 0.93*inch	        # flange thickness
    twCol = 0.57*inch	        # web thickness
    nfdw = 16		            # number of fibers along dw
    nftw = 2		            # number of fibers along tw
    nfbf = 16		            # number of fibers along bf
    nftf = 4			        # number of fibers along tf

    buildFiberSection.BuildSteelWSection(ColSecTag, matIDhard, 
                                         dCol, bfCol, tfCol, twCol, nfdw, nftw, nfbf, nftf, plotSection)
    # BEAM SECTION: W24x94
    dBeam =  24.31*inch	        # depth
    bfBeam =  9.065*inch	    # flange width
    tfBeam =  0.875*inch	    # flange thickness
    twBeam = 0.515*inch	        # web thickness
    nfdw = 16		            # number of fibers along dw
    nftw = 2		            # number of fibers along tw
    nfbf = 16		            # number of fibers along bf
    nftf = 4			        # number of fibers along tf
    buildFiberSection.BuildSteelWSection(BeamSecTag, matIDhard, 
                                         dBeam, bfBeam, tfBeam, twBeam, nfdw, nftw, nfbf, nftf, plotSection)

    ops.geomTransf(ColSecTransfType, ColSecTransf)
    ops.geomTransf('Linear', BeamSecTransf)


ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)
getModel(buildingID=10, NBay=1, NStory=4, LBeam=4, LCol=4, sectionType='RCFiber', startCoor=(0,0))
getModel(buildingID=20, NBay=1, NStory=3, LBeam=4, LCol=3, sectionType='SteelFiber', startCoor=(5,0))
getModel(buildingID=30, NBay=1, NStory=5, LBeam=5, LCol=3, sectionType='SteelFiber', startCoor=(10,0))
#getModel(buildingID=40, NBay=1, NStory=7, LBeam=3, LCol=3, sectionType='RCFiber', startCoor=(14,0))
#getModel(buildingID=50, NBay=1, NStory=3, LBeam=4, LCol=3, sectionType='SteelFiber', startCoor=(18,0))
ovs.plot_model()
plt.show()