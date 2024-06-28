# Import Packages #
from openseespy import opensees as ops
import opsvis as ovs
import numpy as np
import matplotlib.pyplot as plt 

#units conversion
m = 1.0
mm = 0.001 * m

Pa = 1.0
kPa = 1000 * Pa
MPa = 1000 * kPa
GPa = 1000 * MPa

N = 1.0
kN = 1000*N

def analysisLoopCyclicLCD(status, i, stepSize):
    """
    The load control analysis loop.
    """
    
    if status != 0:
        print("Trying 5 times smaller timestep at load factor", i)
        ops.integrator("LoadControl", stepSize/5)
        status = ops.analyze(1)
    
    if status != 0:
        print("Trying 20 times smaller timestep at load factor", i)
        ops.integrator("LoadControl", stepSize/20)
        status = ops.analyze(1)        
        
    if status != 0:
        print("Trying 80 times smaller timestep at load factor", i)
        ops.integrator("LoadControl", stepSize/80)
        status = ops.analyze(1)       
        
    if status != 0:
        print("Trying 160 times smaller timestep at load factor", i)
        ops.integrator("LoadControl", stepSize/160)
        status = ops.analyze(1)
        
    if status != 0:
        print("Trying 200 interations at load factor", i)
        ops.test('NormDispIncr', 1.*10**-8, 200)
        status = ops.analyze(1)
        
    if status != 0:
        print("Trying ModifiedNewton at load factor", i)
        ops.algorithm("ModifiedNewton")
        ops.test('NormDispIncr', 1.*10**-8, 200)
        status = ops.analyze(1)
    
    ops.test('NormDispIncr', 1.*10**-8, 50)
    ops.integrator("LoadControl", stepSize)
    ops.algorithm("Newton")
    return status


def cyclicAnalysis_LCD(dispMax, filename = "LoadProtocols/Ganey_small.thf"):
    #model parameters
    controlNode = 4
    controlNodeDof = 1
    du = 1.0*m
    
    #time series
    ops.timeSeries('Constant', 1)
    ops.timeSeries('Linear', 2)
    #timeSeries('Path', tag, '-dt', dt=0.0, '-values', *values, '-time', *time, '-filePath', filePath='', '-fileTime', fileTime='', '-factor', factor=1.0, '-startTime', startTime=0.0, '-useLast', '-prependZero')
    ops.timeSeries('Path', 3, '-dt', 1.0, '-filePath', filename, '-factor', 1.0, '-perpendZero')
    
    #load pattern
    loadPatternTag = 1                   #timeSeries tag
    ops.pattern('Plain', loadPatternTag, 3)
    ops.sp(controlNode, controlNodeDof, du)
    #ops.load(controlNode, *[P0_x, 0.0, 0.0]

    #using 0.1 time steps , but we have dt of 1.0sec so opensees will interpolate to generate extra 10 points btwn each time sets
    stepSize = 0.1
    
    #analysis commands
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    #integrator('LoadControl', incr, numIter=1, minIncr=incr, maxIncr=incr)
    ops.integrator('LoadControl', stepSize, 1, stepSize, stepSize*100)
    ops.algorithm("Newton")
    ops.analysis("Static")

    #ops.test('NormDispIncr', tol, iter, pFlag=0, nType=2)
    ops.test('NormDispIncr', 1.0e-8, 50)
    
    #recorder('Node', '-file', filename, , '-timeSeries', tsTag, '-time','-closeOnWrite', '-node', *nodeTags=[],'-dof', *dofs=[], respType)
    ops.recorder('Node', '-file', "Outputs/CyclicLCD/displacements.txt",'-closeOnWrite', '-time', '-node', *[4], '-dof', *[1,2,3], 'disp')
    ops.recorder('Node', '-file', "Outputs/CyclicLCD/reactions.txt", '-closeOnWrite', '-time', '-node', *[1], '-dof', *[1,2,3], 'reaction')
    
    #analyze(numIncr=1, dt=0.0, dtMin=0.0, dtMax=0.0, Jd=0)
    ops.record()

    #recording time and loadFactors at each steps
    time = []
    loadFactors = []
    i = 0 #variable to track step of each load factor

    while (ops.nodeDisp(controlNode, controlNodeDof) < dispMax):
        status = ops.analyze(1)
        if status != 0:
            status = analysisLoopCyclicLCD(status, i, stepSize)
        if status != 0:
            print("Analysis Failed at loadFactor ", i)
            break
        i += 1
        
        time.append(ops.getTime())
        loadFactors.append(ops.getLoadFactor(loadPatternTag))
    
    if status == 0:
        print("Analysis Successful")

    return time, loadFactors



def analysisLoopCyclicDCD(status, i, dx, controlNode, controlNodeDof):
    """
    The displacement control analysis loop.
    """
    
    if status != 0:
        print("Trying 5 times smaller timestep at load factor", i)
        ops.integrator('DisplacementControl', controlNode, controlNodeDof, dx/5)
        status = ops.analyze(1)
    
    if status != 0:
        print("Trying 20 times smaller timestep at load factor", i)
        ops.integrator('DisplacementControl', controlNode, controlNodeDof, dx/20)
        status = ops.analyze(1)        
        
    if status != 0:
        print("Trying 80 times smaller timestep at load factor", i)
        ops.integrator('DisplacementControl', controlNode, controlNodeDof, dx/80)
        status = ops.analyze(1)       
        
    if status != 0:
        print("Trying 160 times smaller timestep at load factor", i)
        ops.integrator('DisplacementControl', controlNode, controlNodeDof, dx/160)
        status = ops.analyze(1)
        
    if status != 0:
        print("Trying 1000 times smaller timestep at load factor", i)
        ops.integrator('DisplacementControl', controlNode, controlNodeDof, dx/1000)
        status = ops.analyze(1)   
        
    if status != 0:
        print("Trying ModifiedNewton at load factor", i)
        ops.algorithm("ModifiedNewton")
        ops.test('NormDispIncr', 1.*10**-6, 200)
        status = ops.analyze(1)

    ops.integrator('DisplacementControl', controlNode, controlNodeDof, dx)
    ops.test('NormDispIncr', 1.*10**-10, 50)
    ops.algorithm("Newton")
    return status


def cyclicAnalysis_DCD(loadProtocol = [0.02,0.05], Nrepeat = [2,2], dx = 0.0001*m):
    #model parameters
    controlNode = 4
    controlNodeDof = 1
    dForce = 1*kN
    
    #time series
    ops.timeSeries('Constant', 1)  
    #load pattern
    loadPatternTag = 1                   #timeSeries tag
    ops.pattern('Plain', loadPatternTag, 1)
    ops.load(controlNode, *[dForce, 0.0, 0.0])
    #ops.sp(controlNode, controlNodeDof, du)

    #using 0.1 time steps , but we have dt of 1.0sec so opensees will interpolate to generate extra 10 points btwn each time sets
    stepSize = 0.1
    
    #analysis commands
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.test('NormDispIncr', 1.0e-10, 50)  #ops.test('NormDispIncr', tol, iter, pFlag=0, nType=2)
    
    #recorder('Node', '-file', filename, , '-timeSeries', tsTag, '-time','-closeOnWrite', '-node', *nodeTags=[],'-dof', *dofs=[], respType)
    ops.recorder('Node', '-file', "Outputs/CyclicDCD/displacements.txt",'-closeOnWrite', '-time', '-node', *[4], '-dof', *[1,2,3], 'disp')
    ops.recorder('Node', '-file', "Outputs/CyclicDCD/reactions.txt", '-closeOnWrite', '-time', '-node', *[1], '-dof', *[1,2,3], 'reaction')
    
    #analyze(numIncr=1, dt=0.0, dtMin=0.0, dtMax=0.0, Jd=0)
    ops.record()

    i = 0 #variable to track step of each load factor
    #integrator('DisplacementControl', nodeTag, dof, incr, numIter=1, dUmin=incr, dUmax=incr)
    for x, Ncycle in zip(loadProtocol, Nrepeat):
        for cycle in range(Ncycle):
            ops.integrator('DisplacementControl', controlNode, controlNodeDof, dx, 10, dx/1000, dx*10)
            while (ops.nodeDisp(controlNode, controlNodeDof) < x):
                status = ops.analyze(1)
                i += 1
                
                if status != 0:
                    status = analysisLoopCyclicDCD(status, i, dx, controlNode, controlNodeDof) 
                if status != 0:
                    print("Breaking Analysis! Analysis Failed")
                    return
            #negative cycle
            ops.integrator('DisplacementControl', controlNode, controlNodeDof, -dx, 10, -dx/1000, -dx*10)
            while(ops.nodeDisp(controlNode, controlNodeDof) > -x):
                status = ops.analyze(1)
                i += 1
                if status != 0:
                    status = analysisLoopCyclicDCD(status, i, -dx, controlNode, controlNodeDof)
                if status != 0:
                    print("Breaking Analysis! Analysis Failed")
                    return
            print(x, i)
                
                