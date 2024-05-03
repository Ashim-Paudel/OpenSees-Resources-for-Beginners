# Import Packages #
from openseespy import opensees as ops
import opsvis as ovs
import numpy as np
import matplotlib.pyplot as plt 


def getRecorders(AnalysisType):
    disp_filename = 'Outputs/' + AnalysisType + '/displacements.txt'
    reaction_filename = 'Outputs/' + AnalysisType + '/reactions.txt'
    #recorder('Node', '-file', filename, , '-timeSeries', tsTag, '-time','-closeOnWrite', '-node', *nodeTags=[],'-dof', *dofs=[], respType)
    ops.recorder('Node', '-file', disp_filename,'-closeOnWrite', '-time',
                  '-node', *[4], '-dof', *[1,2,3], 'disp')
    ops.recorder('Node', '-file', reaction_filename, '-closeOnWrite', '-time',
                  '-node', *[1], '-dof', *[1,2,3], 'reaction')


def plotPushover(AnalysisType):
    disp_filename = 'Outputs/' + AnalysisType + '/displacements.txt'
    reaction_filename = 'Outputs/' + AnalysisType + '/reactions.txt'

    displacements = np.loadtxt(disp_filename, delimiter=" ")
    reactions = np.loadtxt(reaction_filename, delimiter=" ")

    disp_x = displacements[:, 1]
    base_shear = -reactions[:, 1]

    fig, ax = plt.subplots()
    ax.plot(disp_x, base_shear)
    fig.savefig(f"BackBone Curve - Cyclic {AnalysisType}")

    
