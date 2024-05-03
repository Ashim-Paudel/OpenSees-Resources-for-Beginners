# Import Packages #
from openseespy import opensees as ops
import opsvis as ovs
import numpy as np
import matplotlib.pyplot as plt 

import model_functions as mf
import analysis_functions as af
import output_functions as of

#Displacement Control
loadProtocol = [0.001547,  0.00221, 0.003315, 0.00442, 0.00663, 0.00884, 0.01326,
                0.01768,   0.02652, 0.03978,  0.05967, 0.0884,  0.1326,  0.1768, 0.221, 
                0.3315,    0.4199,  0.4862]
Nrepeat      = 3*np.ones_like(loadProtocol, int) 

#analysis name for displacement control displacement analysis
DCD = 'CyclicDCD'

#analysis codes
ops.wipe()
mf.getSections()
mf.getModel()
of.getRecorders(DCD)
af.cyclicAnalysis_DCD(loadProtocol, Nrepeat)
of.plotPushover(DCD)
