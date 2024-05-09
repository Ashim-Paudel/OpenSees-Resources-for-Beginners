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


data = np.loadtxt("Pounding Models/linear_springDisp.txt")
plt.plot(-data[1], data[0])
plt.show()
#print(data)