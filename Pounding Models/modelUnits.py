#length conversion
m = 1.0
cm = 100*m # centimeter, needed for displacement input in MultipleSupport excitation
inch = m/39.37
ft = 12.*inch 

#area conversion
sqinch = inch**2
inch4 = inch**4

#ksi to pa
pa = 1.0
kpa = 1000*pa
Mpa = 1000*kpa
ksi = 6.89*1000*kpa  #ksi = kip / sqinch
psi = ksi/1000.


#kip to kN
N = 1.0
kN = 1000 * N
kip = 4.45 * kN

# pound to kg
kg = 1 
pound = 0.453592 * kg
lbf = psi*inch**2  #pounds force
pcf = lbf/pow(ft, 3) #pounds per cubic foot
psf = lbf/pow(ft, 2) #pounds per square foot


# acceleration due to gravity
g = 9.81 #m/s^2

#time
sec = 1.0
