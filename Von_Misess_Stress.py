# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 06:19:17 2020

@author: HJ Hoogendoorn

Maximum stress calculations using Von 
mises stress
"""

l_a = 1.691 #[m] span of aileron
from Moment_of_Inertia import *

from main import *


#integration steps

dx = 0.001 
dy = 0.001
dz = 0.001


sigma_x = [] 
sigma_y = []
sigma_z = []


tau_xy = 0
tau_yz = 1


    

Y = []

   
#print(max(Y))
#Y.index(max(Y))


#stresses in cross section
for x in range(0,int(l_a/dx)):
    for i in range(len(ystr)):
        sigma_x.append((-R1y*Mac(x-x1,1) - R2y*Mac(x-x2,1) - R3y*Mac(x-x3,1) - sin(theta)*Ra1*Mac(x-xa1,1) + sin(theta)*P*Mac(x-xa2,1) )/Izz*ystr[i] + (R1z*Mac(x-x1,1) + R2z*Mac(x-x2,1) + R3z*Mac(x-x3,1) + cos(theta)*Ra1*Mac(x-xa1,1) - cos(theta)*P*Mac(x-xa2,1))/Iyy*zstr[i])
        sigma_y.append((-R1y*Mac(x-x1,1) - R2y*Mac(x-x2,1) - R3y*Mac(x-x3,1) - sin(theta)*Ra1*Mac(x-xa1,1) + sin(theta)*P*Mac(x-xa2,1) )/Izz*x)
        sigma_z.append((R1z*Mac(x-x1,1) + R2z*Mac(x-x2,1) + R3z*Mac(x-x3,1) + cos(theta)*Ra1*Mac(x-xa1,1) - cos(theta)*P*Mac(x-xa2,1))/Iyy*x)

for k in range(0,len(sigma_x)):
    Y.append(np.sqrt(0.5*((sigma_x[k]-sigma_y[k])**2+(sigma_y[k]-sigma_z[k])**2+(sigma_z[k]-sigma_x[k])**2)+3*tau_yz**2))

print(max(Y))   



#print(max(Y))
#Y.index(max(Y))



