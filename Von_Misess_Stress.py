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
import numpy as np



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



#stresses in cross section
for x in range(0,int(l_a/dx)):
    for i in range(len(ystr)):
        sigma_x.append(Mz(x)/Izz*ystr[i] + (My(x))/Iyy*zstr[i])
        sigma_y.append((Mz(x))/Izz*x)
        sigma_z.append((My(x))/Iyy*x)

for k in range(0,len(sigma_x)):
    Y.append(np.sqrt(0.5*((sigma_x[k]-sigma_y[k])**2+(sigma_y[k]-sigma_z[k])**2+(sigma_z[k]-sigma_x[k])**2)+3*tau_yz**2))

print(max(Y))   



#print(max(Y))
#Y.index(max(Y))



