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


sigma_x = 0 #no stresses in x direction
sigma_y = 1
sigma_z = 1
#maximum shear occurs on z' axis, so only 3 possibilities:
#1: z' = 0
#2: z' = ha/2
#3: z' = Ca

tau_xy = 0
tau_yz_1 = 1
tau_yz_2 = 1.1
tau_yz_3 = 1.3


    

Y = []

for x in range(0,int(l_a/dx)):
    sigma_y = (-R1y*Mac(x-x1,1) - R2y*Mac(x-x2,1) - R3y*Mac(x-x3,1) - sin(theta)*Ra1*Mac(x-xa1,1) + sin(theta)*P*Mac(x-xa2,1) )/Izz
    sigma_z = (R1z*Mac(x-x1,1) + R2z*Mac(x-x2,1) + R3z*Mac(x-x3,1) + cos(theta)*Ra1*Mac(x-xa1,1) - cos(theta)*P*Mac(x-xa2,1))/Iyy
    tau_max = max(tau_yz_1*x,tau_yz_2*x,tau_yz_3*x)
    
    
    
    #maximum stress in aileron
    Y.append(np.sqrt(0.5*((sigma_y)**2+(sigma_y-sigma_z)**2+(sigma_z)**2)+3*tau_max**2))
    
print(max(Y))
Y.index(max(Y))












