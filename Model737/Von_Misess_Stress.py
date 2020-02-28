# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 06:19:17 2020

@author: HJ Hoogendoorn

Maximum stress calculations using Von 
mises stress
"""

l_a = 2.661 #[m] span of aileron
from Moment_of_Inertia import *
from main import *
import numpy as np
from scratch_1 import q_b_LC, q_b_RC, t_skin
import base_shearflow as qb
import math



#integration steps

dx = 0.001 
dy = 0.001
dz = 0.001


sigma_x = [] 
sigma_y = []
sigma_z = []
tau_yz_lst = []


def tau_yz(i_str, sy):
    if i_str == 0:
        tau = q_b_LC / t_skin

    elif i_str in (1, 3, 5, 7, 9):
        tau = (qb.qb3(sy, qb.s3_str(zstr[i_str], ystr[i_str])) + q_b_RC) / t_skin

    elif i_str in (2, 4, 6, 8, 10):
        tau = (qb.qb4(sy, qb.s4_str(zstr[i_str], ystr[i_str])) + q_b_RC) / t_skin

    elif i_str == 11:
        tau = (qb.qb1(sy, math.atan(ystr[i_str])/(qb.h - zstr[i_str])) + q_b_LC) / t_skin

    elif i_str == 12:
        tau = (qb.qb6(sy, math.atan(-ystr[i_str]) / (qb.h - zstr[i_str])) + q_b_LC) / t_skin

    else:
        raise ValueError("No stringer with this index exists")

    return tau


    

Y = []
C1 = -7198.75274179719 
C2 = 14607.6260287575 
C3 = 2.27188428673498 
C4 = -0.390764097318417 
C5 = 0.319684787967052
R1y = 12.4791795238527
R1z = -12.6407286520606 
R2y = -43.2303087688269 
R2z = -2.66458088954677
R3y = 20.3679700072242
R3z =  12.4105963431757
Ra1 = 100.678465473878


#stresses in cross section
for x in range(0,int(l_a/dx)):
    sy = Sy(x)
    for i in range(len(ystr)):
        sigma_x.append(Mz(x)/Izz*ystr[i] + (My(x))/Iyy*zstr[i])
        sigma_y.append((Mz(x))/Izz*x)
        sigma_z.append((My(x))/Iyy*x)
        tau_yz_lst.append(tau_yz(i, sy))

for k in range(0,len(sigma_x)):
    Y.append((0.5*((sigma_x[k]-sigma_y[k])**2+(sigma_y[k]-sigma_z[k])**2+(sigma_z[k]-sigma_x[k])**2)+3*tau_yz_lst[k]**2)**0.5)

print(max(Y))   



#print(max(Y))
#Y.index(max(Y))



