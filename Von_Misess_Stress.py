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
C1 = 0.279519165731257
C2 = -0.0384445931088122
C3 = -0.417965528308128
C4 = 0.0619789183565207
C5 = 4587.38387648207
R1y = -3.24159865890804
R1z = -12.3218865035638
R2y = 5.00062052855391
R2z = 24.8263794419833
R3y = -0.607282578573421
R3z = -0.924661188757137
Ra1 = 25.0162567249494


#stresses in cross section
for x in range(0,int(l_a/dx)):
    sy = Sy(x)
    for i in range(len(ystr)):
        sigma_x.append(Mz(x)/Izz*ystr[i] + (My(x))/Iyy*zstr[i])
        sigma_y.append((Mz(x))/Izz*x)
        sigma_z.append((My(x))/Iyy*x)
        tau_yz_lst.append(tau_yz(i, sy))

for k in range(0,len(sigma_x)):
    #Y.append((0.5*((sigma_x[k]-sigma_y[k])**2+(sigma_y[k]-sigma_z[k])**2+(sigma_z[k]-sigma_x[k])**2)+3*tau_yz_lst[k]**2)**0.5)
    Y.append((3*tau_yz_lst[k]**2)**0.5)

print(max(Y))   



#print(max(Y))
#Y.index(max(Y))



