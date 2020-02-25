# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 13:59:34 2020

@author: Daan, Mohammad
"""
#importing modules
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from integration import Polynomial
from interpolation import *
from math import cos, sin
from sympy import *
from Moment_of_Inertia import Iyy, Izz

############################# Data Setup ###################################
#reading and obtaining the data points
data = open('aerodynamicloadcrj700.dat','r')
datalist = []
for i in data:
    j = i.split(",")
    datarow = np.array([])
    for k in j:
        k = float(k)
        datarow = np.append(datarow,k)
    datalist.append(datarow)

datalist = np.array(datalist)
datalist = np.pad(datalist, 1, 'constant')
#In this list (matrix), the rows correspond to the z-coordinate in the coordinate system
#and the columns correspond to the x-coordinate in the coordinate system

#calculating the z and x coordinates:
Nx = len(datalist[0])
Nz = len(datalist)
Ca = 0.484
la = 1.691

thetax = np.array([])
thetaz = np.array([])

for i in np.arange(1,43):
    thetax = np.append(thetax,np.pi*(i-1)/41)
    
for i in np.arange(1,83):
    thetaz = np.append(thetaz,np.pi*(i-1)/81)

zarray = np.array([0])
xarray = np.array([0])
#Spanwise
for i in np.arange(41):
    xi = 0.5*(la/2*(1-np.cos(thetax[i]))+la/2*(1-np.cos(thetax[i+1])))
    xarray = np.append(xarray,xi)
xarray = np.append(xarray,-la)
#Chordwise    
for i in np.arange(81):
    zi = -0.5*(Ca/2*(1-np.cos(thetaz[i]))+Ca/2*(1-np.cos(thetaz[i+1])))
    zarray = np.append(zarray,zi) 
zarray = np.append(zarray,-Ca) 


#Other Data
ha = 0.173 #m
z_hinge = ha/2 #m, z position of hinge line
theta = 26*np.pi/180 #radians
x1 = 0.149
x2 = 0.554
x3 = 1.541
xa = 0.272
xa1 = Ca/2 - xa/2
xa2 = Ca/2 + xa/2
y1 = 0.00681*cos(theta)
y2 = 0
y3 = 0.0203*cos(theta)
z1 = 0
z2 = 0
z3 = 0
a1 = 0
P = 37.9

I_yy = Iyy
I_zz = Izz
z_sc = -0.027624002859342803 + ha
J = 8.275514338203897e-06
E = 1
G = 1

x, R1z, R2z, R3z, Ra1, R1y, R2y, R3y, C1, C2, C3, C4, C5 = symbols('x, R1z, R2z, R3z, Ra1, R1y, R2y, R3y, C1, C2, C3, C4, C5', real=True)

####################### Calculating the Distribution #########################
# Unit test for correct center of pressure
#datalist = np.ones((83,43))
#xarray = np.arange(0,Ca,Ca/43)
#zarray = np.arange(0,la,la/83)

dq = np.array([])
dqm = np.array([])
cop = np.array([])
for i in range(Nx):
    splines = d3splines(datalist, zarray, i)
    dq_temp, dqm_temp, cop_temp = resultantload(splines, zarray, z_sc)
    dq = np.append(dq, dq_temp)
    dqm = np.append(dqm, dqm_temp)
    cop = np.append(cop, cop_temp)
# Find polynomial describibing the loading
#print(dq)
#print(dqm)
#print(cop)


#lagrangePoly = lagrange(xarray,dqm)
#print(lagrangePoly)

#As calculated, we have that the loading along the Span is as follows:
Loading = tuple([-29509481109147.1, 973089896901375., -1.53534963623188e+16, 1.54143203855473e+17, -1.10423627336087e+18, 5.99623419445990e+18, -2.55711908765841e+19, 8.74175881223672e+19, -2.41726851774350e+20, 5.38505232159416e+20, -9.40780774833795e+20, 1.17417598503425e+21, -6.14999921565349e+20, -1.58790426834909e+21, 6.03895442263498e+21, -1.24603828547045e+22, 1.93975666854832e+22, -2.46837715938424e+22, 2.65182517394726e+22, -2.44386699684457e+22, 1.94870693141255e+22, -1.35070631821410e+22, 8.15551918173636e+21, -4.29130301166523e+21, 1.96583402334124e+21, -7.82313899450026e+20, 2.69564466254067e+20, -8.00707322662478e+19, 2.03880947538399e+19, -4.41931008930020e+18, 8.08578364084700e+17, -1.23594063750417e+17, 1.55855421427658e+16, -1.59663105896230e+15, 130364565874585., -8282687446888.95, 397088926359.289, -13796246243.3501, 328870936.981902, -4983645.55435716, 43270.9107925202, -211.492182015877, 0.0])
Moment = tuple([1248505240737.82, -41170141318112.8, 649586894027260., -6.52161095115704e+15, 4.67189454277131e+16, -2.53693943961908e+17, 1.08188957938163e+18, -3.69854798589831e+18, 1.02272243004717e+19, -2.27836485869779e+19, 3.98035960391797e+19, -4.96784320570175e+19, 2.60203555039730e+19, 6.71823822931072e+19, -2.55502293324393e+20, 5.27187407275235e+20, -8.20693949469037e+20, 1.04434917890153e+21, -1.12196490268832e+21, 1.03397985692282e+21, -8.24481922078627e+20, 5.71472855727625e+20, -3.45053416578760e+20, 1.81561566280006e+20, -8.31728505801697e+19, 3.30990679827259e+19, -1.14050533232560e+19, 3.38772705247185e+18, -8.62603394372117e+17, 1.86977292921058e+17, -3.42102567736621e+16, 5.22915589960339e+15, -659410045361003., 67551928253251.3, -5515591233911.65, 350431239696.032, -16800332333.1533, 583698740.281493, -13913912.4919012, 210844.072085794, -1830.43627502361, 8.95028196636646, 0.0])
for i in range(len(Loading)):
    print(Loading[i]*x**(len(Loading)-i-1))

q = np.array([]) # Loading integrals
for i in range(len(Loading)):
    q = np.append(q,Loading[i])
q = Polynomial(q)
qqq = np.array(list(q.definite_integral(0,x,1)))
q1 = np.array([])
for i in range(len(qqq)):
    q1 = np.append(q1,qqq[i]*x**(len(qqq)-i-1))
qqq = np.array(list(q.definite_integral(0,x,2)))
q2 = np.array([])
for i in range(len(qqq)):
    q2 = np.append(q2,qqq[i]*x**(len(qqq)-i-2))    
qqq = np.array(list(q.definite_integral(0,x,4)))
q4 = np.array([])
for i in range(len(qqq)):
    q4 = np.append(q4,qqq[i]*x**(len(qqq)-i-4)) 
    
qm = np.array([]) # Moment integrals
for i in range(len(Loading)):
    qm = np.append(qm,Loading[i])
qm = Polynomial(qm)
mmm = np.array(list(qm.definite_integral(0,x,1)))
qm1 = np.array([])
for i in range(len(mmm)):
    qm1 = np.append(qm1,mmm[i]*x**(len(qqq)-i-1))
mmm = np.array(list(qm.definite_integral(0,x,2)))
qm2 = np.array([])
for i in range(len(mmm)):
    qm2 = np.append(qm2,mmm[i]*x**(len(qqq)-i-2))    

    
    
    
# Macauly Functions & Analysis

def Mac(x,i):
    if i == 0:
        return 1
    else:
        return Heaviside(x)*x**i

Sz = R1z*Mac(x-x1,0) + R2z*Mac(x-x2,0) + R3z*Mac(x-x3,0) + cos(theta)*Ra1*Mac(x-xa1,0) - cos(theta)*P*Mac(x-xa2,0)

Sy = -R1y*Mac(x-x1,0) - R2y*Mac(x-x2,0) - R3y*Mac(x-x3,0) - sin(theta)*Ra1*Mac(x-xa1,0) + sin(theta)*P*Mac(x-xa2,0) 
for i in q1:
    Sy = Sy+i

Mz = -R1y*Mac(x-x1,1) - R2y*Mac(x-x2,1) - R3y*Mac(x-x3,1) - sin(theta)*Ra1*Mac(x-xa1,1) + sin(theta)*P*Mac(x-xa2,1) 
for i in q2:
    Mz = Mz + i

My = R1z*Mac(x-x1,1) + R2z*Mac(x-x2,1) + R3z*Mac(x-x3,1) + cos(theta)*Ra1*Mac(x-xa1,1) - cos(theta)*P*Mac(x-xa2,1)

Torque = -ha/2*cos(theta)*P*Mac(x-xa2,0) + z_sc*sin(theta)*P*Mac(x-xa2,0) + ha/2*cos(theta)*Ra1*Mac(x-xa1,0) - z_sc*sin(theta)*Ra1*Mac(x-xa1,0) - (z_sc-z_hinge)*R1y*Mac(x-x1,0) - (z_sc-z_hinge)*R2y*Mac(x-x2,0) - (z_sc-z_hinge)*R3y*Mac(x-x3,0)
for i in range(len(q1)):
    Torque = Torque + qm1[i]
    
twist = -ha/2*cos(theta)*P*Mac(x-xa2,1) + z_sc*sin(theta)*P*Mac(x-xa2,1) + ha/2*cos(theta)*Ra1*Mac(x-xa1,1) - z_sc*sin(theta)*Ra1*Mac(x-xa1,1) - (z_sc-z_hinge)*R1y*Mac(x-x1,1) - (z_sc-z_hinge)*R2y*Mac(x-x2,1) - (z_sc-z_hinge)*R3y*Mac(x-x3,1) + C5
for i in range(len(q2)):
    twist = twist + q2[i] + qm2[i]
twist = twist/(G*J)

defly = -R1y/6*Mac(x-x1,3) - R2y/6*Mac(x-x2,3) - R3y/6*Mac(x-x3,3) - sin(theta)*Ra1/6*Mac(x-xa1,3) + sin(theta)*P/6*Mac(x-xa2,3) + C1*x + C2
for i in q4:
    defly = defly + i
defly = defly*(-1/(E*I_zz))

deflz = R1z/6*Mac(x-x1,3) + R2z/6*Mac(x-x2,3) + R3z/6*Mac(x-x3,3) + cos(theta)*Ra1/6*Mac(x-xa1,3) - cos(theta)*P/6*Mac(x-xa2,3) + C3*x + C4
deflz = deflz*(-1/(E*I_yy))

Sz.subs(x,la)
# Solving the equations
system = [Sz.subs(x,la), 
       Sy.subs(x,la), 
       Mz.subs(x,la), 
       My.subs(x,la), 
       Torque.subs(x,la), 
       defly.subs(x,x1) + twist.subs(x,x1)*(z_sc-z_hinge) - y1,
       defly.subs(x,x2) + twist.subs(x,x2)*(z_sc-z_hinge) - y2,
       defly.subs(x,x3) + twist.subs(x,x3)*(z_sc-z_hinge) - y3,
       deflz.subs(x,x1) - z1,
       deflz.subs(x,x2) - z2,
       deflz.subs(x,x3) - z3,
       deflz.subs(x,xa1) + twist.subs(x,xa1)*(ha/2) - a1]
print(nonlinsolve(system, [R1z, R2z, R3z, Ra1, R1y, R2y, R3y, C1, C2, C3, C4, C5]))
print(system)