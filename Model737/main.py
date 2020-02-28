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
Ca = 0.605
la = 2.661

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
xarray = np.append(xarray,la)
#Chordwise    
for i in np.arange(81):
    zi = -0.5*(Ca/2*(1-np.cos(thetaz[i]))+Ca/2*(1-np.cos(thetaz[i+1])))
    zarray = np.append(zarray,zi) 
zarray = np.append(zarray,-Ca) 


#Other Data
ha = 0.205 #m
z_hinge = ha/2 #m, z position of hinge line
theta = 28*np.pi/180 #radians
x1 = 0.172
x2 = 1.211
x3 = 2.591
xa = 0.35
xa1 = la/2 - xa/2
xa2 = la/2 + xa/2
y1 = 0.00681*cos(theta)
y2 = 0
y3 = 0.0203*cos(theta)
z1 = 0
z2 = 0
z3 = 0
a1 = 0
P = 97.4

I_yy = Iyy
I_zz = Izz
z_sc = 0.027624002859342803 + z_hinge
J = 8.275514338203897e-06
E = 73.1*10**6 #KPa
G = 28*10**6 #kPa

x, R1z, R2z, R3z, Ra1, R1y, R2y, R3y, C1, C2, C3, C4, C5 = symbols('x R1z R2z R3z Ra1 R1y R2y R3y C1 C2 C3 C4 C5', real=True)

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
dq = np.transpose(np.vstack((dq,dq)))
q0splines = d3splines(dq,xarray,0)
q1splines, q1 = intSpline(q0splines,xarray)
q2splines, q2 = intSpline(q1splines,xarray)
q3splines, q3 = intSpline(q2splines,xarray)
q4splines, q4 = intSpline(q3splines,xarray)

dm = np.transpose(np.vstack((dqm,dqm)))
m0splines = d3splines(dm,xarray,0)
m1splines, m1 = intSpline(m0splines,xarray)
m2splines, m2 = intSpline(m1splines,xarray)
m3splines, m3 = intSpline(m2splines,xarray)
m4splines, m4 = intSpline(m3splines,xarray) 

## Macauly Functions & Analysis
#
def Mac(x,i):
    if i == 0:
        if x < 0:
            return 0
        else:
            return 1
    else:
        return Heaviside(x)*x**i

def Sz(x):
    return R1z*Mac(x-x1,0) + R2z*Mac(x-x2,0) + R3z*Mac(x-x3,0) + cos(theta)*Ra1*Mac(x-xa1,0) - cos(theta)*P*Mac(x-xa2,0)

def Sy(x):
    base = -R1y*Mac(x-x1,0) - R2y*Mac(x-x2,0) - R3y*Mac(x-x3,0) - sin(theta)*Ra1*Mac(x-xa1,0) + sin(theta)*P*Mac(x-xa2,0)
    l = -1
    for j in range(len(xarray)-1):
        if max(xarray[j],x) != x:
            l = j-1
    temp = q1splines[l,0]*x + q1splines[l,1]
    return base+temp

def Mz(x):
    base = -R1y*Mac(x-x1,1) - R2y*Mac(x-x2,1) - R3y*Mac(x-x3,1) - sin(theta)*Ra1*Mac(x-xa1,1) + sin(theta)*P*Mac(x-xa2,1)
    l = -1
    for j in range(len(xarray)-1):
        if max(xarray[j],x) != x:
            l = j-1
    temp = q2splines[l,0]*x + q2splines[l,1]
    return base+temp

def My(x):
    return R1z*Mac(x-x1,1) + R2z*Mac(x-x2,1) + R3z*Mac(x-x3,1) + cos(theta)*Ra1*Mac(x-xa1,1) - cos(theta)*P*Mac(x-xa2,1)

def Torque(x):
    base = -ha/2*cos(theta)*P*Mac(x-xa2,0) + z_sc*sin(theta)*P*Mac(x-xa2,0) + ha/2*cos(theta)*Ra1*Mac(x-xa1,0) - z_sc*sin(theta)*Ra1*Mac(x-xa1,0) - (z_sc-z_hinge)*R1y*Mac(x-x1,0) - (z_sc-z_hinge)*R2y*Mac(x-x2,0) - (z_sc-z_hinge)*R3y*Mac(x-x3,0)
    l = -1
    for j in range(len(xarray)-1):
        if max(xarray[j],x) != x:
            l = j-1
    temp = m1splines[l,0]*x + m1splines[l,1]
    return base+temp

def twist(x):
    base = -ha/2*cos(theta)*P*Mac(x-xa2,1) + z_sc*sin(theta)*P*Mac(x-xa2,1) + ha/2*cos(theta)*Ra1*Mac(x-xa1,1) - z_sc*sin(theta)*Ra1*Mac(x-xa1,1) - (z_sc-z_hinge)*R1y*Mac(x-x1,1) - (z_sc-z_hinge)*R2y*Mac(x-x2,1) - (z_sc-z_hinge)*R3y*Mac(x-x3,1) + C5
    l = -1
    for j in range(len(xarray)-1):
        if max(xarray[j],x) != x:
            l = j-1
    temp = m2splines[l,0]*x + m2splines[l,1]
    return (base+temp)/(G*J)

def defly(x):
    base = (-1/(E*I_zz))*(-R1y/6*Mac(x-x1,3) - R2y/6*Mac(x-x2,3) - R3y/6*Mac(x-x3,3) - sin(theta)*Ra1/6*Mac(x-xa1,3) + sin(theta)*P/6*Mac(x-xa2,3) + C1*x + C2)
    l = -1
    for j in range(len(xarray)-1):
        if max(xarray[j],x) != x:
            l = j-1
    temp = q4splines[l,0]*x + q4splines[l,1]
    return base+temp

def deflz(x):
    return (-1/(E*I_yy))*(R1z/6*Mac(x-x1,3) + R2z/6*Mac(x-x2,3) + R3z/6*Mac(x-x3,3) + cos(theta)*Ra1/6*Mac(x-xa1,3) - cos(theta)*P/6*Mac(x-xa2,3) + C3*x + C4)

system = [Sz(la), 
       Sy(la), 
       Mz(la), 
       My(la), 
       Torque(la), 
       defly(x1) + twist(x1)*(z_sc-z_hinge) - y1,
       defly(x2) + twist(x2)*(z_sc-z_hinge) - y2,
       defly(x3) + twist(x3)*(z_sc-z_hinge) - y3,
       deflz(x1) - z1,
       deflz(x2) - z2,
       deflz(x3) - z3,
       deflz(xa1) + twist(xa1)*(ha/2) - a1]
sol = solve(system)
print(sol)  
# Values
C1 = -856.826048389909 
C2 = 1111.24391723489 
C3 = 0.143408535198873 
C4 = -0.0213678717446321 
C5 = 0.665268525843958 
R1y = 75.0037396172511 
R1z = -5.24585405391400 
R2y = -111.383785896061 
R2z = 2.47872225235871 
R3y = 34.4705541117598 
R3z = 8.07050893266608 
Ra1 = 31.9994523128544

xx = np.linspace(0,la,200)
y1 = np.array([])
y2 = np.array([])
for i in xx:
    y1 = np.append(y1, twist(i))    
    y2 = np.append(y2, Torque(i))    
y1 = y1*180/np.pi
plt.close()
plt.plot(xx,y1,'r')
plt.plot(xx,y2,'b')
plt.show()