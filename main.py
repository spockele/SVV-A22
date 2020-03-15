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

for i in np.arange(1,Nx):
    thetax = np.append(thetax,np.pi*(i-1)/41)
    
for i in np.arange(1,Nz):
    thetaz = np.append(thetaz,np.pi*(i-1)/81)

zarray = np.array([0])
xarray = np.array([0])
#Spanwise
for i in np.arange(Nx-2):
    xi = 0.5*(la/2*(1-np.cos(thetax[i]))+la/2*(1-np.cos(thetax[i+1])))
    xarray = np.append(xarray,xi)
xarray = np.append(xarray,la)
#Chordwise    
for i in np.arange(Nz-2):
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
xa1 = la/2 - xa/2
xa2 = la/2 + xa/2
d1 = 0.00681
d3 = 0.0203
y1 = d1*cos(theta)
y2 = 0
y3 = d3*cos(theta)
z1 = d1*sin(theta)
z2 = 0
z3 = d3*sin(theta)
a1 = 0
P = 37.9

I_yy = Iyy
I_zz = Izz
z_sc = -0.027624002859342803 + z_hinge
J = 8.275514338203897e-06
E = 73.1*10**6 #KPa
G = 20*10**6 #kPa

x, R1z, R2z, R3z, Ra1, R1y, R2y, R3y, C1, C2, C3, C4, C5 = symbols('x R1z R2z R3z Ra1 R1y R2y R3y C1 C2 C3 C4 C5', real=True)

####################### Calculating the Distribution #########################
# Unit test for correct center of pressure
#datalist = np.ones((83,43))
#xarray = np.arange(0,la,la/43)
#zarray = np.arange(0,Ca,Ca/83)

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
    return -(- R1z*Mac(x-x1,0) - R2z*Mac(x-x2,0) - R3z*Mac(x-x3,0) + cos(theta)*Ra1*Mac(x-xa1,0) - cos(theta)*P*Mac(x-xa2,0))

def Sy(x):
    base = - R1y*Mac(x-x1,0) - R2y*Mac(x-x2,0) - R3y*Mac(x-x3,0) - sin(theta)*Ra1*Mac(x-xa1,0) + sin(theta)*P*Mac(x-xa2,0)
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
    return (base+temp)

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
    base = -ha/2*cos(theta)*P*Mac(x-xa2,1) + z_sc*sin(theta)*P*Mac(x-xa2,1) + ha/2*cos(theta)*Ra1*Mac(x-xa1,1) - z_sc*sin(theta)*Ra1*Mac(x-xa1,1) - (z_sc-z_hinge)*R1y*Mac(x-x1,1) - (z_sc-z_hinge)*R2y*Mac(x-x2,1) - (z_sc-z_hinge)*R3y*Mac(x-x3,1) + C5*G*J
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
       deflz(xa1)*cos(theta) + defly(xa1)*sin(theta) + twist(xa1)*(ha/2) - a1]
sol = solve(system)
print(sol)  
# Values
C1 = -856.756808300039; C2 = 1111.23123884183; C3 = -31.7892338942513; C4 = 14.4562135386494; C5 = -0.00485172476832893; R1y = 77.0376772716926; R1z = 284.963068074844; R2y = -120.913937040834; R2z = -423.787107363642; R3y = 34.0593497248714; R3z = 149.733271408323; Ra1 = 50.0376328252782

xx = np.linspace(0,la,200)
dy = np.array([])
dz = np.array([])
sy = np.array([])
sz = np.array([])
my = np.array([])
mz = np.array([])
torq = np.array([])
tw = np.array([])
for i in xx:
    dy = np.append(dy, defly(i))   
    dz = np.append(dz, deflz(i))    
    sy = np.append(sy, Sy(i)*1000)    
    sz = np.append(sz, Sz(i)*1000)    
    my = np.append(my, My(i)*1000)    
    mz = np.append(mz, Mz(i)*1000)    
    torq = np.append(torq, Torque(i)*1000)    
    tw = np.append(tw, twist(i))     
plt.close()
### Deflections ###
plt.figure(0)
plt.xlabel('x [m]')
plt.ylabel('deflection(x) [m]')
plt.title('Deflection along the aileron span')
plt.plot(xx,dy,'r',label='deflection in y')
plt.plot(xx,dz,'b',label='deflection in z')
plt.legend()
plt.show()
### Shears ###
plt.figure(1)
plt.xlabel('x [m]')
plt.ylabel('Shear(x) [kN]')
plt.title('Shear force along the aileron span')
plt.plot(xx,sy,'r',label='shear in y')
plt.plot(xx,sz,'b',label='shear in z')
plt.legend()
plt.show()
### Moments ###
plt.figure(2)
plt.xlabel('x [m]')
plt.ylabel('Moment(x) [kNm]')
plt.title('Moment along the aileron span')
plt.plot(xx,my,'r',label='moment in y')
plt.plot(xx,mz,'b',label='moment in z')
plt.legend()
plt.show()
### Torque & Twist ###
plt.figure(3)
plt.xlabel('x [m]')
plt.ylabel('Torque(x) [kNm]')
plt.title('Torque along the aileron span')
plt.plot(xx,torq,'r',label='Torque')
plt.legend()
plt.show()
### Torque & Twist ###
plt.figure(4)
plt.xlabel('x [m]')
plt.ylabel('Twist(x) [degrees]')
plt.title('Twist along the aileron span')
plt.plot(xx,tw,'b',label='Twist')
plt.legend()
plt.show()