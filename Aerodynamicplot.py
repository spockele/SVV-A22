# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:16:04 2020

@author: daanv
"""

#importing modules
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from integration import *


#reading and obtaining the data points
data = open('aerodynamicloadcrj700.dat','r')
datalist = []
for i in data:
    j = i.split(",")
    datarow = np.array([0])
    for k in j:
        k = float(k)
        datarow = np.append(datarow,k)
    datarow = np.append(datarow,0)
    datalist.append(datarow)

datalist = np.array(datalist)
#In this list (matrix), the rows correspond to the z-coordinate in the coordinate system
#and the columns correspond to the x-coordinate in the coordinate system

#calculating the z and x coordinates:
Nz = len(datalist[0])
Nx = len(datalist)
Ca = 0.484
la = 1.691

thetax = np.array([])
thetaz = np.array([])

for i in np.arange(1,83):
    thetax = np.append(thetax,np.pi*(i-1)/Nx)
    
for i in np.arange(1,43):
    thetaz = np.append(thetaz,np.pi*(i-1)/Nz)

zarray = np.array([0])
xarray = np.array([0])

#Spanwise
for i in np.arange(0,81):
    xi = 0.5*(la/2*(1-np.cos(thetax[i]))+la/2*(1-np.cos(thetax[i+1])))
    xarray = np.append(xarray,xi)
xarray = np.append(xarray,-la)

#Chordwise    
for i in np.arange(0,41):
    zi = -0.5*(Ca/2*(1-np.cos(thetaz[i]))+Ca/2*(1-np.cos(thetaz[i+1])))
    zarray = np.append(zarray,zi) 
zarray = np.append(zarray,-Ca) 

# Linear splines & integrations for each cross section
dq = np.array([])
cop = np.array([])

for x in range(81):
    si = np.array([0, 0])
    for z in range(42):
        fi0 = datalist[x,z]
        fi1 = datalist[x,z+1]
        zi0 = zarray[z]
        zi1 = zarray[z+1]
        spline = np.array([(fi1 - fi0)/(zi1 - zi0), fi0 - (fi1 - fi0)/(zi1 - zi0)*(zi0)]);
        # splines are formatted as a_n, a_n-1, ... , a_0
        si = np.vstack((si, spline))
    si = np.delete(si, 0, 0)
    #all the splines for the cross section are in this array now we integrate
    integral = 0
    moment = 0
    for z in range(len(zarray)-1):
        sp = Polynomial(si[z,0], si[z,1])
        sp = sp.definite_integral(zarray[z], zarray[z+1], 1)
        z_arm = (zarray[z] - zarray[z+1])/2 + zarray[z]
        integral += sp
        moment += sp*z_arm
    dq = np.append(dq, integral)
    cop = np.append(cop, moment/integral)

print(cop)
print(len(dq))   
    
#avgline = ax.plot(xarray,cavg,15)
