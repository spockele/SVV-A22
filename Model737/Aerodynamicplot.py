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
from integration import Polynomial
from interpolation import Interpolate

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

# Unit test for constant data
#datalist = np.ones((81,43))
#zarray = np.linspace(0,Ca,43)
#xarray = np.linspace(0,la,83)


# Linear splines & integrations for each cross section
dq = np.array([])
cop = np.array([])

max_h = 0
for i in range(len(zarray)-1):
    hi = abs(zarray[i+1]-zarray[i])
    if hi > max_h:
        max_h = hi

#for x in range(81):
#    si = np.array([0, 0])
#    for z in range(42):
#        fi0 = datalist[x,z]
#        fi1 = datalist[x,z+1]
#        zi0 = zarray[z]
#        zi1 = zarray[z+1]
#        spline = np.array([(fi1 - fi0)/(zi1 - zi0), fi0 - (fi1 - fi0)/(zi1 - zi0)*(zi0)]);
#        # splines are formatted as a_n, a_n-1, ... , a_0
#        si = np.vstack((si, spline))
#    si = np.delete(si, 0, 0)
#    
#    #all the splines for the cross section are in this array now we integrate
#    integral = 0
#    weighted = 0
#    for z in range(len(zarray)-1):
#        sp = Polynomial(si[z,0], si[z,1])
#        sp = sp.definite_integral(zarray[z], zarray[z+1], 1)
#        integral += sp
#        spw = Polynomial(si[z,0], si[z,1], 0)
#        spw = spw.definite_integral( zarray[z], zarray[z+1], 1)
#        weighted += spw
#    dq = np.append(dq, integral)
#    cop = np.append(cop, weighted/integral)


   
zarray = np.arange(0,10.01,0.01)
datalist = [zarray**3]
datalist = np.array(datalist)
datalist = np.tile(datalist, (2,1))
datalist = np.transpose(datalist)
io = Interpolate()
io = io.d3spline(datalist, zarray, 1)

p = Polynomial(1,0,0,0)
p_int = p.definite_integral(0,10, 7)
print(p_int)