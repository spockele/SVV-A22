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

#setting up the plot
fig = plt.figure()
ax = fig.gca(projection='3d')

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

#In this list (matrix), the rows correspond to the z-coordinate in the coordinate system
#and the columns correspond to the x-coordinate in the coordinate system

#calculating the z and x coordinates:
Nz = len(datalist)
Nx = len(datalist[0])
Ca = 0.484
la = 1.691

thetax = np.array([])
thetaz = np.array([])

for i in np.arange(1,Nx+2):
    thetax = np.append(thetax,np.pi*(i-1)/Nx)
    
for i in np.arange(1,Nz+2):
    thetaz = np.append(thetaz,np.pi*(i-1)/Nz)

zarray = np.array([])
xarray = np.array([])

#Spanwise
for i in np.arange(0,Nx):
    xi = 0.5*(la/2*(1-np.cos(thetax[i]))+la/2*(1-np.cos(thetax[i+1])))
    xarray = np.append(xarray,xi)

#Chordwise    
for i in np.arange(0,Nz):
    zi = -0.5*(Ca/2*(1-np.cos(thetaz[i]))+Ca/2*(1-np.cos(thetaz[i+1])))
    zarray = np.append(zarray,zi) 

#Interpolation
for z in range(81):
    #Complete spline derivative conditions
    df0 = -(datalist[z,1] - datalist[z,0])/(zarray[1] - zarray[0]);
    dfn = -(datalist[z,40] - datalist[z,39])/(zarray[40] - zarray[39]);
    
    # contructing the linalg equation Az=d to solve for z
    d = np.array([]);
    h = np.array([]);
    A = np.zeros((40,40));
    s = np.zeros((41, 4))
    # finding the step sizes and making A
    for x in range(40):
            h = np.append(h, zarray[x+1]-zarray[x]);
            A[x:,x:] = h[x];
    A[0,0] = 2*h[0];
    A[39,39] = 2*h[39];
    # A note about h and A (which is constructed from h values). h is the step, so for
    # n given nodes/knots there are n-1 steps, and so h[0] is the first step, from x0 to x1,
    # while h[39] is the final step, from x39 to x40.
    
    # finding d
    for x in range(41):
        if x == 0:
            d = np.append(d, 6/h[x]*(datalist[z,1]-datalist[z,0] - 6*df0));
        else:
            if x == 40:
                d = np.append(d, 6*dfn - 6/h[x-1]*(datalist[z,40]-datalist[z,39]));
            else:
                d = np.append(d, 6/h[x]*(datalist[z,x+1]-datalist[z,x] - 6/h[x-1]*(datalist[z,x]-datalist[z,x-1])));
                A[x-1,x-1] = 2*(h[x-1] + h[x]);
    # solving the equation
    Mi = np.linalg.solve(A,d);
    for x in range(40):
        
    

#Calculating the mean pressure location
dataT = np.transpose(datalist)
cavg = np.array([])
for i in dataT:
    cavg = np.append(cavg,sum(i*zarray)/sum(i))

#avgline = ax.plot(xarray,cavg,15)
