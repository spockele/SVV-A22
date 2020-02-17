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

#Not allowed in final code, but good enough for quick calculations
from scipy.integrate import simps

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

#Plotting the graph
x,y = np.meshgrid(xarray,zarray)
ax.set_xlabel('Spanwise [m]')
ax.set_ylabel('Chordwise [m]')
ax.set_zlabel('Aerodynamic load [kN/m²]')

surf = ax.plot_surface(x, y, datalist, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

#Calculating the mean pressure location
dataT = np.transpose(datalist)
cavg = np.array([])
for i in dataT:
    cavg = np.append(cavg,sum(i*zarray)/sum(i))

#avgline = ax.plot(xarray,cavg,15)