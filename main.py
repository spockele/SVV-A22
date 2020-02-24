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

####################### Calculating the Distribution #########################
# Unit test for correct center of pressure
#datalist = np.ones((83,43))
#xarray = np.arange(0,Ca,Ca/43)
#zarray = np.arange(0,la,la/83)

dq = np.array([])
dqm = np.array([])
cop = np.array([])
#for i in range(Nx):
#    splines = d3splines(datalist, zarray, i)
#    dq_temp, dqm_temp, cop_temp = resultantload(splines, zarray, z_hinge)
#    dq = np.append(dq, dq_temp)
#    dqm = np.append(dqm, dqm_temp)
#    cop = np.append(cop, cop_temp)
## Find p[ol]
print(dq)
print(dqm)
print(cop)


#lagrangePoly = lagrange(xarray,dq)
#print(lagrangePoly)

Loading = [-29509481109147.1, 973089896901375., -1.53534963623188e+16, 1.54143203855473e+17, -1.10423627336087e+18, 5.99623419445990e+18, -2.55711908765841e+19, 8.74175881223672e+19, -2.41726851774350e+20, 5.38505232159416e+20, -9.40780774833795e+20, 1.17417598503425e+21, -6.14999921565349e+20, -1.58790426834909e+21, 6.03895442263498e+21, -1.24603828547045e+22, 1.93975666854832e+22, -2.46837715938424e+22, 2.65182517394726e+22, -2.44386699684457e+22, 1.94870693141255e+22, -1.35070631821410e+22, 8.15551918173636e+21, -4.29130301166523e+21, 1.96583402334124e+21, -7.82313899450026e+20, 2.69564466254067e+20, -8.00707322662478e+19, 2.03880947538399e+19, -4.41931008930020e+18, 8.08578364084700e+17, -1.23594063750417e+17, 1.55855421427658e+16, -1.59663105896230e+15, 130364565874585., -8282687446888.95, 397088926359.289, -13796246243.3501, 328870936.981902, -4983645.55435716, 43270.9107925202, -211.492182015877, 0.0]
print(Loading)