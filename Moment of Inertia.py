# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:21:48 2020

@author: daanv
"""

#modules
import numpy as np
import matplotlib.pyplot as plt

#Data of the aileron cross section and stringers
Ca = 1.691 #[m]
ha = 0.173 #[m]
nstr = 13
tstr = 1.2e-3 #[m]
hstr = 1.4e-2 #[m]
wstr = 1.8e-2 #[m]
tsk = 1.1e-3 #[m]

#Calculations on stringer placements
#Only done on one side due to symmetry
nstringhalf = int(np.ceil(nstr/2))
circlen = np.pi*ha/2 #[m]
halfcirclen = circlen/2
straightlen = np.sqrt((Ca-ha/2)**2 + (ha/2)**2) #[m]
totlen = halfcirclen + straightlen #[m]
lenperstr = totlen/nstringhalf  #[m]
#The distance between the stringers is larger than the length of the circle part of the aileron
#i.e. except for the stringer at the leading edge, all stringers are placed on the straigth part of the aileron
ang = np.arctan((ha/2)/(Ca-ha/2)) #[rad]
zstr = [0]
ystr = [0]

for i in range(1,nstringhalf):
    zstr.append(Ca-np.cos(ang)*i*lenperstr) #[m]
    ystr.append(np.sin(ang)*i*lenperstr) #[m]
    #and for the other side
    zstr.append(Ca-np.cos(ang)*i*lenperstr) #[m]
    ystr.append(-np.sin(ang)*i*lenperstr) #[m]

    
#Calculating stringer and skin areas
Astr = wstr*tstr+(hstr-tstr)*tstr #[m^2]
Askcirc = circlen*tsk
Askstraight = straightlen*tsk

#for a halfcircle with radius R, its centroid lies at 2R/pi from its origin
#i.e. for this quartercircle, at 2*(ha/2)/pi
zcirc = ha - ha/np.pi

#for a straigt line, its centroid lies halfway
zstraight = Ca - np.cos(ang)*straightlen/2

#Calculating the centroid of the aileron
#In y-direction, it is 0 due to symmetry
Cz = (sum(zstr)*Astr + zcirc*Askcirc + zstraight*Askstraight)/(nstr*Astr + Askcirc + Askstraight)



