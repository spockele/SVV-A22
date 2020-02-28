# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:21:48 2020

@author: daanv
"""

#modules
import numpy as np
import matplotlib.pyplot as plt

#toggles
plot = True 

#Data of the aileron cross section and stringers
Ca = 0.605 #[m]
ha = 0.205 #[m]
nstr = 13
tstr = 1.2e-3 #[m]
hstr = 1.6e-2 #[m]
wstr = 1.9e-2 #[m]
tsk = 1.1e-3 #[m]
tsp = 2.8e-3 #[m]

#Verification values
Izzver = 5.8159389575991465e-06
Iyyver = 4.363276766019503e-05
Czver = -0.19406263838748938
strver = np.array([[-0.        ,  0.        ],
       [-0.03725877,  0.07111646],
       [-0.11689227,  0.07988634],
       [-0.19847177,  0.06213382],
       [-0.28005126,  0.0443813 ],
       [-0.36163076,  0.02662878],
       [-0.44321025,  0.00887626],
       [-0.44321025, -0.00887626],
       [-0.36163076, -0.02662878],
       [-0.28005126, -0.0443813 ],
       [-0.19847177, -0.06213382],
       [-0.11689227, -0.07988634],
       [-0.03725877, -0.07111646]])
strver = np.transpose(strver)
Aver = 0.0021255886520793153


#Calculations on stringer placements
nstringhalf = int(np.ceil(nstr/2))
circlen = np.pi*ha/2 #[m]
straightlen = np.sqrt((Ca-ha/2)**2 + (ha/2)**2) #[m]
totlen = circlen + 2*straightlen #[m]
lenperstr = totlen/nstr  #[m]

ang = np.arctan((ha/2)/(Ca-ha/2)) 

#first stringer is located at the origin
zstr = [0]
ystr = [0]

for i in range(0,nstringhalf-2):
    zstr.append(Ca-np.cos(ang)*lenperstr*(1/2+i)) #[m]
    ystr.append(np.sin(ang)*lenperstr*(1/2+i)) #[m]
    #and for the other side
    zstr.append(Ca-np.cos(ang)*lenperstr*(1/2+i)) #[m]
    ystr.append(-np.sin(ang)*lenperstr*(1/2+i)) #[m]

#The final stringers are placed on the semicircular part 
strangle = lenperstr/(ha/2)
zstr.append(ha/2*(1-np.cos(strangle)))
ystr.append(ha/2*np.sin(strangle))
#and for the other side
zstr.append(ha/2*(1-np.cos(strangle)))
ystr.append(-ha/2*np.sin(strangle))
 
#Calculating stringer and skin areas
Astr = wstr*tstr+hstr*tstr #[m^2]
Askcirc = circlen*tsk #[m^2]
Askstraight = straightlen*tsk #[m^2]
Asp = tsp*ha #[m^2]

#for a semicircle with radius R, its centroid lies at 2R/pi from its origin
#i.e. for this semicircle, at 2*(ha/2)/pi
zcirc = ha/2 - ha/np.pi #[m]

#for a straight line, its centroid lies halfway
zstraight = Ca - np.cos(ang)*straightlen/2 #[m]
ystraight = np.sin(ang)*straightlen/2 #[m]

zsp = ha/2

Atot = nstr*Astr + Askcirc + 2*Askstraight + Asp

#Switching to the correct coordinate system:
zstr = np.array(zstr)*-1
zcirc = zcirc*-1
zsp = zsp*-1
zstraight = zstraight*-1


#Calculating the centroid of the aileron
#In y-direction, it is 0 due to symmetry
Cz = (sum(zstr)*Astr + zcirc*Askcirc + 2*zstraight*Askstraight + zsp*Asp)/Atot #[m]


#For a thin walled semicircle, the moment of Inertia about its own centroid is pi/2*R^3*t
Icirc = np.pi/2*(ha/2)**3*tsk #[m^4]
#About the Y axis: we have to subtract a Steiner term
Icirc0 = Icirc - Askcirc*(zsp-zcirc)**2

#For a thin walled plate at an angle, its moment of Inertia about its own centroid is t*L^3*sin^2(a)/12 or t*L^3*cos^2(a)/12
#About Z axis: sin^2
#About Y axis: cos^2
Istraightz = tsk*straightlen**3*(np.sin(ang))**2/12 #[m^4]
Istraighty = tsk*straightlen**3*(np.cos(ang))**2/12 #[m^4]

Ispz = tsp*ha**3/12
#Ispy = 0 due to thin walled assumption

#For the total MoI, sum the individual parts and their corresponding Steiner Terms
#For the stringers, only Steiner terms will be considered

Izz = Icirc + 2*(Istraightz + Askstraight*ystraight**2) + sum((np.array(ystr))**2*Astr) + Ispz#[m^4]
Iyy = Icirc + Askcirc*(zcirc-Cz)**2 + 2*(Istraighty + Askstraight*(zstraight-Cz)**2) + sum((np.array(zstr-Cz))**2*Astr) + Asp*(zsp-Cz)**2 #[m^4]

if plot:
    xcirc = []
    ycirc = []
    for i in np.arange(-np.pi/2,np.pi/2,0.1):
        xcirc.append(np.cos(i))
        ycirc.append(np.sin(i))
    xcirc = np.array(xcirc)*ha/2- ha/2
    ycirc = np.array(ycirc)*ha/2
    xcirc = np.insert(xcirc,0, -Ca)
    xcirc = np.append(xcirc, -Ca)
    ycirc = np.insert(ycirc,0, 0)
    ycirc = np.append(ycirc, 0)
    fig = plt.gca()
    fig.set_aspect('equal')    
    plt.plot(xcirc,ycirc)
    plt.plot([zsp,zsp],[ha/2,-ha/2])
    plt.scatter(zstr,ystr)
    plt.scatter(zcirc,0)
    plt.scatter([zstraight,zstraight],[ystraight,-ystraight])
    plt.scatter(Cz,0)
    plt.scatter(zsp, 0)
    plt.show()
print(straightlen)

