#we find shear centre
#input data
import math as m
import numpy as np
from Moment_of_Inertia import Izz
Sy=1
t_skin=1.1/1000   #m
t_spar=2.8/1000   #m
h=0.205      #m
r=h/2
ca=0.605     #m
lsk=((ca-r)**2+(r**2))**0.5
alpha=m.acos((ca-r)/lsk)
l=r*m.sin(alpha) #moment arm for shear force for slanting bit
#converting the equation into matrix form
q_b_LC=((Sy*r**3)/Izz)*(1-2/6+1-m.pi/2)
q_b_RC=(-Sy/Izz)*((r**3/3)+(r*lsk**2/6)+(r**2*t_spar*lsk/t_skin))
x1=((2*r*m.pi)/t_skin+(2*r)/t_spar)
y1=(-2*r/t_spar)
x2=(-2*r/t_spar)
y2=(2*lsk/t_skin+(2*r/t_spar))
A=np.array([[x1,y1],[x2,y2]])
B=np.array([q_b_LC,q_b_RC])
C=np.linalg.solve(A,B)  # finding qs01 and qs02
# finding MB about the middle of the spar
# results in finding the shear centre
shear_centre=(-1/Izz)*(2*t_skin*r**4*(1-m.pi/2)+(t_skin*l*r*lsk**2)/6+(r**2*t_spar*lsk*l)-(2*t_skin*r**2*lsk*l)+(t_skin*r*lsk**2*l)/2)+2*C[0]*r**2+2*C[1]*l*lsk
print(shear_centre)

