from scratch_1 import *
x=1.691   #m
T=1
G=28*10**9  #pa
theta=m.radians(26)   #radians
p=37900       #N
A1= 0.5*m.pi*r**2  # Area of left cell
A2=lsk*r # area of right cell
# making a 3 by 3 matrix
A1z=-0.884   #N
R1y=4253.957   #N
R2y=-5.2  #N
R3y=24.02 #N
moment=  0.278   #moment due to aerodynamic loading in Nm
mat_A=np.array([[2*A1,2*A2,0],[x1/(2*A1),y1/(2*A1),-1],[x2/(2*A2),y2/(2*A2),-1]])
mat_B=np.array([1,-(1*q_b_LC)/(2*A1),-(1*q_b_RC)/(2*A2)])
mat_C=np.linalg.solve(mat_A,mat_B)
J=T/(mat_C[2])
print(J)
#finding torque, clockwise positive
T1=r*(p*m.cos(theta)-A1z)-p*m.sin(theta)*(shear_centre+r)+(R1y+R2y+R3y)*(shear_centre)+(moment)
# integrating the torque in order to get the deflection along the span of the wing
twist=T1*x/(G*J)
print(twist)