import numpy as np
import matplotlib.pyplot as plt
from integration import Polynomial

class Interpolate:

    def d3spline(self, data, zarray, x):
        # contructing the linalg equation Az=d to solve for z
        d = np.array([]);
        h = np.array([]);
        # finding the step sizes
        for z in range(len(zarray)-1):
            h = np.append(h, zarray[z+1]-zarray[z]);
        # finding A
        A = np.zeros((len(zarray),len(zarray)));
        A[0,0] = 1
        A[-1,-1] = 1
        A[0,1] = 0
        A[-1,-2] = 0
        # finding d
        for z in range(len(zarray)):
            if z == 0:
                d = np.append(d, 0);
            else:
                if z == (len(zarray)-1):
                    d = np.append(d, 0);
                else:
                    d = np.append(d, 6*(((data[z+1,x]-data[z,x])/(zarray[z+1]-zarray[z]))-((data[z,x]-data[z-1,x])/(zarray[z]-zarray[z-1])))/(zarray[z+1]-zarray[z-1]));
                    A[z,z-1] = h[z-1]/(h[z]+h[z-1])
                    A[z,z+1] = h[z]/(h[z]+h[z-1])
                    A[z,z] = 2

        # solving the equation
        M = np.linalg.solve(A, np.transpose(d))
        print("Matrix A \n",A)
        print("Matrix d \n",d)
        print("Matrix M, from AM=d \n",M)
        splines = np.array([0,0,0,0]) #splines baby
        for i in range(len(zarray)-1):
            si = np.array([])
            Ma = M[i+1]/(6*h[i])
            Mb = M[i]/(6*h[i])
            Mc = data[i,x]/h[i] - M[i]*h[i]/6
            Md = data[i+1,x]/h[i] - M[i+1]*h[i]/6
            si = np.append(si, Ma-Mb)
            si = np.append(si, -3*zarray[i]*Ma + 3*zarray[i+1]*Mb)
            si = np.append(si, 3*zarray[i]**2*Ma - 3*zarray[i+1]**2*Mb - Mc + Md)
            si = np.append(si, -zarray[i]**3*Ma + zarray[i+1]**3*Mb + zarray[i+1]*Mc - zarray[i]*Md)
            splines = np.vstack((splines, si))
        splines = np.delete(splines, 0,0)
        print("Spline functions: \n",splines)
        print("Number of splines: \n",len(splines))
        integral = 0
        for z in range(len(zarray)-1):
            sp = Polynomial(splines[z,0], splines[z,1], splines[z,2], splines[z,3])
            sp = sp.definite_integral(zarray[z], zarray[z+1], 2)
            integral += sp
        print("Integral: ",integral)
            
