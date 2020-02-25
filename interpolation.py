import numpy as np
from sympy import Symbol, simplify, Poly
from functools import reduce
import operator
from integration import Polynomial
import matplotlib.pyplot as plt

def d3splines(data, zarray, x):
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
#    print("Matrix A \n",A)
#    print("Matrix d \n",d)
#    print("Matrix M, from AM=d \n",M)
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
#    print("Spline functions: \n",splines)
#    print("Number of splines: \n",len(splines))
   
    return splines

def intSpline(splines, xarray):
    newPoints = np.array([0])
    integral = 0
    for x in range(len(xarray)-1):
        if len(splines[0]) ==  4:
            p = Polynomial(splines[x,0], splines[x,1], splines[x,2], splines[x,3])
        else:
            p = Polynomial(splines[x,0], splines[x,1])            
        p = p.definite_integral(xarray[x], xarray[x+1], 1)
        integral += p
        newPoints = np.append(newPoints, integral)
    newPoints = np.append(newPoints,0)
    si = np.array([0, 0])
    for z in range(len(xarray)-1):
        fi0 = newPoints[z]
        fi1 = newPoints[z+1]
        zi0 = xarray[z]
        zi1 = xarray[z+1]
        spline = np.array([(fi1 - fi0)/(zi1 - zi0), fi0 - (fi1 - fi0)/(zi1 - zi0)*(zi0)]);
        # splines are formatted as a_n, a_n-1, ... , a_0
        si = np.vstack((si, spline))
    newSplines = np.delete(si, 0, 0)
    return newSplines, integral
    

def resultantload(splines, zarray, z_hinge):
    #all the splines for the cross section are in this array now we integrate
    integral = 0
    weighted = 0
    for z in range(len(zarray)-1):
        p = Polynomial(splines[z,0], splines[z,1], splines[z,2], splines[z,3])
        p = p.definite_integral(zarray[z], zarray[z+1], 1)
        integral += p
        pz = Polynomial(splines[z,0], splines[z,1], splines[z,2], splines[z,3], 0)
        pz = pz.definite_integral(zarray[z], zarray[z+1], 1)
        weighted += pz
    dq = integral
    if integral != 0:
        cop = abs(weighted/integral)
    else:
        cop = 0
#    print(cop)
#    print(dq)
    dqm = abs(dq*(z_hinge-cop))
    
    return dq, dqm, cop

X = [ 0.00000000e+00,  1.24043255e-03,  6.19488340e-03,  1.60747104e-02,
        3.08219348e-02,  5.03500141e-02,  7.45443496e-02,  1.03262959e-01,
        1.36337311e-01,  1.73573312e-01,  2.14752446e-01,  2.59633057e-01,
        3.07951768e-01,  3.59425026e-01,  4.13750765e-01,  4.70610180e-01,
        5.29669598e-01,  5.90582435e-01,  6.52991229e-01,  7.16529741e-01,
        7.80825101e-01,  8.45500000e-01,  9.10174899e-01,  9.74470259e-01,
        1.03800877e+00,  1.10041757e+00,  1.16133040e+00,  1.22038982e+00,
        1.27724924e+00,  1.33157497e+00,  1.38304823e+00,  1.43136694e+00,
        1.47624755e+00,  1.51742669e+00,  1.55466269e+00,  1.58773704e+00,
        1.61645565e+00,  1.64064999e+00,  1.66017807e+00,  1.67492529e+00,
        1.68480512e+00,  1.68975957e+00, 1.69100000e+00]
Y = [ 0.        , -0.20453439, -0.45610082, -0.73152864, -1.00626597,
       -1.2745728 , -1.53301278, -1.7787306 , -2.00911488, -2.22178707,
       -2.41469215, -2.58625995, -2.73543299, -2.86190878, -2.96620537,
       -3.04955116, -3.11399564, -3.16205113, -3.19657674, -3.22017462,
       -3.23515695, -3.24279353, -3.24331433, -3.23565556, -3.21761604,
       -3.18604238, -3.13733611, -3.06766209, -2.97359769, -2.85247519,
       -2.70267228, -2.5237452 , -2.31673031, -2.08407758, -1.8296829 ,
       -1.55883877, -1.27815243, -0.99552932, -0.72027252, -0.46354811,
       -0.24019092, -0.07830771,  0.        ]

