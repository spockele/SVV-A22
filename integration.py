# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:04:15 2020

@author: Mohammad
"""
# Analytic integration for a given polynomial
import numpy as np
class Polynomial:
    
    def __init__(self, *coefficients):
        self.coefficients = list(coefficients)
        self.degree = len(self.coefficients) - 1
        
    def order(self):
        self.degree = len(self.coefficients) - 1
    
    def integrated(self, integrals):
        for integral in range(integrals):
            integr_coeffs= []
            exponent = self.degree
            for i in range (len(self.coefficients)):
                integr_coeffs.append(self.coefficients[i]/(exponent+1))
                exponent -= 1
                
            integr_coeffs.append(0.0)
            self.coefficients = integr_coeffs
            self.order()
        #print('Polynomial of Degree', self.degree)
        #print('Polynomial given by', self.coefficients)

    
    def definite_integral(self, a, b, integrals):
        self.integrated(integrals)
        value = 0
        for i in range(self.degree):
            value += self.coefficients[i]*(b**(self.degree-i) - a**(self.degree-i))
        #print('Value of integral on (', a, ',', b, ') is', value)
        return value
    
p = Polynomial(1)


p_int = p.definite_integral(0, 1, 3)

zarray = np.arange(0,10,0.001)
datalist = zarray**2
si = np.array([0, 0])
for z in range(len(zarray)-1):
    fi0 = datalist[z]
    fi1 = datalist[z+1]
    zi0 = zarray[z]
    zi1 = zarray[z+1]
    spline = np.array([(fi1 - fi0)/(zi1 - zi0), fi0 - (fi1 - fi0)/(zi1 - zi0)*(zi0)]);
    # splines are formatted as a_n, a_n-1, ... , a_0
    si = np.vstack((si, spline))
si = np.delete(si, 0, 0)

integral = 0
for z in range(len(zarray)-1):
    sp = Polynomial(si[z,0], si[z,1])
    sp = sp.definite_integral(zarray[z], zarray[z+1], 1)
    integral += sp
print(integral)