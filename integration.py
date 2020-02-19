# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:04:15 2020

@author: Mohammad
"""
import numpy as np

# Analytic integration for a given polynomial
class Polynomial:
    
    def __init__(self, *coefficients):
        self.coefficients = list(coefficients)
       
    def __call__(self,x):
        layer = 0
        for coeff in self.coefficients:
            layer = coeff + x*layer
        return layer
    
    def definite_integral(self, a, b):
        integr_coeffs= []
        exponent = len(self.coefficients) - 1
        for i in range (len(self.coefficients)):
            integr_coeffs.append(self.coefficients[i]/(exponent+1))
            exponent -= 1
        
        value = 0
        degree = len(integr_coeffs)
        for i in range(degree):
            value += integr_coeffs[i]*(b^(degree-i) - a^(degree-i))
        
        return value
    
p = Polynomial(1, 0, 0, 0, 0)

p_int = p.definite_integral(0, 1)

print(p_int)

