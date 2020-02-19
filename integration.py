# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:04:15 2020

@author: Mohammad
"""
import numpy as np

# Analytic integration for a given polynomial
class Polynomial:
    def _init_(self, *coefficients):
        self.coefficients = list(coefficients)
        
    def _call_(self,x):
        layer = 0
        for coeff in self.coefficients:
            layer = coeff + x*layer
        return layer
    
    def definite_integral(self):
        integr_coeffs= []
        degree = len(self.coefficients)
        exponent = len(self.coefficients) + 1
        for i in range (len(self.coefficients) - 1):
            integr_coeffs.append(self.coefficients[i]/exponent)
            exponent -= 1
        return Polynomial(*integr_coeffs)
    
        