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
