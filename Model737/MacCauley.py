# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:01:36 2020

@author: daanv
"""

def MacCauley(x, value):
    i = 1 
    if x < value:
        i = 0 
    return i 
        
for i in range(10):
    print(MacCauley(i,6 ))
    