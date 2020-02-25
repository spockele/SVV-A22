# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:06:55 2020

@author: HJ Hoogendoorn


Calculating shear flows in aileron cross section
"""

from Moment_of_Inertia import *
import numpy as np
#import math as m

#To calculate the shear center a Sy with value 1 will be assumed
Sy = 1

#shear forces
Sy = 1


#Underneath the shear flows of the different elements of the aileron will be calculated
qb1 = -Sy/Izz