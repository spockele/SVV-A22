"""
Calculations for the shear flows in all 6 elements of the aileron
"""
import math

def qb1(sy, theta):
    return 0 - (sy / I_zz) * tsk * h**2 * (math.cos(theta) - 1)