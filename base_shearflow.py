"""
Calculations for the shear flows in all 6 elements of the aileron
"""
import math
from Moment_of_Inertia import tsk, ha, Izz, straightlen, Astr, zstr, ystr, Ca

h = ha / 2


def qb1(sy, theta1):
    """
    Calculates the base shear flow in the top half of the semicircle
    :param sy: The internal shear force
    :param theta1: The angle coordinate running upwards from 0 to pi/2
    :return: The base shear flow at theta1
    """
    if not 0 <= theta1 <= math.pi / 2:
        raise ValueError("theta is outside the range for qb1")
    qb = 0 - (sy / Izz) * tsk * h**2 * (math.cos(theta1) - 1)
    qb += Astr * (ystr[11] * (h * (1 - math.cos(theta1)) >= 0 - zstr[11]))
    return qb


def qb2(sy, s2):
    """
    Calculates the base shear flow in the top half of the spar
    :param sy: The internal shear force
    :param s2: The coordinate running upwards along the spar from 0 to h
    :return: The base shear flow at s2
    """
    if not 0 <= s2 <= h:
        raise ValueError("s2 is outside the range for qb2")
    return 0 - (sy / Izz) * tsk * (0.5 * s2**2)


def s3_str(z, y):
    return math.sqrt((z + h)**2 + y**2)


def qb3(sy, s3):
    """
    Calculates the base shear flow along the top skin
    :param sy: The internal shear force
    :param s3: The coordinate running towards the TE along the top skin from 0 to straightlen
    :return: The base shear flow at s3
    """
    if not 0 <= s3 <= straightlen:
        raise ValueError("s3 is outside the range for qb3")
    qb = 0 - (sy / Izz) * tsk * (h*s3 - h/(2*straightlen)*s3**2) + 0.5 * (sy / Izz) * tsk * h**2
    qb += Astr * ystr[11]
    qb += Astr * ystr[1] * (s3 >= s3_str(zstr[1], ystr[1]))
    qb += Astr * ystr[3] * (s3 >= s3_str(zstr[1], ystr[1]))
    qb += Astr * ystr[5] * (s3 >= s3_str(zstr[1], ystr[1]))
    qb += Astr * ystr[7] * (s3 >= s3_str(zstr[1], ystr[1]))
    qb += Astr * ystr[9] * (s3 >= s3_str(zstr[1], ystr[1]))
    return qb


def s4_str(z, y):
    return math.sqrt((z-Ca)**2 + y**2)


def qb4(sy, s4):
    """
    Calculates the base shear flow along the bottom skin
    :param sy: The internal shear force
    :param s4: The coordinate running towards the LE along the top skin from 0 to straightlen
    :return: The base shear flow at s4
    """
    if not 0 <= s4 <= straightlen:
        raise ValueError("s4 is outside the range for qb4")
    qb = 0 - (sy / Izz) * tsk * (h*straightlen/2 - h/(2*straightlen)*s4**2) + 0.5 * (sy / Izz) * tsk * h**2
    qb += Astr * (ystr[11] + ystr[1] + ystr[3] + ystr[5] + ystr[7] + ystr[9])
    qb += Astr * ystr[2] * (s4 >= s4_str(zstr[2], ystr[2]))
    qb += Astr * ystr[4] * (s4 >= s4_str(zstr[4], ystr[4]))
    qb += Astr * ystr[6] * (s4 >= s4_str(zstr[6], ystr[6]))
    qb += Astr * ystr[8] * (s4 >= s4_str(zstr[8], ystr[8]))
    qb += Astr * ystr[10] * (s4 >= s4_str(zstr[10], ystr[10]))
    return qb


def qb5(sy, s5):
    """
    Calculates the base shear flow in the bottom half of the spar
    :param sy: The internal shear force
    :param s5: The coordinate running downwards along the spar from 0 to h
    :return: The base shear flow at s5
    """
    if not 0 <= s5 <= h:
        raise ValueError("s5 is outside the range for qb5")
    return (sy / Izz) * tsk * (0.5 * s5**2)


def qb6(sy, theta6):
    """
    Calculates the base shear flow in the bottom half of the semicircle
    :param sy: The internal shear force
    :param theta6: The angle coordinate running downwards from 0 to pi/2
    :return: The base shear flow at theta6
    """
    if not 0 <= theta6 <= math.pi / 2:
        raise ValueError("theta is outside the range for qb6")
    qb = 0 - (sy / Izz) * tsk * h**2 * (math.cos(theta6) - 1)
    qb -= Astr * (ystr[12] * (h * (1 - math.sin(theta6)) <= 0 - zstr[12]))
    return qb
