"""
Fused silica according to wikipedia (Sellmeier equation)

"""
import numpy as np
from scipy.constants import pi, c

materialclass = "dielectric"

B1 = 0.696166300
B2 = 0.407942600
B3 = 0.897479400
C1 = 4.67914826e-3*1.e-12
C2 = 1.35120631e-2*1.e-12
C3 = 97.9325*1.e-12

def epsilon(xi):
    l = 2*pi*c/xi
    return 1 + B1/(1+C1/l**2) + B2/(1+C2/l**2) + B3/(1+C3/l**2)
