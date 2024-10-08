"""
We calculate the Casimir pressure between a gold halfspace and a halfspace made of a user-specified material.
The whole system is immersed in ethanol at 300K and the two plates are separated by 100nm.
"""
from califorcia import system
from califorcia.materials import gold, ethanol

class user_material:
    def __init__(self):
        self.materialclass = 'dielectric'
    def epsilon(self, xi):
        wj = 1.911e15   # = 1.26 eV/hbar
        cj = 1.282      # dimensionless
        return 1. + cj * wj ** 2 / (wj ** 2 + xi ** 2)
u = user_material()

T = 300         # in K
d = 100.e-9     # in m
s = system(T, d, gold, u, ethanol)
print(s.pressure())


