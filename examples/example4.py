"""
We caculate the Casimir pressure between two gold halfspaces in vacuum at 300K and separation of 1 micron.
We calculate the Casimir pressure between surfaces of varying plasma frequencies.
"""
import numpy as np
from califorcia import system
from califorcia.materials import vacuum
from scipy.constants import hbar
from scipy.constants import e as eV

class metal:
    def __init__(self):
        self.materialclass = "drude"
        self.wp = 9 * eV / hbar  # plasma frequency
        self.gamma = 0.035 * eV / hbar  # damping

    def epsilon(self, xi):
        return 1.+self.wp**2/xi/(xi+self.gamma)

T = 300     # in K
d = 1.e-6   # in m
m = metal()
s = system(T, d, m, m, vacuum)
print('pressure (gold):',s.pressure())

WP = np.linspace(1., 20., 9)
for wp in WP:
    m.wp = wp * eV / hbar
    s = system(T, d, m, m, vacuum)
    print('pressure (wp = ', wp, 'eV):', s.pressure())


