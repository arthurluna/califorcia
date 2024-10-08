'''
We calculate the Casimir pressure between two gold halfspaces in vacuum at 300K, where one surfaces is coated with a
50nm of Teflon film. The plates are separated by a distance of one micron.
'''

from califorcia import system
from califorcia.materials import gold, teflon, vacuum

T = 300     # in K
d = 1.e-6   # in m
s = system(300, d, [teflon, gold], gold, vacuum, deltaL=[50.e-9])
print(s.pressure())
