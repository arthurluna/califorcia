from califorcia.compute import system
from califorcia.materials import gold, pec, vacuum
from scipy.constants import hbar, c, pi, k
from numpy.testing import assert_allclose

zeta3 = 1.2020569031595942

def cas_energy(d):
    # Casmir's result for energy at T=0
    return -hbar*c*pi**2/720./d**3

def cas_pressure(d):
    # Casmir's result for pressure at T=0
    return -hbar*c*pi**2/240./d**4

def cas_pressuregradient(d):
    # Casmir's result for pressure gradient at T=0
    return hbar*c*pi**2/60./d**5

def test_pec():
    # test Casimir interaction for pec-pec in vacuum against Casimir's exact results
    d = 1.e-6
    s = system(0., d, pec, pec, vacuum)
    assert_allclose(s.calculate('energy'), cas_energy(d))
    assert_allclose(s.calculate('pressure'), cas_pressure(d))
    assert_allclose(s.calculate('pressuregradient'), cas_pressuregradient(d))

def test_zero_frequency():
    d = 100.e-9
    s = system(1/k, d, pec, pec, vacuum)
    # analytical result from: Bordag et al., Advances in the Casimir effect (Oxford University press, 2009), Eq. (14.5)
    assert_allclose(s.calculate('energy', ht_limit=True), -zeta3/8/pi/d**2)
    assert_allclose(s.calculate('pressure', ht_limit=True), -zeta3/4/pi/d**3)
    assert_allclose(s.calculate('pressuregradient', ht_limit=True), 3*zeta3/4/pi/d**4)

    d = 100.e-9
    s = system(1 / k, d, gold, gold, vacuum)
    # analytical result from: Bordag et al., Advances in the Casimir effect (Oxford University press, 2009), Eq. (14.7)
    assert_allclose(s.calculate('energy', ht_limit=True), -zeta3 / 16 / pi / d ** 2)
    assert_allclose(s.calculate('pressure', ht_limit=True), -zeta3 / 8 / pi / d ** 3)
    assert_allclose(s.calculate('pressuregradient', ht_limit=True), 3 * zeta3 / 8 / pi / d ** 4)

if __name__ == '__main__':
    test_pec()
