import pytest
from numpy.testing import assert_allclose

from califorcia.compute import system
from califorcia.materials import ethanol, gold, gold_plasma, teflon, vacuum


class LorentzDielectric:
    def __init__(self):
        self.materialclass = "dielectric"

    def epsilon(self, xi):
        return 1.0 + 1.5e30 / (1.0e30 + xi**2)


class PlasmaMetal:
    def __init__(self, wp):
        self.materialclass = "plasma"
        self.wp = wp

    def epsilon(self, xi):
        return 1.0 + self.wp**2 / xi**2


def test_psd_and_msd_agree_for_dielectric_system():
    s = system(300.0, 200e-9, teflon, teflon, vacuum)
    pressure_psd = s.pressure(fs="psd", epsrel=1e-7)
    pressure_msd = s.pressure(fs="msd", epsrel=1e-7, N=200)
    assert_allclose(pressure_psd, pressure_msd, rtol=1e-5)


def test_psd_and_msd_agree_for_metallic_system():
    s = system(300.0, 1e-6, gold, gold, vacuum)
    energy_psd = s.energy(fs="psd", epsrel=1e-7)
    energy_msd = s.energy(fs="msd", epsrel=1e-7, N=250)
    assert_allclose(energy_psd, energy_msd, rtol=1e-5)


def test_swapping_left_and_right_planes_keeps_pressure_invariant():
    left = system(300.0, 120e-9, gold, teflon, ethanol)
    right = system(300.0, 120e-9, teflon, gold, ethanol)
    assert_allclose(left.pressure(), right.pressure(), rtol=1e-12)


def test_zero_thickness_coating_reduces_to_bare_halfspace():
    coated = system(300.0, 1e-6, [teflon, gold], gold, vacuum, deltaL=[0.0])
    bare = system(300.0, 1e-6, gold, gold, vacuum)
    assert_allclose(coated.pressure(), bare.pressure(), rtol=1e-10)


def test_identical_coating_and_substrate_match_uncoated_case():
    coated = system(300.0, 1e-6, [gold, gold], gold, vacuum, deltaL=[50e-9])
    bare = system(300.0, 1e-6, gold, gold, vacuum)
    assert_allclose(coated.energy(), bare.energy(), rtol=1e-10)


def test_left_and_right_coating_are_symmetric():
    left = system(300.0, 1e-6, [teflon, gold], gold, vacuum, deltaL=[50e-9])
    right = system(300.0, 1e-6, gold, [teflon, gold], vacuum, deltaR=[50e-9])
    assert_allclose(left.pressuregradient(), right.pressuregradient(), rtol=1e-10)


def test_pressure_matches_finite_difference_of_energy():
    d = 500e-9
    step = 1e-9
    s = system(300.0, d, gold, gold, vacuum)
    s_plus = system(300.0, d + step, gold, gold, vacuum)
    s_minus = system(300.0, d - step, gold, gold, vacuum)
    finite_difference = -(s_plus.energy() - s_minus.energy()) / (2 * step)
    assert_allclose(s.pressure(), finite_difference, rtol=5e-4)


def test_custom_dielectric_material_can_be_used_in_system():
    s = system(300.0, 100e-9, gold, LorentzDielectric(), ethanol)
    assert pytest.approx(0.22099374769016233, rel=1e-12) == s.pressure()


def test_custom_plasma_material_matches_builtin_gold_plasma():
    s_builtin = system(300.0, 1e-6, gold_plasma, gold_plasma, vacuum)
    s_custom = system(300.0, 1e-6, PlasmaMetal(gold_plasma.wp), PlasmaMetal(gold_plasma.wp), vacuum)
    assert_allclose(s_custom.pressure(), s_builtin.pressure(), rtol=1e-12)
