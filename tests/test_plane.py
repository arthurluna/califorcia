from math import sqrt, exp

from numpy.testing import assert_allclose

from califorcia.plane import def_fresnel_coefficients, def_reflection_coeff, kappa
from califorcia.materials import gold_drude, gold_plasma, pec, teflon, vacuum


def legacy_reflection_coeff(medium, materials, thicknesses):
    nlayers = len(materials)

    if nlayers == 1:
        return def_fresnel_coefficients(medium, materials[0])

    if nlayers == 2:
        mat1 = materials[0]
        mat2 = materials[1]
        d1 = thicknesses[0]
        fresnelm1 = def_fresnel_coefficients(medium, mat1)
        fresnel12 = def_fresnel_coefficients(mat1, mat2)

        def reflection_coeff(k0, k):
            rTMm1, rTEm1 = fresnelm1(k0, k)
            rTM12, rTE12 = fresnel12(k0, k)
            kappa1 = kappa(mat1, k0, k)
            rTM = (rTMm1 + rTM12 * exp(-2 * kappa1 * d1)) / (1 + rTMm1 * rTM12 * exp(-2 * kappa1 * d1))
            rTE = (rTEm1 + rTE12 * exp(-2 * kappa1 * d1)) / (1 + rTEm1 * rTE12 * exp(-2 * kappa1 * d1))
            return rTM, rTE

        return reflection_coeff

    if nlayers == 3:
        mat1 = materials[0]
        mat2 = materials[1]
        mat3 = materials[2]
        d1 = thicknesses[0]
        d2 = thicknesses[1]
        fresnelm1 = def_fresnel_coefficients(medium, mat1)
        fresnel12 = def_fresnel_coefficients(mat1, mat2)
        fresnel23 = def_fresnel_coefficients(mat2, mat3)

        def reflection_coeff(k0, k):
            rTMm1, rTEm1 = fresnelm1(k0, k)
            rTM12, rTE12 = fresnel12(k0, k)
            rTM23, rTE23 = fresnel23(k0, k)
            kappa1 = kappa(mat1, k0, k)
            kappa2 = kappa(mat2, k0, k)

            rTM12_eff = (rTM12 + rTM23 * exp(-2 * kappa2 * d2)) / (1 + rTM12 * rTM23 * exp(-2 * kappa2 * d2))
            rTE12_eff = (rTE12 + rTE23 * exp(-2 * kappa2 * d2)) / (1 + rTE12 * rTE23 * exp(-2 * kappa2 * d2))

            rTM = (rTMm1 + rTM12_eff * exp(-2 * kappa1 * d1)) / (1 + rTMm1 * rTM12_eff * exp(-2 * kappa1 * d1))
            rTE = (rTEm1 + rTE12_eff * exp(-2 * kappa1 * d1)) / (1 + rTEm1 * rTE12_eff * exp(-2 * kappa1 * d1))
            return rTM, rTE

        return reflection_coeff

    if nlayers == 4:
        mat1 = materials[0]
        mat2 = materials[1]
        mat3 = materials[2]
        mat4 = materials[3]
        d1 = thicknesses[0]
        d2 = thicknesses[1]
        d3 = thicknesses[2]
        fresnelm1 = def_fresnel_coefficients(medium, mat1)
        fresnel12 = def_fresnel_coefficients(mat1, mat2)
        fresnel23 = def_fresnel_coefficients(mat2, mat3)
        fresnel34 = def_fresnel_coefficients(mat3, mat4)

        def reflection_coeff(k0, k):
            rTMm1, rTEm1 = fresnelm1(k0, k)
            rTM12, rTE12 = fresnel12(k0, k)
            rTM23, rTE23 = fresnel23(k0, k)
            rTM34, rTE34 = fresnel34(k0, k)
            kappa1 = kappa(mat1, k0, k)
            kappa2 = kappa(mat2, k0, k)
            kappa3 = kappa(mat3, k0, k)

            rTM23_eff = (rTM23 + rTM34 * exp(-2 * kappa3 * d3)) / (1 + rTM23 * rTM34 * exp(-2 * kappa3 * d3))
            rTE23_eff = (rTE23 + rTE34 * exp(-2 * kappa3 * d3)) / (1 + rTE23 * rTE34 * exp(-2 * kappa3 * d3))

            rTM12_eff = (rTM12 + rTM23_eff * exp(-2 * kappa2 * d2)) / (1 + rTM12 * rTM23_eff * exp(-2 * kappa2 * d2))
            rTE12_eff = (rTE12 + rTE23_eff * exp(-2 * kappa2 * d2)) / (1 + rTE12 * rTE23_eff * exp(-2 * kappa2 * d2))

            rTM = (rTMm1 + rTM12_eff * exp(-2 * kappa1 * d1)) / (1 + rTMm1 * rTM12_eff * exp(-2 * kappa1 * d1))
            rTE = (rTEm1 + rTE12_eff * exp(-2 * kappa1 * d1)) / (1 + rTEm1 * rTE12_eff * exp(-2 * kappa1 * d1))
            return rTM, rTE

        return reflection_coeff

    raise ValueError("legacy_reflection_coeff only supports up to 4 materials.")


def test_kappa_zero_frequency_depends_on_materialclass():
    kval = 2.0
    assert_allclose(kappa(teflon, 0.0, kval), kval)
    assert_allclose(kappa(gold_drude, 0.0, kval), kval)
    assert_allclose(kappa(gold_plasma, 0.0, kval), sqrt((gold_plasma.wp / 299792458.0) ** 2 + kval**2))


def test_fresnel_coefficients_for_pec_are_exact():
    fresnel = def_fresnel_coefficients(vacuum, pec)
    assert fresnel(0.0, 1.0) == (1.0, -1.0)
    assert fresnel(1.0, 1.0) == (1.0, -1.0)


def test_dielectric_dielectric_zero_frequency_te_mode_vanishes():
    fresnel = def_fresnel_coefficients(vacuum, teflon)
    r_tm, r_te = fresnel(0.0, 1.0)
    assert_allclose(r_te, 0.0, atol=0.0)
    assert r_tm > 0.0


def test_zero_thickness_multilayer_reflection_reduces_to_single_interface():
    multilayer = def_reflection_coeff(vacuum, [teflon, gold_drude], [0.0])
    single = def_reflection_coeff(vacuum, [gold_drude], [])
    assert_allclose(multilayer(1.2, 0.7), single(1.2, 0.7), rtol=1e-12)


def test_identical_material_layers_reduce_to_same_reflection():
    multilayer = def_reflection_coeff(vacuum, [gold_drude, gold_drude, gold_drude, gold_drude], [5e-9, 7e-9, 9e-9])
    single = def_reflection_coeff(vacuum, [gold_drude], [])
    assert_allclose(multilayer(1.2, 0.7), single(1.2, 0.7), rtol=1e-12)


def test_many_identical_material_layers_reduce_to_same_reflection():
    multilayer = def_reflection_coeff(
        vacuum,
        [gold_drude, gold_drude, gold_drude, gold_drude, gold_drude, gold_drude],
        [5e-9, 7e-9, 9e-9, 11e-9, 13e-9],
    )
    single = def_reflection_coeff(vacuum, [gold_drude], [])
    assert_allclose(multilayer(1.2, 0.7), single(1.2, 0.7), rtol=1e-12)


def test_pec_top_layer_ignores_underlying_stack():
    coated = def_reflection_coeff(vacuum, [pec, gold_drude, teflon, gold_drude], [5e-9, 7e-9, 9e-9])
    bare = def_reflection_coeff(vacuum, [pec], [])
    assert_allclose(coated(1.2, 0.7), bare(1.2, 0.7), rtol=1e-12)


def test_recursive_reflection_matches_legacy_formula_for_one_material():
    recursive = def_reflection_coeff(vacuum, [teflon], [])
    legacy = legacy_reflection_coeff(vacuum, [teflon], [])
    assert_allclose(recursive(1.2, 0.7), legacy(1.2, 0.7), rtol=1e-12)


def test_recursive_reflection_matches_legacy_formula_for_two_materials():
    recursive = def_reflection_coeff(vacuum, [teflon, gold_drude], [5e-9])
    legacy = legacy_reflection_coeff(vacuum, [teflon, gold_drude], [5e-9])
    assert_allclose(recursive(1.2, 0.7), legacy(1.2, 0.7), rtol=1e-12)


def test_recursive_reflection_matches_legacy_formula_for_three_materials():
    recursive = def_reflection_coeff(vacuum, [teflon, gold_drude, teflon], [5e-9, 7e-9])
    legacy = legacy_reflection_coeff(vacuum, [teflon, gold_drude, teflon], [5e-9, 7e-9])
    assert_allclose(recursive(1.2, 0.7), legacy(1.2, 0.7), rtol=1e-12)


def test_recursive_reflection_matches_legacy_formula_for_four_materials():
    recursive = def_reflection_coeff(vacuum, [teflon, gold_drude, teflon, gold_drude], [5e-9, 7e-9, 9e-9])
    legacy = legacy_reflection_coeff(vacuum, [teflon, gold_drude, teflon, gold_drude], [5e-9, 7e-9, 9e-9])
    assert_allclose(recursive(1.2, 0.7), legacy(1.2, 0.7), rtol=1e-12)
