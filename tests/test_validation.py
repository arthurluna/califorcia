import pytest

from califorcia.compute import system
from califorcia.materials import gold, vacuum


def test_left_coating_length_mismatch_raises_value_error():
    with pytest.raises(ValueError, match="len\\(matL\\)=len\\(deltaL\\)\\+1"):
        system(300.0, 1e-6, [gold, gold], gold, vacuum, deltaL=[])


def test_right_coating_length_mismatch_raises_value_error():
    with pytest.raises(ValueError, match="len\\(matR\\)=len\\(deltaR\\)\\+1"):
        system(300.0, 1e-6, gold, [gold, gold], vacuum, deltaR=[])


def test_invalid_observable_raises_value_error():
    s = system(300.0, 1e-6, gold, gold, vacuum)
    with pytest.raises(ValueError, match="Supported values for 'observable'"):
        s.calculate("force")


def test_invalid_frequency_summation_method_raises_value_error():
    s = system(300.0, 1e-6, gold, gold, vacuum)
    with pytest.raises(ValueError, match="Supported values for fs"):
        s.calculate("energy", fs="invalid")


def test_non_positive_separation_raises_value_error():
    with pytest.raises(ValueError, match="must be positive"):
        system(300.0, 0.0, gold, gold, vacuum)


def test_negative_left_coating_thickness_raises_value_error():
    with pytest.raises(ValueError, match="plate L must be non-negative"):
        system(300.0, 1e-6, [gold, gold], gold, vacuum, deltaL=[-1e-9])


def test_negative_right_coating_thickness_raises_value_error():
    with pytest.raises(ValueError, match="plate R must be non-negative"):
        system(300.0, 1e-6, gold, [gold, gold], vacuum, deltaR=[-1e-9])


def test_more_than_four_materials_on_left_plate_raises_value_error():
    with pytest.raises(ValueError, match="At most 4 materials"):
        system(300.0, 1e-6, [gold, gold, gold, gold, gold], gold, vacuum, deltaL=[1e-9, 1e-9, 1e-9, 1e-9])


class UnsupportedMaterial:
    materialclass = "mystery"

    @staticmethod
    def epsilon(xi):
        return 1.0


class PlasmaMaterialWithoutWp:
    materialclass = "plasma"

    @staticmethod
    def epsilon(xi):
        return 2.0


def test_unsupported_materialclass_raises_value_error():
    with pytest.raises(ValueError, match="Unsupported materialclass"):
        system(300.0, 1e-6, UnsupportedMaterial(), gold, vacuum)


def test_plasma_material_without_wp_raises_value_error():
    with pytest.raises(ValueError, match="must define 'wp'"):
        system(300.0, 1e-6, PlasmaMaterialWithoutWp(), gold, vacuum)
