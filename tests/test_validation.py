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
