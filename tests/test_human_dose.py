"""Tests for human dose projection module."""

import pytest
from doseprojection.human_dose import (
    calc_hed,
    calc_mrsd,
    hed_from_noael,
    bsa_conversion,
)


def test_calc_hed_rat():
    # Rat NOAEL = 50 mg/kg, Km_rat = 6, Km_human = 37
    # HED = 50 * (6/37) = 8.108 mg/kg
    result = calc_hed(50, "rat")
    assert result["hed_mg_kg"] == pytest.approx(50 * 6 / 37, rel=0.001)
    assert result["hed_total_mg"] == pytest.approx(50 * 6 / 37 * 60, rel=0.001)


def test_calc_hed_mouse():
    # Mouse dose 100 mg/kg, Km_mouse = 3
    # HED = 100 * (3/37) = 8.108 mg/kg
    result = calc_hed(100, "mouse")
    assert result["hed_mg_kg"] == pytest.approx(100 * 3 / 37, rel=0.001)


def test_calc_hed_dog():
    # Dog NOAEL 10 mg/kg, Km_dog = 20
    # HED = 10 * (20/37) = 5.405 mg/kg
    result = calc_hed(10, "dog")
    assert result["hed_mg_kg"] == pytest.approx(10 * 20 / 37, rel=0.001)


def test_calc_hed_unknown_species():
    with pytest.raises(ValueError):
        calc_hed(50, "fish")


def test_calc_mrsd():
    result = calc_mrsd(8.108, safety_factor=10)
    assert result["mrsd_mg_kg"] == pytest.approx(0.8108, rel=0.001)
    assert result["mrsd_total_mg"] == pytest.approx(0.8108 * 60, rel=0.001)


def test_hed_from_noael_rat():
    result = hed_from_noael(50, "rat", safety_factor=10)
    assert result["hed_mg_kg"] == pytest.approx(50 * 6 / 37, rel=0.001)
    assert result["mrsd_mg_kg"] == pytest.approx(50 * 6 / 37 / 10, rel=0.001)
    assert result["noael_mg_kg"] == 50
    assert result["species"] == "rat"


def test_bsa_conversion_rat_to_human():
    result = bsa_conversion(50, "rat", "human")
    assert result["dose_mg_kg"] == pytest.approx(50 * 6 / 37, rel=0.001)


def test_bsa_conversion_mouse_to_rat():
    result = bsa_conversion(100, "mouse", "rat")
    # 100 * (3/6) = 50
    assert result["dose_mg_kg"] == pytest.approx(50.0, rel=0.001)


def test_bsa_conversion_identity():
    result = bsa_conversion(50, "rat", "rat")
    assert result["dose_mg_kg"] == pytest.approx(50.0)
