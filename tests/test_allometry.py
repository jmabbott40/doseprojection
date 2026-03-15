"""Tests for allometric scaling module."""

import pytest
from doseprojection.allometry import (
    simple_allometry,
    predict_human_cl,
    predict_human_vss,
    predict_human_thalf,
    rule_of_exponents,
)


def test_simple_allometry_clearance():
    # Rat CL = 25 mL/min/kg, BW = 0.15 kg
    # Human = 25 * (60/0.15)^0.75
    result = simple_allometry(25, 0.15, 60.0, 0.75)
    assert result > 25  # Should be larger (total)
    # Actually this scales total parameter, not per-kg
    # 25 * (400)^0.75 = 25 * 89.4 = 2236 (total human value)
    # This function returns the scaled value, not per-kg


def test_simple_allometry_identity():
    # Scaling to same weight should return same value
    result = simple_allometry(10.0, 5.0, 5.0, 0.75)
    assert result == pytest.approx(10.0)


def test_predict_human_cl_single_species():
    result = predict_human_cl({"rat": 25})
    assert "cl_human_mL_min_kg" in result
    assert result["method"] == "single_species_allometry"
    assert result["cl_human_mL_min_kg"] > 0


def test_predict_human_cl_multi_species():
    result = predict_human_cl({
        "mouse": 55,
        "rat": 25,
        "dog": 8,
    })
    assert "cl_human_mL_min_kg" in result
    assert "exponent" in result
    assert result["cl_human_mL_min_kg"] > 0


def test_predict_human_cl_with_custom_bw():
    result = predict_human_cl({"rat": (25, 0.25)})
    assert result["cl_human_mL_min_kg"] > 0


def test_predict_human_vss_single():
    result = predict_human_vss({"rat": 1.5})
    assert result["vss_human_L_kg"] > 0


def test_predict_human_thalf():
    # CL = 5 mL/min/kg = 0.3 L/h/kg, Vss = 1.0 L/kg
    # t1/2 = 0.693 * 1.0 / 0.3 = 2.31 h
    result = predict_human_thalf(5.0, 1.0)
    expected = 0.693 * 1.0 / (5.0 * 60 / 1000)
    assert result == pytest.approx(expected, rel=0.01)


def test_rule_of_exponents():
    assert rule_of_exponents(0.50) == "simple_allometry (may overestimate)"
    assert rule_of_exponents(0.65) == "simple_allometry"
    assert rule_of_exponents(0.85) == "MLP_correction_recommended"
    assert rule_of_exponents(1.1) == "brain_weight_correction_recommended"
