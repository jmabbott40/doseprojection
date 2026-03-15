"""Tests for IVIVE module."""

import math
import pytest
from doseprojection.ivive import (
    calc_clint_invitro,
    scale_clint,
    predict_hepatic_clearance,
    predict_extraction_ratio,
    predict_fh,
    ivive_workflow,
)


def test_calc_clint_invitro():
    # t1/2 = 30 min, volume = 1000 uL, protein = 0.5 mg
    # CLint = 0.693 / 30 * (1000 / 0.5) = 0.0231 * 2000 = 46.2
    result = calc_clint_invitro(30, 1000, 0.5)
    assert abs(result - 46.2) < 0.1


def test_calc_clint_invitro_negative_thalf():
    with pytest.raises(ValueError):
        calc_clint_invitro(-5)


def test_scale_clint_human():
    # CLint = 46.2 uL/min/mg, human: MPPGL=32, liver=21.4
    # scaled = 46.2 * 32 * 21.4 / 1000 = 31.64
    result = scale_clint(46.2, "human")
    expected = 46.2 * 32 * 21.4 / 1000
    assert abs(result - expected) < 0.01


def test_scale_clint_rat():
    result = scale_clint(46.2, "rat")
    expected = 46.2 * 45 * 40 / 1000
    assert abs(result - expected) < 0.01


def test_scale_clint_unknown_species():
    with pytest.raises(ValueError):
        scale_clint(46.2, "fish")


def test_predict_hepatic_clearance_human():
    # CLint_scaled = 31.64, fu = 0.1, human QH = 20.7
    clint = 31.64
    fu = 0.1
    qh = 20.7
    expected = (qh * fu * clint) / (qh + fu * clint)
    result = predict_hepatic_clearance(clint, fu, "human")
    assert abs(result - expected) < 0.01


def test_predict_extraction_ratio():
    clint = 31.64
    fu = 0.1
    qh = 20.7
    expected = (fu * clint) / (qh + fu * clint)
    result = predict_extraction_ratio(clint, fu, "human")
    assert abs(result - expected) < 0.001


def test_predict_fh():
    assert predict_fh(0.3) == pytest.approx(0.7)
    assert predict_fh(0.0) == pytest.approx(1.0)
    assert predict_fh(1.0) == pytest.approx(0.0)


def test_ivive_workflow():
    result = ivive_workflow(t_half_min=30, fu=0.1, species="human")
    assert "clint_invitro_uL_min_mg" in result
    assert "cl_hepatic_mL_min_kg" in result
    assert "fh" in result
    assert 0 < result["cl_hepatic_mL_min_kg"] < 20.7  # Can't exceed QH
    assert 0 < result["fh"] <= 1.0
