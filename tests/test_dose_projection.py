"""Tests for dose projection module."""

import pytest
from doseprojection.dose_projection import (
    efficacious_dose_mg_kg,
    dose_from_target_css,
    unbound_concentration,
    steady_state_css,
    project_animal_dose,
    _resolve_f,
)


# --- Route resolution tests ---

def test_resolve_f_iv():
    """IV route should always return F=1.0."""
    assert _resolve_f(None, "iv") == 1.0
    assert _resolve_f(0.5, "iv") == 1.0  # input ignored


def test_resolve_f_po():
    assert _resolve_f(0.5, "po") == 0.5


def test_resolve_f_sc():
    assert _resolve_f(0.7, "sc") == 0.7


def test_resolve_f_po_requires_f():
    with pytest.raises(ValueError, match="required"):
        _resolve_f(None, "po")


def test_resolve_f_sc_requires_f():
    with pytest.raises(ValueError, match="required"):
        _resolve_f(None, "sc")


def test_resolve_f_invalid_route():
    with pytest.raises(ValueError, match="Route"):
        _resolve_f(0.5, "topical")


def test_resolve_f_invalid_value():
    with pytest.raises(ValueError, match="between 0"):
        _resolve_f(0, "po")


# --- Core function tests (PO route) ---

def test_unbound_concentration():
    assert unbound_concentration(100, 0.1) == pytest.approx(10.0)
    assert unbound_concentration(100, 1.0) == pytest.approx(100.0)
    assert unbound_concentration(100, 0.0) == pytest.approx(0.0)


def test_steady_state_css():
    # 10 mg/kg, F=0.5, CL=25 mL/min/kg, tau=24h
    css = steady_state_css(10, 25, 24, f=0.5, route="po")
    # css = 10 * 0.5 * 1e6 / (25*60 * 24) = 5e6 / 36000 = 138.9 ng/mL
    assert css == pytest.approx(138.9, rel=0.01)


def test_efficacious_dose_mg_kg():
    # IC50 = 100 nM, MW = 400, CL = 25, fu = 0.1, F = 0.5, tau = 24h
    dose = efficacious_dose_mg_kg(
        ic50_nm=100, mw=400, cl_mL_min_kg=25, fu=0.1, tau_h=24, f=0.5
    )
    # IC50 in mg/L = 100 * 400 / 1e6 = 0.04 mg/L
    # CL in mL/h/kg = 25 * 60 = 1500
    # dose = 0.04 * 1500 * 24 / (0.1 * 0.5 * 1000) = 1440 / 50 = 28.8 mg/kg
    assert dose == pytest.approx(28.8, rel=0.01)


def test_efficacious_dose_with_coverage():
    dose_1x = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, f=0.5, coverage_multiple=1.0)
    dose_3x = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, f=0.5, coverage_multiple=3.0)
    assert dose_3x == pytest.approx(dose_1x * 3, rel=0.01)


def test_dose_from_target_css():
    # Target 500 ng/mL, CL=25, F=0.5, tau=24
    dose = dose_from_target_css(500, 25, 24, f=0.5)
    # 500 * 1500 * 24 / (0.5 * 1e6) = 18e6 / 5e5 = 36 mg/kg
    assert dose == pytest.approx(36.0, rel=0.01)


def test_project_animal_dose():
    result = project_animal_dose(
        ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
        fu_animal=0.1, tau_h=24, f_animal=0.5, route="po"
    )
    assert "dose_mg_kg" in result
    assert "target_cu_ng_mL" in result
    assert "predicted_css_total_ng_mL" in result
    assert result["route"] == "po"
    assert result["dose_mg_kg"] > 0


# --- IV route tests ---

def test_efficacious_dose_iv():
    """IV dose should be lower than PO dose (F=1 vs F<1)."""
    dose_po = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, f=0.5, route="po")
    dose_iv = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, route="iv")
    assert dose_iv == pytest.approx(dose_po * 0.5, rel=0.01)


def test_efficacious_dose_iv_f_ignored():
    """F parameter should be ignored for IV route."""
    dose1 = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, f=0.3, route="iv")
    dose2 = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, f=0.9, route="iv")
    assert dose1 == pytest.approx(dose2)


def test_steady_state_css_iv():
    """IV Css should be higher than PO for same dose (F=1 vs F<1)."""
    css_iv = steady_state_css(10, 25, 24, route="iv")
    css_po = steady_state_css(10, 25, 24, f=0.5, route="po")
    assert css_iv == pytest.approx(css_po * 2, rel=0.01)


def test_dose_from_target_css_iv():
    """IV dose should be lower than PO dose for same target Css."""
    dose_iv = dose_from_target_css(500, 25, 24, route="iv")
    dose_po = dose_from_target_css(500, 25, 24, f=0.5, route="po")
    assert dose_iv == pytest.approx(dose_po * 0.5, rel=0.01)


def test_project_animal_dose_iv():
    result = project_animal_dose(
        ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
        fu_animal=0.1, tau_h=24, route="iv"
    )
    assert result["route"] == "iv"
    assert result["dose_mg_kg"] > 0
    result_po = project_animal_dose(
        ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
        fu_animal=0.1, tau_h=24, f_animal=0.5, route="po"
    )
    assert result["dose_mg_kg"] == pytest.approx(
        result_po["dose_mg_kg"] * 0.5, rel=0.01
    )


# --- SC route tests ---

def test_efficacious_dose_sc():
    """SC dose with F=0.7 should be between IV and PO (F=0.5)."""
    dose_iv = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, route="iv")
    dose_sc = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, f=0.7, route="sc")
    dose_po = efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, f=0.5, route="po")
    assert dose_iv < dose_sc < dose_po


def test_project_animal_dose_sc():
    result = project_animal_dose(
        ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
        fu_animal=0.1, tau_h=24, f_animal=0.8, route="sc"
    )
    assert result["route"] == "sc"
    assert result["dose_mg_kg"] > 0


def test_sc_requires_f():
    with pytest.raises(ValueError, match="required"):
        efficacious_dose_mg_kg(100, 400, 25, 0.1, 24, route="sc")
