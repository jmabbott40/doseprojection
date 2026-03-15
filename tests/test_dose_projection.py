"""Tests for dose projection module."""

import pytest
from doseprojection.dose_projection import (
    efficacious_dose_mg_kg,
    dose_from_target_css,
    unbound_concentration,
    steady_state_css,
    project_animal_dose,
)


def test_unbound_concentration():
    assert unbound_concentration(100, 0.1) == pytest.approx(10.0)
    assert unbound_concentration(100, 1.0) == pytest.approx(100.0)
    assert unbound_concentration(100, 0.0) == pytest.approx(0.0)


def test_steady_state_css():
    # 10 mg/kg, F=0.5, CL=25 mL/min/kg, tau=24h
    css = steady_state_css(10, 0.5, 25, 24)
    # css = 10 * 0.5 * 1e6 / (25*60 * 24) = 5e6 / 36000 = 138.9 ng/mL
    assert css == pytest.approx(138.9, rel=0.01)


def test_efficacious_dose_mg_kg():
    # IC50 = 100 nM, MW = 400, CL = 25, fu = 0.1, F = 0.5, tau = 24h
    dose = efficacious_dose_mg_kg(
        ic50_nm=100, mw=400, cl_mL_min_kg=25, fu=0.1, f=0.5, tau_h=24
    )
    # IC50 in mg/L = 100 * 400 / 1e6 = 0.04 mg/L
    # CL in mL/h/kg = 25 * 60 = 1500
    # dose = 0.04 * 1500 * 24 / (0.1 * 0.5 * 1000) = 1440 / 50 = 28.8 mg/kg
    assert dose == pytest.approx(28.8, rel=0.01)


def test_efficacious_dose_with_coverage():
    dose_1x = efficacious_dose_mg_kg(100, 400, 25, 0.1, 0.5, 24, 1.0)
    dose_3x = efficacious_dose_mg_kg(100, 400, 25, 0.1, 0.5, 24, 3.0)
    assert dose_3x == pytest.approx(dose_1x * 3, rel=0.01)


def test_dose_from_target_css():
    # Target 500 ng/mL, CL=25, F=0.5, tau=24
    dose = dose_from_target_css(500, 25, 0.5, 24)
    # 500 * 1500 * 24 / (0.5 * 1e6) = 18e6 / 5e5 = 36 mg/kg
    assert dose == pytest.approx(36.0, rel=0.01)


def test_project_animal_dose():
    result = project_animal_dose(
        ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
        fu_animal=0.1, f_animal=0.5, tau_h=24
    )
    assert "dose_mg_kg" in result
    assert "target_cu_ng_mL" in result
    assert "predicted_css_total_ng_mL" in result
    assert result["dose_mg_kg"] > 0
