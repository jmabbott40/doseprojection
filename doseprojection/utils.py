"""
Unit conversion helpers and utility functions for dose projection.
"""

import math


def nm_to_um(conc_nm):
    """Convert nanomolar to micromolar."""
    return conc_nm / 1000.0


def um_to_nm(conc_um):
    """Convert micromolar to nanomolar."""
    return conc_um * 1000.0


def nm_to_mg_per_ml(conc_nm, mw):
    """Convert nanomolar concentration to mg/mL given molecular weight (g/mol)."""
    return conc_nm * mw / 1e9


def mg_per_ml_to_nm(conc_mg_ml, mw):
    """Convert mg/mL to nanomolar given molecular weight (g/mol)."""
    return conc_mg_ml * 1e9 / mw


def um_to_mg_per_ml(conc_um, mw):
    """Convert micromolar to mg/mL given molecular weight (g/mol)."""
    return conc_um * mw / 1e6


def mg_per_ml_to_um(conc_mg_ml, mw):
    """Convert mg/mL to micromolar given molecular weight (g/mol)."""
    return conc_mg_ml * 1e6 / mw


def ml_min_kg_to_l_h(cl_ml_min_kg):
    """Convert clearance from mL/min/kg to L/h/kg."""
    return cl_ml_min_kg * 60 / 1000


def l_h_to_ml_min_kg(cl_l_h_kg):
    """Convert clearance from L/h/kg to mL/min/kg."""
    return cl_l_h_kg * 1000 / 60


def mg_kg_to_total_mg(dose_mg_kg, body_weight_kg):
    """Convert mg/kg dose to total mg for a given body weight."""
    return dose_mg_kg * body_weight_kg


def total_mg_to_mg_kg(dose_mg, body_weight_kg):
    """Convert total mg dose to mg/kg for a given body weight."""
    return dose_mg / body_weight_kg


def half_life_to_kel(t_half):
    """Convert half-life to elimination rate constant (same time units)."""
    return math.log(2) / t_half


def kel_to_half_life(kel):
    """Convert elimination rate constant to half-life (same time units)."""
    return math.log(2) / kel
