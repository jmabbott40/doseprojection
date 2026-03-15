"""
Allometric scaling of pharmacokinetic parameters across species.

References:
- Mahmood & Balian (1996): Rule of exponents
- FDA Guidance (2005): BSA normalization
"""

import math

import numpy as np

from doseprojection.constants import BODY_WEIGHTS, ALLOMETRIC_EXPONENTS


def simple_allometry(animal_value, animal_bw_kg, human_bw_kg=60.0, exponent=0.75):
    """Scale a PK parameter from one species to another using simple allometry.

    Y_human = Y_animal * (BW_human / BW_animal)^exponent

    Parameters
    ----------
    animal_value : float
        PK parameter value in the animal species.
    animal_bw_kg : float
        Animal body weight in kg.
    human_bw_kg : float
        Human body weight in kg (default: 60).
    exponent : float
        Allometric exponent (default: 0.75 for clearance).

    Returns
    -------
    float
        Predicted human PK parameter value.
    """
    return animal_value * (human_bw_kg / animal_bw_kg) ** exponent


def predict_human_cl(animal_cl_data, human_bw_kg=60.0):
    """Predict human clearance from single or multi-species data.

    For single species: uses simple allometry with exponent 0.75.
    For multi-species (>=3): fits log-log regression to find a and b,
    then applies rule of exponents.

    Parameters
    ----------
    animal_cl_data : dict
        Dict of {species: CL_mL_min_kg} or {species: (CL_mL_min_kg, BW_kg)}.
        If BW not provided, standard weights from constants are used.
    human_bw_kg : float
        Human body weight in kg (default: 60).

    Returns
    -------
    dict
        Keys: cl_human_mL_min_kg, method, exponent, coefficient
    """
    # Parse input
    species_data = []
    for species, val in animal_cl_data.items():
        sp = species.lower()
        if isinstance(val, (tuple, list)):
            cl, bw = val
        else:
            cl = val
            if sp not in BODY_WEIGHTS:
                raise ValueError(f"No default body weight for '{sp}'. Provide (CL, BW) tuple.")
            bw = BODY_WEIGHTS[sp]
        # Convert CL from mL/min/kg to total mL/min for allometry
        cl_total = cl * bw
        species_data.append((sp, bw, cl_total))

    if len(species_data) == 1:
        sp, bw, cl_total = species_data[0]
        exponent = ALLOMETRIC_EXPONENTS["clearance"]
        human_cl_total = cl_total * (human_bw_kg / bw) ** exponent
        human_cl = human_cl_total / human_bw_kg
        return {
            "cl_human_mL_min_kg": human_cl,
            "method": "single_species_allometry",
            "exponent": exponent,
        }

    # Multi-species: log-log linear regression
    log_bw = np.array([math.log(d[1]) for d in species_data])
    log_cl = np.array([math.log(d[2]) for d in species_data])
    # Fit: log(CL) = log(a) + b * log(BW)
    b, log_a = np.polyfit(log_bw, log_cl, 1)
    a = math.exp(log_a)

    method = rule_of_exponents(b)
    human_cl_total = a * human_bw_kg ** b
    human_cl = human_cl_total / human_bw_kg

    return {
        "cl_human_mL_min_kg": human_cl,
        "method": method,
        "exponent": b,
        "coefficient": a,
    }


def predict_human_vss(animal_vss_data, human_bw_kg=60.0):
    """Predict human volume of distribution from animal data.

    Parameters
    ----------
    animal_vss_data : dict
        Dict of {species: Vss_L_kg} or {species: (Vss_L_kg, BW_kg)}.
    human_bw_kg : float
        Human body weight in kg (default: 60).

    Returns
    -------
    dict
        Keys: vss_human_L_kg, method, exponent
    """
    species_data = []
    for species, val in animal_vss_data.items():
        sp = species.lower()
        if isinstance(val, (tuple, list)):
            vss, bw = val
        else:
            vss = val
            if sp not in BODY_WEIGHTS:
                raise ValueError(f"No default body weight for '{sp}'. Provide (Vss, BW) tuple.")
            bw = BODY_WEIGHTS[sp]
        vss_total = vss * bw
        species_data.append((sp, bw, vss_total))

    if len(species_data) == 1:
        sp, bw, vss_total = species_data[0]
        exponent = ALLOMETRIC_EXPONENTS["volume"]
        human_vss_total = vss_total * (human_bw_kg / bw) ** exponent
        human_vss = human_vss_total / human_bw_kg
        return {
            "vss_human_L_kg": human_vss,
            "method": "single_species_allometry",
            "exponent": exponent,
        }

    log_bw = np.array([math.log(d[1]) for d in species_data])
    log_vss = np.array([math.log(d[2]) for d in species_data])
    b, log_a = np.polyfit(log_bw, log_vss, 1)
    a = math.exp(log_a)
    human_vss_total = a * human_bw_kg ** b
    human_vss = human_vss_total / human_bw_kg

    return {
        "vss_human_L_kg": human_vss,
        "method": "multi_species_allometry",
        "exponent": b,
        "coefficient": a,
    }


def predict_human_thalf(cl_human_mL_min_kg, vss_human_L_kg):
    """Predict human half-life from clearance and volume of distribution.

    t1/2 = (0.693 * Vss) / CL

    Parameters
    ----------
    cl_human_mL_min_kg : float
        Predicted human clearance (mL/min/kg).
    vss_human_L_kg : float
        Predicted human volume of distribution (L/kg).

    Returns
    -------
    float
        Predicted half-life in hours.
    """
    # Convert CL to L/h/kg for consistent units
    cl_L_h_kg = cl_human_mL_min_kg * 60 / 1000
    return (0.693 * vss_human_L_kg) / cl_L_h_kg


def rule_of_exponents(exponent):
    """Apply Mahmood & Balian rule of exponents to select scaling method.

    Parameters
    ----------
    exponent : float
        Observed allometric exponent from multi-species regression.

    Returns
    -------
    str
        Recommended scaling method.
    """
    if exponent < 0.55:
        return "simple_allometry (may overestimate)"
    elif exponent <= 0.70:
        return "simple_allometry"
    elif exponent <= 0.99:
        return "MLP_correction_recommended"
    else:
        return "brain_weight_correction_recommended"
