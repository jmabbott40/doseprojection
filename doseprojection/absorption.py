"""
Absorption assessment: solubility-limited absorption, permeability
classification, BCS/DCS classification, and maximum absorbable dose (MAD).

References:
- Johnson & Swindell (1996): Maximum absorbable dose
- Amidon et al. (1995): BCS classification
- Butler & Dressman (2010): Developability Classification System
"""

from doseprojection.constants import (
    PERMEABILITY_THRESHOLDS,
    BCS_GI_VOLUME_ML,
    SI_TRANSIT_TIME_MIN,
)


def dose_number(dose_mg, solubility_mg_mL, volume_mL=BCS_GI_VOLUME_ML):
    """Calculate dimensionless dose number.

    D0 = Dose / (Solubility * Volume)

    D0 <= 1 indicates the dose can fully dissolve in GI fluid.
    D0 > 1 suggests solubility-limited absorption is possible.

    Parameters
    ----------
    dose_mg : float
        Dose in mg.
    solubility_mg_mL : float
        Aqueous solubility at intestinal pH (mg/mL).
    volume_mL : float
        GI fluid volume (default: 250 mL per FDA/BCS convention).

    Returns
    -------
    dict
        Keys: dose_number, solubility_limited (bool), interpretation
    """
    if solubility_mg_mL <= 0:
        raise ValueError("Solubility must be positive.")

    d0 = dose_mg / (solubility_mg_mL * volume_mL)
    limited = d0 > 1.0

    interpretation = (
        "Dose exceeds solubility capacity; enabling formulation may be needed"
        if limited
        else "Dose can dissolve in GI fluid volume"
    )

    return {
        "dose_number": d0,
        "solubility_limited": limited,
        "interpretation": interpretation,
    }


def max_absorbable_dose(solubility_mg_mL, ka_per_min,
                         siwv_mL=BCS_GI_VOLUME_ML,
                         transit_time_min=SI_TRANSIT_TIME_MIN):
    """Calculate Maximum Absorbable Dose (MAD).

    MAD = S * ka * SIWV * T

    Parameters
    ----------
    solubility_mg_mL : float
        Drug solubility at intestinal pH (mg/mL).
    ka_per_min : float
        Absorption rate constant (1/min). Can be estimated from Peff.
    siwv_mL : float
        Small intestinal water volume (default: 250 mL).
    transit_time_min : float
        Small intestinal transit time (default: 270 min).

    Returns
    -------
    dict
        Keys: mad_mg, solubility_mg_mL, ka_per_min, transit_time_min
    """
    mad = solubility_mg_mL * ka_per_min * siwv_mL * transit_time_min
    return {
        "mad_mg": mad,
        "solubility_mg_mL": solubility_mg_mL,
        "ka_per_min": ka_per_min,
        "transit_time_min": transit_time_min,
    }


def classify_permeability(papp_cm_s):
    """Classify permeability from Caco-2 Papp value.

    Parameters
    ----------
    papp_cm_s : float
        Apparent permeability coefficient (cm/s).

    Returns
    -------
    str
        "low", "moderate", or "high"
    """
    if papp_cm_s < PERMEABILITY_THRESHOLDS["low_upper"]:
        return "low"
    elif papp_cm_s < PERMEABILITY_THRESHOLDS["moderate_upper"]:
        return "moderate"
    else:
        return "high"


def classify_bcs(solubility_mg_mL, papp_cm_s, dose_mg=None):
    """Classify compound according to BCS (Biopharmaceutics Classification System).

    Uses permeability thresholds and solubility relative to dose.
    If dose is provided, high solubility means dose dissolves in 250 mL.
    If dose is not provided, a threshold of 0.1 mg/mL is used as a rough guide.

    Parameters
    ----------
    solubility_mg_mL : float
        Aqueous solubility (mg/mL).
    papp_cm_s : float
        Caco-2 apparent permeability (cm/s).
    dose_mg : float, optional
        Highest anticipated dose (mg). Used for dose-dependent
        solubility classification.

    Returns
    -------
    dict
        Keys: bcs_class, solubility_class, permeability_class, interpretation
    """
    perm_class = classify_permeability(papp_cm_s)
    high_perm = perm_class in ("moderate", "high")

    if dose_mg is not None:
        high_sol = (dose_mg / solubility_mg_mL) <= BCS_GI_VOLUME_ML
    else:
        high_sol = solubility_mg_mL >= 0.1  # rough threshold

    if high_sol and high_perm:
        bcs = "I"
        interp = "High solubility, high permeability — well absorbed"
    elif not high_sol and high_perm:
        bcs = "II"
        interp = "Low solubility, high permeability — dissolution rate-limited"
    elif high_sol and not high_perm:
        bcs = "III"
        interp = "High solubility, low permeability — permeability-limited"
    else:
        bcs = "IV"
        interp = "Low solubility, low permeability — poor oral candidate"

    return {
        "bcs_class": bcs,
        "solubility_class": "high" if high_sol else "low",
        "permeability_class": perm_class,
        "interpretation": interp,
    }


def predict_fa(papp_cm_s):
    """Estimate fraction absorbed from Caco-2 permeability.

    Uses a simplified sigmoidal relationship. This is an approximation;
    PBPK models provide more accurate predictions.

    Parameters
    ----------
    papp_cm_s : float
        Caco-2 Papp (cm/s).

    Returns
    -------
    float
        Estimated fraction absorbed (0-1).
    """
    # Simplified logistic model calibrated to general BCS boundaries
    # Fa ~ 1 / (1 + 10^(log10(2e-6) - log10(Papp)))
    import math
    if papp_cm_s <= 0:
        return 0.0
    log_ratio = math.log10(2e-6) - math.log10(papp_cm_s)
    fa = 1.0 / (1.0 + 10 ** log_ratio)
    return min(fa, 1.0)
