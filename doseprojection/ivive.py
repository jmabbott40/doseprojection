"""
In Vitro to In Vivo Extrapolation (IVIVE) using the well-stirred model.

Predicts hepatic clearance from microsomal stability data.

References:
- Obach (1999): Prediction of human clearance from hepatic microsomal data
- Barter et al. (2007): Scaling factors for IVIVE
"""

from doseprojection.constants import MPPGL, LIVER_WEIGHT, HEPATIC_BLOOD_FLOW


def calc_clint_invitro(t_half_min, volume_ul=1000, protein_mg=0.5):
    """Calculate in vitro intrinsic clearance from substrate depletion.

    CLint = (0.693 / t1/2) * (incubation volume / mg protein)

    Parameters
    ----------
    t_half_min : float
        In vitro half-life from microsomal stability assay (minutes).
    volume_ul : float
        Incubation volume in microliters (default: 1000 uL).
    protein_mg : float
        Microsomal protein in the incubation (default: 0.5 mg).

    Returns
    -------
    float
        Intrinsic clearance in uL/min/mg protein.
    """
    if t_half_min <= 0:
        raise ValueError("Half-life must be positive.")
    return (0.693 / t_half_min) * (volume_ul / protein_mg)


def scale_clint(clint_invitro, species="human"):
    """Scale in vitro CLint to whole-liver in vivo CLint.

    CLint_scaled = CLint_invitro * MPPGL * liver_weight

    Parameters
    ----------
    clint_invitro : float
        In vitro intrinsic clearance (uL/min/mg protein).
    species : str
        Species name (default: "human").

    Returns
    -------
    float
        Scaled intrinsic clearance in mL/min/kg body weight.
    """
    species = species.lower()
    if species not in MPPGL:
        raise ValueError(
            f"Unknown species '{species}'. Available: {list(MPPGL.keys())}"
        )
    # uL/min/mg * mg/g * g/kg = uL/min/kg -> convert to mL/min/kg
    return clint_invitro * MPPGL[species] * LIVER_WEIGHT[species] / 1000.0


def predict_hepatic_clearance(clint_scaled, fu, species="human"):
    """Predict hepatic clearance using the well-stirred model.

    CLh = (QH * fu * CLint) / (QH + fu * CLint)

    Parameters
    ----------
    clint_scaled : float
        Scaled in vivo intrinsic clearance (mL/min/kg).
    fu : float
        Fraction unbound in blood (0-1).
    species : str
        Species name (default: "human").

    Returns
    -------
    float
        Predicted hepatic clearance (mL/min/kg).
    """
    species = species.lower()
    if species not in HEPATIC_BLOOD_FLOW:
        raise ValueError(
            f"Unknown species '{species}'. "
            f"Available: {list(HEPATIC_BLOOD_FLOW.keys())}"
        )
    qh = HEPATIC_BLOOD_FLOW[species]
    return (qh * fu * clint_scaled) / (qh + fu * clint_scaled)


def predict_extraction_ratio(clint_scaled, fu, species="human"):
    """Calculate hepatic extraction ratio from the well-stirred model.

    EH = (fu * CLint) / (QH + fu * CLint)

    Parameters
    ----------
    clint_scaled : float
        Scaled in vivo intrinsic clearance (mL/min/kg).
    fu : float
        Fraction unbound in blood (0-1).
    species : str
        Species name (default: "human").

    Returns
    -------
    float
        Hepatic extraction ratio (0-1).
    """
    species = species.lower()
    qh = HEPATIC_BLOOD_FLOW[species]
    return (fu * clint_scaled) / (qh + fu * clint_scaled)


def predict_fh(extraction_ratio):
    """Calculate hepatic bioavailability from extraction ratio.

    Fh = 1 - EH

    Parameters
    ----------
    extraction_ratio : float
        Hepatic extraction ratio (0-1).

    Returns
    -------
    float
        Fraction escaping hepatic first-pass metabolism (0-1).
    """
    return 1.0 - extraction_ratio


def ivive_workflow(t_half_min, fu, species="human",
                   volume_ul=1000, protein_mg=0.5):
    """Run full IVIVE workflow from microsomal t1/2 to hepatic clearance.

    Parameters
    ----------
    t_half_min : float
        Microsomal stability half-life (minutes).
    fu : float
        Fraction unbound in plasma (0-1).
    species : str
        Species (default: "human").
    volume_ul : float
        Incubation volume (uL).
    protein_mg : float
        Microsomal protein amount (mg).

    Returns
    -------
    dict
        Keys: clint_invitro, clint_scaled, cl_hepatic, extraction_ratio, fh
    """
    clint_iv = calc_clint_invitro(t_half_min, volume_ul, protein_mg)
    clint_sc = scale_clint(clint_iv, species)
    eh = predict_extraction_ratio(clint_sc, fu, species)
    cl_h = predict_hepatic_clearance(clint_sc, fu, species)
    fh = predict_fh(eh)

    return {
        "clint_invitro_uL_min_mg": clint_iv,
        "clint_scaled_mL_min_kg": clint_sc,
        "cl_hepatic_mL_min_kg": cl_h,
        "extraction_ratio": eh,
        "fh": fh,
    }
