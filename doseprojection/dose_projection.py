"""
Core dose projection calculations.

Projects efficacious doses for preclinical studies and clinical candidates
based on in vitro potency (IC50) and pharmacokinetic parameters.

Supports oral (PO), intravenous (IV), and subcutaneous (SC) routes.
For IV, F=1.0 (complete bioavailability). For SC, F is the SC-specific
bioavailability. For PO, F is oral bioavailability.

References:
- Smith et al. (2010): PK/PD in dose selection
- Lowe et al. (2007): Anticipation of human dose (AHD)
"""

VALID_ROUTES = ("po", "iv", "sc")


def _resolve_f(f, route):
    """Resolve bioavailability based on route of administration.

    For IV, F is always 1.0 regardless of input.
    For PO and SC, F must be provided and valid.

    Parameters
    ----------
    f : float or None
        Bioavailability as fraction (0-1). Ignored for IV route.
    route : str
        Route of administration: "po", "iv", or "sc".

    Returns
    -------
    float
        Resolved bioavailability value.
    """
    route = route.lower()
    if route not in VALID_ROUTES:
        raise ValueError(f"Route must be one of {VALID_ROUTES}, got '{route}'")

    if route == "iv":
        return 1.0

    # PO and SC both require F
    if f is None:
        raise ValueError(f"Bioavailability (f) is required for route '{route}'")
    if f <= 0 or f > 1:
        raise ValueError(f"F must be between 0 (exclusive) and 1, got {f}")
    return f


def efficacious_dose(ic50_um, cl_mL_min_kg, fu, tau_h, f=None,
                     coverage_multiple=1.0, route="po"):
    """Project efficacious dose targeting unbound Css,avg >= IC50.

    Dose = (IC50 * coverage * CL * tau) / (fu * F)

    Parameters
    ----------
    ic50_um : float
        IC50 in micromolar (uM). This is the target unbound concentration.
    cl_mL_min_kg : float
        Clearance (mL/min/kg).
    fu : float
        Fraction unbound in plasma (0-1).
    tau_h : float
        Dosing interval in hours.
    f : float or None
        Bioavailability as fraction (0-1). Required for PO and SC routes.
        Ignored for IV (set to 1.0 automatically).
    coverage_multiple : float
        Multiple of IC50 to target (default: 1.0). Use higher values
        (e.g., 3-10) for targets requiring near-complete inhibition.
    route : str
        Route of administration: "po", "iv", or "sc" (default: "po").

    Returns
    -------
    float
        Projected dose in concentration-consistent units.

    Notes
    -----
    The IC50 is used as a concentration target. For dose in mg/kg,
    use efficacious_dose_mg_kg() which handles unit conversion via MW.
    """
    if fu <= 0 or fu > 1:
        raise ValueError(f"fu must be between 0 (exclusive) and 1, got {fu}")
    f = _resolve_f(f, route)

    cl_mL_h_kg = cl_mL_min_kg * 60
    target_conc = ic50_um * coverage_multiple
    dose = (target_conc * cl_mL_h_kg * tau_h) / (fu * f)
    return dose


def efficacious_dose_mg_kg(ic50_nm, mw, cl_mL_min_kg, fu, tau_h, f=None,
                            coverage_multiple=1.0, route="po"):
    """Project efficacious dose in mg/kg from IC50 in nM.

    Converts IC50 from nM to mg/L (= ug/mL), then calculates dose.

    Parameters
    ----------
    ic50_nm : float
        IC50 in nanomolar.
    mw : float
        Molecular weight (g/mol).
    cl_mL_min_kg : float
        Clearance (mL/min/kg).
    fu : float
        Fraction unbound in plasma (0-1).
    tau_h : float
        Dosing interval in hours.
    f : float or None
        Bioavailability as fraction (0-1). Required for PO and SC.
        Ignored for IV (set to 1.0 automatically).
    coverage_multiple : float
        Multiple of IC50 to target (default: 1.0).
    route : str
        Route of administration: "po", "iv", or "sc" (default: "po").

    Returns
    -------
    float
        Projected dose in mg/kg.
    """
    f = _resolve_f(f, route)

    ic50_mg_L = ic50_nm * mw / 1e6
    target_conc_mg_L = ic50_mg_L * coverage_multiple

    cl_mL_h_kg = cl_mL_min_kg * 60
    dose_mg_kg = (target_conc_mg_L * cl_mL_h_kg * tau_h) / (fu * f * 1000)
    return dose_mg_kg


def dose_from_target_css(css_target_ng_mL, cl_mL_min_kg, tau_h, f=None,
                          route="po"):
    """Calculate dose to achieve a target average steady-state concentration.

    Dose (mg/kg) = (Css * CL * tau) / F

    Parameters
    ----------
    css_target_ng_mL : float
        Target average steady-state plasma concentration (ng/mL).
    cl_mL_min_kg : float
        Clearance (mL/min/kg).
    tau_h : float
        Dosing interval in hours.
    f : float or None
        Bioavailability as fraction (0-1). Required for PO and SC.
        Ignored for IV (set to 1.0 automatically).
    route : str
        Route of administration: "po", "iv", or "sc" (default: "po").

    Returns
    -------
    float
        Dose in mg/kg.
    """
    f = _resolve_f(f, route)
    cl_mL_h_kg = cl_mL_min_kg * 60
    dose_mg_kg = (css_target_ng_mL * cl_mL_h_kg * tau_h) / (f * 1e6)
    return dose_mg_kg


def unbound_concentration(c_total, fu):
    """Calculate unbound (free) drug concentration.

    Cu = fu * C_total

    Parameters
    ----------
    c_total : float
        Total plasma concentration (any units).
    fu : float
        Fraction unbound in plasma (0-1).

    Returns
    -------
    float
        Unbound concentration (same units as c_total).
    """
    return fu * c_total


def steady_state_css(dose_mg_kg, cl_mL_min_kg, tau_h, f=None, route="po"):
    """Calculate average steady-state plasma concentration.

    Css,avg = (Dose * F) / (CL * tau)

    Parameters
    ----------
    dose_mg_kg : float
        Dose in mg/kg.
    cl_mL_min_kg : float
        Clearance (mL/min/kg).
    tau_h : float
        Dosing interval in hours.
    f : float or None
        Bioavailability as fraction (0-1). Required for PO and SC.
        Ignored for IV (set to 1.0 automatically).
    route : str
        Route of administration: "po", "iv", or "sc" (default: "po").

    Returns
    -------
    float
        Average steady-state concentration in ng/mL.
    """
    f = _resolve_f(f, route)
    cl_mL_h_kg = cl_mL_min_kg * 60
    css = (dose_mg_kg * f * 1e6) / (cl_mL_h_kg * tau_h)
    return css


def project_animal_dose(ic50_nm, mw, cl_animal_mL_min_kg, fu_animal, tau_h,
                         f_animal=None, coverage_multiple=1.0, route="po"):
    """Project dose for a preclinical animal study.

    Convenience wrapper for efficacious_dose_mg_kg using animal-specific
    PK parameters. Used for setting dose levels in efficacy studies.

    Parameters
    ----------
    ic50_nm : float
        IC50 in nanomolar.
    mw : float
        Molecular weight (g/mol).
    cl_animal_mL_min_kg : float
        Animal clearance (mL/min/kg).
    fu_animal : float
        Fraction unbound in animal plasma (0-1).
    tau_h : float
        Dosing interval in hours.
    f_animal : float or None
        Animal bioavailability as fraction (0-1). Required for PO and SC.
        Ignored for IV (set to 1.0 automatically).
    coverage_multiple : float
        Multiple of IC50 to target at Css,avg (default: 1.0).
    route : str
        Route of administration: "po", "iv", or "sc" (default: "po").

    Returns
    -------
    dict
        Keys: dose_mg_kg, target_cu_ng_mL, predicted_css_total_ng_mL, route
    """
    dose = efficacious_dose_mg_kg(
        ic50_nm, mw, cl_animal_mL_min_kg, fu_animal, tau_h,
        f=f_animal, coverage_multiple=coverage_multiple, route=route
    )
    target_cu_ng_mL = ic50_nm * mw / 1000 * coverage_multiple
    f_resolved = _resolve_f(f_animal, route)
    css_total = steady_state_css(dose, cl_animal_mL_min_kg, tau_h,
                                  f=f_resolved, route=route)

    return {
        "dose_mg_kg": dose,
        "target_cu_ng_mL": target_cu_ng_mL,
        "predicted_css_total_ng_mL": css_total,
        "route": route.lower(),
    }
