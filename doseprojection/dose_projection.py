"""
Core dose projection calculations.

Projects efficacious doses for preclinical studies and clinical candidates
based on in vitro potency (IC50) and pharmacokinetic parameters.

References:
- Smith et al. (2010): PK/PD in dose selection
- Lowe et al. (2007): Anticipation of human dose (AHD)
"""


def efficacious_dose(ic50_um, cl_mL_min_kg, fu, f, tau_h, coverage_multiple=1.0):
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
    f : float
        Oral bioavailability as fraction (0-1).
    tau_h : float
        Dosing interval in hours.
    coverage_multiple : float
        Multiple of IC50 to target (default: 1.0). Use higher values
        (e.g., 3-10) for targets requiring near-complete inhibition.

    Returns
    -------
    float
        Projected dose in mg/kg (assumes IC50 in ug/mL-equivalent units
        after conversion; see notes).

    Notes
    -----
    The IC50 is used as a concentration target. The dose is in mass/kg units
    consistent with the concentration and clearance units. For exact mass
    units, the IC50 should be converted to ug/mL (= mg/L) using MW:
        IC50_ug_mL = IC50_uM * MW / 1000
    Then dose will be in mg/kg.
    """
    if fu <= 0 or fu > 1:
        raise ValueError(f"fu must be between 0 (exclusive) and 1, got {fu}")
    if f <= 0 or f > 1:
        raise ValueError(f"F must be between 0 (exclusive) and 1, got {f}")

    # Convert CL to mL/h/kg
    cl_mL_h_kg = cl_mL_min_kg * 60
    # IC50 in ug/mL needs MW. Here we work in concentration-consistent units.
    # Dose (ug/kg) = IC50 (ug/mL) * CL (mL/h/kg) * tau (h) / (fu * F)
    target_conc = ic50_um * coverage_multiple  # uM target
    dose = (target_conc * cl_mL_h_kg * tau_h) / (fu * f)
    return dose


def efficacious_dose_mg_kg(ic50_nm, mw, cl_mL_min_kg, fu, f, tau_h,
                            coverage_multiple=1.0):
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
    f : float
        Oral bioavailability as fraction (0-1).
    tau_h : float
        Dosing interval in hours.
    coverage_multiple : float
        Multiple of IC50 to target (default: 1.0).

    Returns
    -------
    float
        Projected dose in mg/kg.
    """
    # Convert IC50 nM -> mg/L (= ug/mL)
    ic50_mg_L = ic50_nm * mw / 1e6
    target_conc_mg_L = ic50_mg_L * coverage_multiple

    cl_mL_h_kg = cl_mL_min_kg * 60
    # Dose (mg/kg) = conc (mg/L) * CL (mL/h/kg) * tau (h) / (fu * F) / 1000
    # Note: mg/L * mL/h/kg = mg*mL/(L*h*kg) = mg/(1000*h*kg) per unit time
    # Actually: mg/L = ug/mL, and CL in mL/h/kg:
    # ug/mL * mL/h/kg * h = ug/kg -> divide by 1000 for mg/kg
    dose_mg_kg = (target_conc_mg_L * cl_mL_h_kg * tau_h) / (fu * f * 1000)
    return dose_mg_kg


def dose_from_target_css(css_target_ng_mL, cl_mL_min_kg, f, tau_h):
    """Calculate dose to achieve a target average steady-state concentration.

    Dose (mg/kg) = (Css * CL * tau) / F

    Parameters
    ----------
    css_target_ng_mL : float
        Target average steady-state plasma concentration (ng/mL).
    cl_mL_min_kg : float
        Clearance (mL/min/kg).
    f : float
        Oral bioavailability as fraction (0-1).
    tau_h : float
        Dosing interval in hours.

    Returns
    -------
    float
        Dose in mg/kg.
    """
    cl_mL_h_kg = cl_mL_min_kg * 60
    # ng/mL * mL/h/kg * h = ng/kg -> /1e6 for mg/kg
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


def steady_state_css(dose_mg_kg, f, cl_mL_min_kg, tau_h):
    """Calculate average steady-state plasma concentration.

    Css,avg = (Dose * F) / (CL * tau)

    Parameters
    ----------
    dose_mg_kg : float
        Dose in mg/kg.
    f : float
        Oral bioavailability as fraction (0-1).
    cl_mL_min_kg : float
        Clearance (mL/min/kg).
    tau_h : float
        Dosing interval in hours.

    Returns
    -------
    float
        Average steady-state concentration in ng/mL.
    """
    cl_mL_h_kg = cl_mL_min_kg * 60
    # mg/kg * 1e6 ng/mg / (mL/h/kg * h) = ng/mL
    css = (dose_mg_kg * f * 1e6) / (cl_mL_h_kg * tau_h)
    return css


def project_animal_dose(ic50_nm, mw, cl_animal_mL_min_kg, fu_animal, f_animal,
                         tau_h, coverage_multiple=1.0):
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
    f_animal : float
        Animal oral bioavailability as fraction (0-1).
    tau_h : float
        Dosing interval in hours.
    coverage_multiple : float
        Multiple of IC50 to target at Css,avg (default: 1.0).

    Returns
    -------
    dict
        Keys: dose_mg_kg, target_cu_ng_mL, predicted_css_total_ng_mL
    """
    dose = efficacious_dose_mg_kg(
        ic50_nm, mw, cl_animal_mL_min_kg, fu_animal, f_animal,
        tau_h, coverage_multiple
    )
    # Target unbound concentration
    target_cu_ng_mL = ic50_nm * mw / 1000 * coverage_multiple  # nM -> ng/mL approx
    # Predicted total Css
    css_total = steady_state_css(dose, f_animal, cl_animal_mL_min_kg, tau_h)

    return {
        "dose_mg_kg": dose,
        "target_cu_ng_mL": target_cu_ng_mL,
        "predicted_css_total_ng_mL": css_total,
    }
