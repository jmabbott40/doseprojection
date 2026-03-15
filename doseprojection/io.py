"""
Data input/output for dose projection tool.

Supports CSV and Excel (.xlsx) file formats with automatic detection.
Validates required columns and flags missing data.
"""

import os
import warnings

import pandas as pd


# Minimum required columns for in vitro data
INVITRO_REQUIRED = ["compound_id", "IC50_nM"]

# Minimum required columns for PK data
PK_REQUIRED = ["compound_id", "species"]


def _read_file(path):
    """Read a CSV or Excel file based on extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        return pd.read_csv(path)
    elif ext in (".xlsx", ".xls"):
        return pd.read_excel(path)
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .csv or .xlsx")


def _validate_columns(df, required, file_desc):
    """Check that required columns are present, warn about missing ones."""
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"{file_desc} is missing required columns: {missing}. "
            f"Found columns: {list(df.columns)}"
        )


def load_invitro_data(path):
    """Load in vitro pharmacology data from CSV or Excel.

    Expected columns (minimum: compound_id, IC50_nM):
        compound_id         - Unique compound identifier
        target              - Pharmacological target name
        IC50_nM             - IC50 in nanomolar
        IC50_uM             - IC50 in micromolar (optional, calculated if absent)
        solubility_ug_mL    - Aqueous solubility (ug/mL)
        formulation         - Formulation/vehicle used for solubility
        fu_plasma_human     - Fraction unbound in human plasma
        fu_plasma_rat       - Fraction unbound in rat plasma
        fu_plasma_mouse     - Fraction unbound in mouse plasma
        microsomal_t_half_min_human - Microsomal stability half-life (min), human
        microsomal_t_half_min_rat   - Microsomal stability half-life (min), rat
        papp_caco2_cm_s     - Caco-2 apparent permeability (cm/s)
        papp_pampa_cm_s     - PAMPA apparent permeability (cm/s)
        MW                  - Molecular weight (g/mol)
        study_id            - Study identifier
        notes               - Free-text notes

    Parameters
    ----------
    path : str
        Path to CSV or Excel file.

    Returns
    -------
    pd.DataFrame
    """
    df = _read_file(path)
    _validate_columns(df, INVITRO_REQUIRED, "In vitro data file")

    # Auto-calculate IC50_uM if not provided
    if "IC50_uM" not in df.columns and "IC50_nM" in df.columns:
        df["IC50_uM"] = df["IC50_nM"] / 1000.0

    # Flag compounds with missing key data
    optional_cols = [
        "fu_plasma_human", "fu_plasma_rat", "solubility_ug_mL",
        "microsomal_t_half_min_human", "MW",
    ]
    for col in optional_cols:
        if col in df.columns:
            n_missing = df[col].isna().sum()
            if n_missing > 0:
                warnings.warn(
                    f"{n_missing} compound(s) missing '{col}' values.",
                    stacklevel=2,
                )

    return df


def load_pk_data(path):
    """Load PK data from CSV or Excel.

    Expected columns (minimum: compound_id, species):
        compound_id     - Unique compound identifier
        species         - Animal species (mouse, rat, dog, monkey)
        route           - Route of administration (IV, PO)
        dose_mg_kg      - Dose administered (mg/kg)
        CL_mL_min_kg    - Clearance (mL/min/kg)
        Vss_L_kg        - Volume of distribution at steady state (L/kg)
        F_pct           - Oral bioavailability (%)
        t_half_h        - Terminal half-life (hours)
        Cmax_ng_mL      - Peak plasma concentration (ng/mL)
        AUC_ng_h_mL     - Area under the curve (ng*h/mL)
        study_id        - Study identifier
        notes           - Free-text notes

    Parameters
    ----------
    path : str
        Path to CSV or Excel file.

    Returns
    -------
    pd.DataFrame
    """
    df = _read_file(path)
    _validate_columns(df, PK_REQUIRED, "PK data file")

    # Normalize species names to lowercase
    df["species"] = df["species"].str.lower().str.strip()

    # Convert F_pct to fraction if present
    if "F_pct" in df.columns:
        df["F_fraction"] = df["F_pct"] / 100.0

    return df
