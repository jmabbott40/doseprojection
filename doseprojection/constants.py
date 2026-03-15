"""
Physiological constants and reference tables for dose projection calculations.

Sources:
- FDA Guidance (2005): "Estimating the Maximum Safe Starting Dose in Initial
  Clinical Trials for Therapeutics in Adult Healthy Volunteers"
- Barter et al. (2007): Scaling factors for IVIVE
- Nair & Jacob (2016): Dose conversion guide between animals and human
"""

# FDA Km factors (body weight / body surface area) for BSA-based dose conversion
# Km = BW (kg) / BSA (m^2)
KM_FACTORS = {
    "mouse": 3,
    "hamster": 5,
    "rat": 6,
    "guinea_pig": 8,
    "rabbit": 12,
    "monkey": 12,
    "dog": 20,
    "baboon": 20,
    "child": 25,
    "human": 37,
}

# Standard reference body weights (kg) used in allometric scaling
BODY_WEIGHTS = {
    "mouse": 0.02,
    "hamster": 0.08,
    "rat": 0.15,
    "guinea_pig": 0.4,
    "rabbit": 1.8,
    "monkey": 3.0,
    "dog": 10.0,
    "baboon": 12.0,
    "child": 20.0,
    "human": 60.0,
}

# Body surface area (m^2)
BSA = {
    "mouse": 0.007,
    "hamster": 0.02,
    "rat": 0.025,
    "guinea_pig": 0.05,
    "rabbit": 0.15,
    "monkey": 0.24,
    "dog": 0.50,
    "baboon": 0.60,
    "child": 0.80,
    "human": 1.62,
}

# Microsomal protein per gram of liver (mg protein / g liver)
MPPGL = {
    "mouse": 45,
    "rat": 45,
    "dog": 55,
    "monkey": 50,
    "human": 32,
}

# Hepatocellularity (10^6 cells / g liver)
HPGL = {
    "mouse": 120,
    "rat": 120,
    "dog": 120,
    "monkey": 120,
    "human": 99,
}

# Liver weight (g liver / kg body weight)
LIVER_WEIGHT = {
    "mouse": 55,
    "rat": 40,
    "dog": 32,
    "monkey": 30,
    "human": 21.4,
}

# Hepatic blood flow (mL/min/kg body weight)
HEPATIC_BLOOD_FLOW = {
    "mouse": 90,
    "rat": 55,
    "dog": 31,
    "monkey": 44,
    "human": 20.7,
}

# Standard allometric exponents for PK parameter scaling
ALLOMETRIC_EXPONENTS = {
    "clearance": 0.75,
    "volume": 1.0,
    "half_life": 0.25,
}

# Caco-2 permeability classification thresholds (Papp in cm/s)
PERMEABILITY_THRESHOLDS = {
    "low_upper": 1e-6,       # < 1e-6 = low permeability
    "moderate_upper": 10e-6,  # 1e-6 to 10e-6 = moderate
    # > 10e-6 = high permeability
}

# BCS solubility threshold: highest dose strength dissolves in 250 mL
# across pH 1.0-6.8 at 37C. Practical threshold often ~1 mg/mL for
# a 250 mg dose (dose/solubility ratio <= 250 mL).
BCS_GI_VOLUME_ML = 250

# Small intestinal transit time (minutes) for MAD calculation
SI_TRANSIT_TIME_MIN = 270
