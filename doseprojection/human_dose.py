"""
Human dose projection calculations.

Includes Human Equivalent Dose (HED) from animal NOAEL, Maximum Recommended
Starting Dose (MRSD), and interspecies BSA-based dose conversion.

References:
- FDA Guidance (2005): Estimating the Maximum Safe Starting Dose
- Nair & Jacob (2016): Dose conversion between animals and human
"""

from doseprojection.constants import KM_FACTORS, BODY_WEIGHTS


def calc_hed(animal_dose_mg_kg, animal_species):
    """Calculate Human Equivalent Dose using FDA BSA normalization.

    HED (mg/kg) = Animal_dose (mg/kg) * (Animal_Km / Human_Km)

    Parameters
    ----------
    animal_dose_mg_kg : float
        Dose in the animal species (mg/kg).
    animal_species : str
        Species name (e.g., "rat", "mouse", "dog", "monkey").

    Returns
    -------
    dict
        Keys: hed_mg_kg, hed_total_mg (for 60 kg human),
              animal_km, human_km, conversion_factor
    """
    species = animal_species.lower()
    if species not in KM_FACTORS:
        raise ValueError(
            f"Unknown species '{species}'. Available: {list(KM_FACTORS.keys())}"
        )
    animal_km = KM_FACTORS[species]
    human_km = KM_FACTORS["human"]

    hed_mg_kg = animal_dose_mg_kg * (animal_km / human_km)
    hed_total_mg = hed_mg_kg * BODY_WEIGHTS["human"]

    return {
        "hed_mg_kg": hed_mg_kg,
        "hed_total_mg": hed_total_mg,
        "animal_km": animal_km,
        "human_km": human_km,
        "conversion_factor": animal_km / human_km,
    }


def calc_mrsd(hed_mg_kg, safety_factor=10):
    """Calculate Maximum Recommended Starting Dose from HED.

    MRSD = HED / Safety Factor

    Parameters
    ----------
    hed_mg_kg : float
        Human Equivalent Dose (mg/kg).
    safety_factor : float
        Safety factor (default: 10, per FDA guidance).

    Returns
    -------
    dict
        Keys: mrsd_mg_kg, mrsd_total_mg (for 60 kg human), safety_factor
    """
    mrsd_mg_kg = hed_mg_kg / safety_factor
    mrsd_total_mg = mrsd_mg_kg * BODY_WEIGHTS["human"]

    return {
        "mrsd_mg_kg": mrsd_mg_kg,
        "mrsd_total_mg": mrsd_total_mg,
        "safety_factor": safety_factor,
    }


def hed_from_noael(noael_mg_kg, species, safety_factor=10):
    """Calculate HED and MRSD from animal NOAEL in one step.

    Parameters
    ----------
    noael_mg_kg : float
        No Observed Adverse Effect Level (mg/kg) from toxicology study.
    species : str
        Animal species.
    safety_factor : float
        Safety factor for MRSD (default: 10).

    Returns
    -------
    dict
        Keys: noael_mg_kg, species, hed_mg_kg, hed_total_mg,
              mrsd_mg_kg, mrsd_total_mg, safety_factor
    """
    hed = calc_hed(noael_mg_kg, species)
    mrsd = calc_mrsd(hed["hed_mg_kg"], safety_factor)

    return {
        "noael_mg_kg": noael_mg_kg,
        "species": species,
        "hed_mg_kg": hed["hed_mg_kg"],
        "hed_total_mg": hed["hed_total_mg"],
        "mrsd_mg_kg": mrsd["mrsd_mg_kg"],
        "mrsd_total_mg": mrsd["mrsd_total_mg"],
        "safety_factor": safety_factor,
    }


def bsa_conversion(dose_mg_kg, from_species, to_species):
    """Convert dose between any two species using BSA normalization.

    Dose_to = Dose_from * (Km_from / Km_to)

    Parameters
    ----------
    dose_mg_kg : float
        Dose in the source species (mg/kg).
    from_species : str
        Source species.
    to_species : str
        Target species.

    Returns
    -------
    dict
        Keys: dose_mg_kg, from_species, to_species, from_km, to_km
    """
    from_sp = from_species.lower()
    to_sp = to_species.lower()

    for sp in [from_sp, to_sp]:
        if sp not in KM_FACTORS:
            raise ValueError(
                f"Unknown species '{sp}'. Available: {list(KM_FACTORS.keys())}"
            )

    from_km = KM_FACTORS[from_sp]
    to_km = KM_FACTORS[to_sp]
    converted_dose = dose_mg_kg * (from_km / to_km)

    return {
        "dose_mg_kg": converted_dose,
        "from_species": from_sp,
        "to_species": to_sp,
        "from_km": from_km,
        "to_km": to_km,
    }
