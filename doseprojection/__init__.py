"""
doseprojection - Preclinical dose projection tool for animal studies and human dose estimation.

Provides functions for:
- In vitro to in vivo extrapolation (IVIVE) using the well-stirred model
- Allometric scaling of PK parameters across species
- Efficacious dose projection from IC50 and PK data
- Human equivalent dose (HED) and maximum recommended starting dose (MRSD)
- Absorption assessment (MAD, BCS classification)
"""

from doseprojection.ivive import (
    calc_clint_invitro,
    scale_clint,
    predict_hepatic_clearance,
    predict_extraction_ratio,
    predict_fh,
)
from doseprojection.allometry import (
    simple_allometry,
    predict_human_cl,
    predict_human_vss,
    predict_human_thalf,
)
from doseprojection.dose_projection import (
    efficacious_dose,
    dose_from_target_css,
    unbound_concentration,
    steady_state_css,
    project_animal_dose,
)
from doseprojection.human_dose import (
    calc_hed,
    calc_mrsd,
    hed_from_noael,
    bsa_conversion,
)
from doseprojection.absorption import (
    dose_number,
    max_absorbable_dose,
    classify_bcs,
    classify_permeability,
)

__version__ = "0.1.0"
