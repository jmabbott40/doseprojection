# doseprojection

A Python toolkit for preclinical dose projection in drug discovery. Estimates efficacious doses for animal studies and projects human equivalent doses using established pharmacokinetic methods.

## Features

- **IVIVE (In Vitro to In Vivo Extrapolation)**: Predict hepatic clearance from microsomal stability data using the well-stirred model
- **Allometric Scaling**: Scale PK parameters across species using simple allometry and the rule of exponents
- **Dose Projection from IC50**: Estimate efficacious doses from in vitro potency, clearance, protein binding, and bioavailability
- **Human Dose Estimation**: Calculate Human Equivalent Dose (HED) and Maximum Recommended Starting Dose (MRSD) from animal NOAEL
- **Absorption Assessment**: Maximum Absorbable Dose (MAD), BCS classification, permeability classification

## Installation

```bash
git clone https://github.com/jmabbott40/doseprojection.git
cd doseprojection
pip install -e .
```

## Quick Start

```python
from doseprojection import ivive, dose_projection, human_dose, absorption

# --- IVIVE: Predict hepatic clearance from microsomal t1/2 ---
result = ivive.ivive_workflow(t_half_min=30, fu=0.1, species="human")
print(f"Predicted human CLh: {result['cl_hepatic_mL_min_kg']:.1f} mL/min/kg")
print(f"Hepatic bioavailability (Fh): {result['fh']:.2f}")

# --- Dose projection from IC50 ---
dose = dose_projection.efficacious_dose_mg_kg(
    ic50_nm=100,       # IC50 in nM
    mw=400,            # Molecular weight
    cl_mL_min_kg=25,   # Rat clearance
    fu=0.1,            # Fraction unbound
    f=0.5,             # Oral bioavailability
    tau_h=24,          # Once daily dosing
    coverage_multiple=1.0
)
print(f"Projected efficacious dose: {dose:.1f} mg/kg")

# --- Human equivalent dose from rat NOAEL ---
result = human_dose.hed_from_noael(noael_mg_kg=50, species="rat")
print(f"HED: {result['hed_mg_kg']:.1f} mg/kg")
print(f"MRSD: {result['mrsd_mg_kg']:.2f} mg/kg ({result['mrsd_total_mg']:.0f} mg total)")

# --- Check absorption feasibility ---
d0 = absorption.dose_number(dose_mg=500, solubility_mg_mL=0.1)
print(f"Dose number: {d0['dose_number']:.1f} — {d0['interpretation']}")
```

## Modules

| Module | Description |
|---|---|
| `constants.py` | FDA Km factors, MPPGL, liver weights, hepatic blood flow, BCS thresholds |
| `io.py` | Load in vitro and PK data from CSV/Excel with validation |
| `ivive.py` | Microsomal CLint → hepatic clearance (well-stirred model) |
| `allometry.py` | Allometric scaling, rule of exponents |
| `dose_projection.py` | IC50-based and PK-based dose projection |
| `human_dose.py` | HED from NOAEL, MRSD, interspecies BSA conversion |
| `absorption.py` | MAD, dose number, BCS/permeability classification |
| `utils.py` | Unit conversions (nM↔uM↔mg/mL, etc.) |

## Key Equations

### Efficacious Dose from IC50
```
Dose (mg/kg) = (IC50 × CL × τ) / (fu × F)
```
Where IC50 is the target unbound concentration, CL is clearance, τ is dosing interval, fu is fraction unbound, and F is oral bioavailability.

### IVIVE — Well-Stirred Model
```
CLint,invitro = (0.693 / t½) × (Volume / mg protein)
CLint,scaled  = CLint,invitro × MPPGL × Liver Weight
CLh           = (QH × fu × CLint) / (QH + fu × CLint)
```

### Human Equivalent Dose (FDA BSA Method)
```
HED (mg/kg) = Animal Dose (mg/kg) × (Animal Km / Human Km)
MRSD        = HED / Safety Factor (default: 10)
```

### Allometric Scaling
```
Y_human = a × BW^b
```
Standard exponents: CL (b=0.75), Vss (b=1.0), t½ (b=0.25).

### Maximum Absorbable Dose
```
MAD = Solubility × ka × SIWV × Transit Time
```

## Data Input Format

### In Vitro Data (CSV/Excel)
Required: `compound_id`, `IC50_nM`

Optional: `target`, `solubility_ug_mL`, `formulation`, `fu_plasma_human`, `fu_plasma_rat`, `fu_plasma_mouse`, `microsomal_t_half_min_human`, `microsomal_t_half_min_rat`, `papp_caco2_cm_s`, `papp_pampa_cm_s`, `MW`, `study_id`, `notes`

### PK Data (CSV/Excel)
Required: `compound_id`, `species`

Optional: `route`, `dose_mg_kg`, `CL_mL_min_kg`, `Vss_L_kg`, `F_pct`, `t_half_h`, `Cmax_ng_mL`, `AUC_ng_h_mL`, `study_id`, `notes`

See `examples/` for sample data files.

## Jupyter Notebooks

Interactive examples in `notebooks/`:

1. **01_data_input_template** — Loading and validating data
2. **02_ivive_clearance** — IVIVE workflow from microsomal t½ to CLh
3. **03_dose_projection** — IC50-based dose projection for animal studies
4. **04_human_dose_estimation** — HED and MRSD from animal NOAEL
5. **05_full_workflow** — End-to-end dose projection pipeline

## Physiological Constants

Species-specific constants from published references:

| Species | Km | BW (kg) | MPPGL (mg/g) | Liver (g/kg) | QH (mL/min/kg) |
|---|---|---|---|---|---|
| Mouse | 3 | 0.02 | 45 | 55 | 90 |
| Rat | 6 | 0.15 | 45 | 40 | 55 |
| Dog | 20 | 10.0 | 55 | 32 | 31 |
| Monkey | 12 | 3.0 | 50 | 30 | 44 |
| Human | 37 | 60.0 | 32 | 21.4 | 20.7 |

## References

1. **FDA (2005)**. Guidance for Industry: Estimating the Maximum Safe Starting Dose in Initial Clinical Trials for Therapeutics in Adult Healthy Volunteers. CDER/CBER.
2. **Nair AB, Jacob S (2016)**. A simple practice guide for dose conversion between animals and human. *J Basic Clin Pharma* 7(2):27-31. [PMC4804402](https://pmc.ncbi.nlm.nih.gov/articles/PMC4804402/)
3. **Mahmood I, Balian JD (1996)**. Interspecies scaling: predicting clearance of drugs in humans. Three different approaches. *Xenobiotica* 26(9):887-895.
4. **Obach RS (1999)**. Prediction of human clearance of twenty-nine drugs from hepatic microsomal intrinsic clearance data. *Drug Metab Dispos* 27(11):1350-1359.
5. **Barter ZE et al. (2007)**. Scaling factors for the extrapolation of in vivo metabolic drug clearance from in vitro data. *Curr Drug Metab* 8(1):33-45. [PMC1884378](https://pmc.ncbi.nlm.nih.gov/articles/PMC1884378/)
6. **Smith DA, Di L, Kerns EH (2010)**. The effect of plasma protein binding on in vivo efficacy: misconceptions in drug discovery. *Nat Rev Drug Discov* 9(12):929-939.
7. **Johnson KC, Swindell AC (1996)**. Guidance in the setting of drug particle size specifications to minimize variability in absorption. *Pharm Res* 13(12):1795-1798.
8. **Butler JM, Dressman JB (2010)**. The developability classification system: application of biopharmaceutics concepts to formulation development. *J Pharm Sci* 99(12):4940-4954.
9. **Lowe PJ et al. (2007)**. On the anticipation of the human dose in first-in-man trials from preclinical and prior clinical information in early drug development. *Xenobiotica* 37(10-11):1331-1354. [PMC2758129](https://pmc.ncbi.nlm.nih.gov/articles/PMC2758129/)

### Additional Recommended Reading

- **IVIVE underprediction**: [Riley et al. (2005)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7325626/) — systematic underprediction of clearance and empirical scaling factors
- **Allometric scaling review**: [Huh et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4181675/) — comparison of interspecies scaling methods
- **Plasma protein binding**: [Bohnert & Gan (2013)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6609267/) — PPB in drug discovery context
- **BCS classification**: [Larregieu & Benet (2014)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2844511/) — permeability classification
- **MAD and formulation**: [Curatolo (2012)](https://pubs.acs.org/doi/10.1021/mp200452r) — MAD model for formulation decisions
- **Human PK prediction**: [Jones et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3326168/) — FIH dose estimation from preclinical data
- **IC50 to clinical dose**: [Bhatt & Bhatt (2020)](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.1846) — correlation of in vitro potency to clinical concentrations

### Areas Where Additional Information Would Be Helpful

The following data/context improves dose projection accuracy when available:

- **Hepatocyte stability data** (in addition to microsomes) — reduces IVIVE underprediction ~2-fold
- **Blood-to-plasma ratio** — needed for accurate well-stirred model (fu,blood vs fu,plasma)
- **CYP reaction phenotyping** — identifies metabolic pathways for species-specific scaling
- **Transporter substrate assessment** (P-gp, BCRP, OATP) — affects Fa and tissue distribution
- **In vivo PK from a second species** (e.g., dog or monkey) — enables multi-species allometry
- **Pharmacodynamic biomarker data** — allows PK/PD-based dose projection rather than IC50 alone
- **Dose-response data from efficacy models** — provides empirical target exposure
- **Formulation composition details** — critical for compounds with low solubility (BCS II/IV)
- **Food effect assessment** — can significantly alter oral absorption

## Running Tests

```bash
pip install pytest
python -m pytest tests/ -v
```

## License

MIT
