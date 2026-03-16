# doseprojection

A Python toolkit for preclinical dose projection in drug discovery. Estimates efficacious doses for animal studies and projects human equivalent doses using established pharmacokinetic methods. Supports **oral (PO)**, **intravenous (IV)**, and **subcutaneous (SC)** routes of administration.

## Features

- **Multi-Route Dose Projection**: Project doses for PO, IV, and SC routes with route-appropriate bioavailability handling
- **IVIVE (In Vitro to In Vivo Extrapolation)**: Predict hepatic clearance from microsomal stability data using the well-stirred model
- **Allometric Scaling**: Scale PK parameters across species using simple allometry and the rule of exponents
- **Dose Projection from IC50**: Estimate efficacious doses from in vitro potency, clearance, protein binding, and bioavailability
- **Human Dose Estimation**: Calculate Human Equivalent Dose (HED) and Maximum Recommended Starting Dose (MRSD) from animal NOAEL
- **Absorption Assessment**: Maximum Absorbable Dose (MAD), BCS classification, permeability classification (PO route only)

## Installation

```bash
git clone https://github.com/jmabbott40/doseprojection.git
cd doseprojection
pip install -e .
```

## Command-Line Usage (Step by Step)

For users who want to run projections without writing Python code. All you need are two CSV files (in vitro data and PK data) — see [Data Input Format](#data-input-format) below for column requirements.

### Step 1: Install

Open your terminal (Mac: Terminal app; Windows: Command Prompt or PowerShell) and run:

```bash
git clone https://github.com/jmabbott40/doseprojection.git
cd doseprojection
pip install -e .
```

If `pip` is not found, try `pip3 install -e .` instead.

### Step 2: Run with Example Data

Test that everything works using the included example files:

```bash
python run_projection.py examples/example_invitro_data.csv examples/example_pk_data.csv
```

This projects oral (PO) doses for all compounds in the example data using rat PK, once-daily dosing, at 1x IC50 coverage.

### Step 3: Run with Your Own Data

Replace the file paths with your own CSV files:

```bash
python run_projection.py my_invitro_data.csv my_pk_data.csv
```

### Step 4: Customize the Run

Add flags to change route, species, dosing interval, or coverage:

```bash
# IV dosing (no bioavailability needed):
python run_projection.py my_invitro.csv my_pk.csv --route iv

# Subcutaneous dosing:
python run_projection.py my_invitro.csv my_pk.csv --route sc

# BID dosing (every 12 hours) with 3x IC50 coverage:
python run_projection.py my_invitro.csv my_pk.csv --tau 12 --coverage 3

# Include HED/MRSD from a 50 mg/kg NOAEL:
python run_projection.py my_invitro.csv my_pk.csv --noael 50

# Mouse instead of rat:
python run_projection.py my_invitro.csv my_pk.csv --species mouse
```

### Step 5: Save Results to a File

Add `--output` (or `-o`) to save a CSV you can open in Excel:

```bash
python run_projection.py my_invitro.csv my_pk.csv --output results.csv
```

### All Options

| Flag | Default | Description |
|---|---|---|
| `--route` | `po` | Route of administration: `po`, `iv`, or `sc` |
| `--species` | `rat` | Species for PK data: `rat`, `mouse`, `dog`, `monkey` |
| `--tau` | `24` | Dosing interval in hours (24=QD, 12=BID, 8=TID) |
| `--coverage` | `1.0` | IC50 coverage multiple (e.g., 3 for 3x IC50) |
| `--noael` | *(none)* | NOAEL in mg/kg — adds HED/MRSD to output |
| `--output` / `-o` | *(none)* | Save results to this CSV file |
| `--quiet` / `-q` | off | Suppress warning messages |

Run `python run_projection.py --help` for the full help text.

## Quick Start (Python API)

```python
from doseprojection import ivive, dose_projection, human_dose, absorption

# --- IVIVE: Predict hepatic clearance from microsomal t1/2 ---
result = ivive.ivive_workflow(t_half_min=30, fu=0.1, species="human")
print(f"Predicted human CLh: {result['cl_hepatic_mL_min_kg']:.1f} mL/min/kg")
print(f"Hepatic bioavailability (Fh): {result['fh']:.2f}")

# --- Oral dose projection from IC50 ---
dose_po = dose_projection.efficacious_dose_mg_kg(
    ic50_nm=100,       # IC50 in nM
    mw=400,            # Molecular weight
    cl_mL_min_kg=25,   # Clearance
    fu=0.1,            # Fraction unbound
    tau_h=24,          # Once daily dosing
    f=0.5,             # Oral bioavailability
    route="po"         # Route of administration
)
print(f"Projected PO dose: {dose_po:.1f} mg/kg")

# --- IV dose projection (F is automatically 1.0) ---
dose_iv = dose_projection.efficacious_dose_mg_kg(
    ic50_nm=100, mw=400, cl_mL_min_kg=25, fu=0.1, tau_h=24,
    route="iv"         # No F needed — set to 1.0 automatically
)
print(f"Projected IV dose: {dose_iv:.1f} mg/kg")

# --- SC dose projection (uses SC-specific bioavailability) ---
dose_sc = dose_projection.efficacious_dose_mg_kg(
    ic50_nm=100, mw=400, cl_mL_min_kg=25, fu=0.1, tau_h=24,
    f=0.7,             # SC bioavailability (not oral F)
    route="sc"
)
print(f"Projected SC dose: {dose_sc:.1f} mg/kg")

# --- Human equivalent dose from rat NOAEL ---
result = human_dose.hed_from_noael(noael_mg_kg=50, species="rat")
print(f"HED: {result['hed_mg_kg']:.1f} mg/kg")
print(f"MRSD: {result['mrsd_mg_kg']:.2f} mg/kg ({result['mrsd_total_mg']:.0f} mg total)")

# --- Check absorption feasibility (PO only) ---
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
| `dose_projection.py` | IC50-based and PK-based dose projection (PO, IV, SC routes) |
| `human_dose.py` | HED from NOAEL, MRSD, interspecies BSA conversion |
| `absorption.py` | MAD, dose number, BCS/permeability classification (PO route only) |
| `utils.py` | Unit conversions (nM↔uM↔mg/mL, etc.) |

## Key Equations

### Efficacious Dose from IC50
```
Dose (mg/kg) = (IC50 × CL × τ) / (fu × F)
```
Where IC50 is the target unbound concentration, CL is clearance, τ is dosing interval, fu is fraction unbound, and F is route-dependent bioavailability (see table below).

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

### Maximum Absorbable Dose (PO only)
```
MAD = Solubility × ka × SIWV × Transit Time
```

## Route of Administration

All dose projection functions accept a `route` parameter (`"po"`, `"iv"`, or `"sc"`). The core equation is the same across routes — only the bioavailability term (F) changes.

### How F is handled per route

| Route | `route=` | F parameter | What F represents |
|---|---|---|---|
| **Oral (PO)** | `"po"` | Required | Oral bioavailability (Fa × Fg × Fh) |
| **Intravenous (IV)** | `"iv"` | Ignored (auto = 1.0) | Complete bioavailability by definition |
| **Subcutaneous (SC)** | `"sc"` | Required | SC bioavailability (typically 0.5–1.0) |

### What applies to each route

| Module / Concept | PO | IV | SC |
|---|---|---|---|
| Dose projection (`dose_projection.py`) | ✅ | ✅ | ✅ |
| IVIVE — clearance prediction (`ivive.py`) | ✅ | ✅ | ✅ |
| Allometric scaling (`allometry.py`) | ✅ | ✅ | ✅ |
| HED / MRSD (`human_dose.py`) | ✅ | ✅ | ✅ |
| Protein binding / fu | ✅ | ✅ | ✅ |
| Absorption: MAD, dose number (`absorption.py`) | ✅ | ❌ Not applicable | ❌ Not applicable |
| BCS / permeability classification | ✅ | ❌ Not applicable | ❌ Not applicable |
| Caco-2 / PAMPA permeability data | ✅ | ❌ Not applicable | ❌ Not applicable |

### Route-specific usage examples

```python
from doseprojection.dose_projection import project_animal_dose

# PO — requires oral bioavailability
po = project_animal_dose(
    ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
    fu_animal=0.1, tau_h=24,
    f_animal=0.45, route="po"
)

# IV — F is automatically 1.0, no f_animal needed
iv = project_animal_dose(
    ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
    fu_animal=0.1, tau_h=24,
    route="iv"
)

# SC — requires SC-specific bioavailability
sc = project_animal_dose(
    ic50_nm=100, mw=400, cl_animal_mL_min_kg=25,
    fu_animal=0.1, tau_h=24,
    f_animal=0.70, route="sc"
)

# IV always gives the lowest dose (F=1 means no absorption loss)
print(f"IV: {iv['dose_mg_kg']:.1f} mg/kg")
print(f"SC: {sc['dose_mg_kg']:.1f} mg/kg")
print(f"PO: {po['dose_mg_kg']:.1f} mg/kg")
```

### Notes on parenteral routes

- **IV**: Clearance (CL) from IV PK studies is the systemic clearance. No absorption phase exists, so Caco-2 permeability, solubility-limited absorption, MAD, and BCS classification are not relevant. Solubility still matters for formulation (solution concentration for injection volume), but this is a formulation constraint rather than an absorption assessment.
- **SC**: Bioavailability (F_sc) must be determined from SC PK studies (comparing SC AUC to IV AUC). SC absorption is from the injection site into systemic circulation — it is not GI absorption, so Caco-2/PAMPA permeability and BCS classification do not apply. SC absorption rate can affect Cmax and time-to-peak but the steady-state dose equation remains the same.
- **Clearance, protein binding, IVIVE, allometric scaling, and HED/MRSD** are all route-independent — they describe drug disposition after the drug reaches systemic circulation.

## Data Input Format

Two input files are used: one for in vitro data (one row per compound) and one for PK data (one row per compound × species × route).

### In Vitro Data (CSV/Excel) — one row per compound

| Column | Required? | Units | Description |
|---|---|---|---|
| `compound_id` | **Required** | — | Unique compound identifier |
| `IC50_nM` | **Required** | nM | IC50 value |
| `target` | Optional | — | Pharmacological target name |
| `MW` | Optional* | g/mol | Molecular weight |
| `fu_plasma_human` | Optional* | 0–1 | Fraction unbound, human plasma |
| `fu_plasma_rat` | Optional* | 0–1 | Fraction unbound, rat plasma |
| `fu_plasma_mouse` | Optional | 0–1 | Fraction unbound, mouse plasma |
| `solubility_ug_mL` | Optional* | µg/mL | Aqueous solubility |
| `formulation` | Optional | — | Vehicle/formulation used |
| `microsomal_t_half_min_human` | Optional* | min | Human microsomal stability t½ |
| `microsomal_t_half_min_rat` | Optional | min | Rat microsomal stability t½ |
| `papp_caco2_cm_s` | Optional | cm/s | Caco-2 permeability (PO route only) |
| `papp_pampa_cm_s` | Optional | cm/s | PAMPA permeability (PO route only) |
| `IC50_uM` | Optional | µM | Auto-calculated from IC50_nM if absent |
| `study_id` | Optional | — | Study identifier |
| `notes` | Optional | — | Free text |

*\*Needed for specific calculations: MW for dose in mg/kg, fu for dose projection, microsomal t½ for IVIVE, solubility for absorption assessment (PO only).*

### PK Data (CSV/Excel) — one row per compound × species × route

| Column | Required? | Units | Description |
|---|---|---|---|
| `compound_id` | **Required** | — | Must match in vitro file |
| `species` | **Required** | — | `rat`, `mouse`, `dog`, `monkey` |
| `route` | Optional | — | `IV`, `PO`, or `SC` |
| `dose_mg_kg` | Optional | mg/kg | Dose administered |
| `CL_mL_min_kg` | Optional | mL/min/kg | Clearance (typically from IV studies) |
| `Vss_L_kg` | Optional | L/kg | Volume of distribution (typically from IV) |
| `F_pct` | Optional | % | Bioavailability (PO or SC, relative to IV) |
| `t_half_h` | Optional | hours | Terminal half-life |
| `Cmax_ng_mL` | Optional | ng/mL | Peak concentration |
| `AUC_ng_h_mL` | Optional | ng·h/mL | Area under the curve |
| `study_id` | Optional | — | Study identifier |
| `notes` | Optional | — | Free text |

**Row structure:** A compound typically has **multiple rows** — one for each route tested. IV rows provide CL and Vss (leave F_pct blank). PO and SC rows provide F_pct, Cmax, AUC (leave CL/Vss blank). Example:

```
compound_id,species,route,dose_mg_kg,CL_mL_min_kg,Vss_L_kg,F_pct,t_half_h,AUC_ng_h_mL,study_id
CPD-001,rat,IV,2,25,1.5,,1.2,4800,PK-001
CPD-001,rat,PO,10,,,45,1.8,13500,PK-001
CPD-001,rat,SC,5,,,72,1.5,8640,PK-001
```

See `examples/` for complete sample data files.

## Jupyter Notebooks

Interactive examples in `notebooks/`:

1. **01_data_input_template** — Loading and validating data
2. **02_ivive_clearance** — IVIVE workflow from microsomal t½ to CLh
3. **03_dose_projection** — IC50-based dose projection for animal studies
4. **04_human_dose_estimation** — HED and MRSD from animal NOAEL
5. **05_full_workflow** — End-to-end dose projection pipeline with PO/IV/SC route comparison

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
- **SC bioavailability data** — SC F must be measured (SC AUC / IV AUC); cannot be predicted from in vitro data like oral F
- **Injection site absorption rate (SC)** — affects Cmax and Tmax for SC dosing; not captured by steady-state equations

## Running Tests

```bash
pip install pytest
python -m pytest tests/ -v
```

## License

MIT
