#!/usr/bin/env python3
"""
Dose Projection CLI — Run the full pipeline from the terminal.

Usage:
    python run_projection.py examples/example_invitro_data.csv examples/example_pk_data.csv
    python run_projection.py my_invitro.csv my_pk.csv --route iv
    python run_projection.py my_invitro.csv my_pk.csv --route sc --species rat --tau 12
    python run_projection.py my_invitro.csv my_pk.csv --output results.csv

Accepts CSV or Excel (.xlsx) files. Outputs a summary table to the terminal
and optionally saves results to a CSV file.
"""

import argparse
import sys
import warnings

import pandas as pd

from doseprojection.io import load_invitro_data, load_pk_data
from doseprojection.dose_projection import efficacious_dose_mg_kg, steady_state_css
from doseprojection.ivive import ivive_workflow
from doseprojection.absorption import classify_bcs, classify_permeability, predict_fa
from doseprojection.human_dose import hed_from_noael


def parse_args():
    parser = argparse.ArgumentParser(
        description="Preclinical dose projection from in vitro and PK data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with example data (oral, rat, QD):
  python run_projection.py examples/example_invitro_data.csv examples/example_pk_data.csv

  # IV dosing in rat:
  python run_projection.py invitro.csv pk.csv --route iv --species rat

  # SC dosing, BID (every 12h), save to file:
  python run_projection.py invitro.csv pk.csv --route sc --tau 12 --output results.csv

  # 3x IC50 coverage with custom NOAEL:
  python run_projection.py invitro.csv pk.csv --coverage 3 --noael 50
        """,
    )
    parser.add_argument("invitro_file", help="Path to in vitro data file (CSV or Excel)")
    parser.add_argument("pk_file", help="Path to PK data file (CSV or Excel)")
    parser.add_argument("--route", default="po", choices=["po", "iv", "sc"],
                        help="Route of administration (default: po)")
    parser.add_argument("--species", default="rat",
                        help="Species for dose projection (default: rat)")
    parser.add_argument("--tau", type=float, default=24,
                        help="Dosing interval in hours (default: 24 = once daily)")
    parser.add_argument("--coverage", type=float, default=1.0,
                        help="IC50 coverage multiple (default: 1.0)")
    parser.add_argument("--noael", type=float, default=None,
                        help="NOAEL (mg/kg) for HED/MRSD calculation (optional)")
    parser.add_argument("--output", "-o", default=None,
                        help="Save results to this CSV file (optional)")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Suppress warnings")
    return parser.parse_args()


def run(args):
    if args.quiet:
        warnings.filterwarnings("ignore")

    # --- Load data ---
    print(f"Loading in vitro data from: {args.invitro_file}")
    invitro = load_invitro_data(args.invitro_file)
    print(f"  Found {len(invitro)} compound(s)")

    print(f"Loading PK data from: {args.pk_file}")
    pk = load_pk_data(args.pk_file)
    print(f"  Found {len(pk)} PK record(s)")

    species = args.species.lower()
    route = args.route.lower()
    tau = args.tau
    coverage = args.coverage

    print(f"\nSettings: route={route.upper()}, species={species}, "
          f"tau={tau}h, coverage={coverage}x IC50")
    print("=" * 70)

    # --- Get PK data for the chosen species ---
    pk_iv = pk[(pk["species"] == species) & (pk["route"] == "IV")]
    pk_extravascular = None
    if route in ("po", "sc"):
        route_label = "PO" if route == "po" else "SC"
        pk_extravascular = pk[(pk["species"] == species) &
                              (pk["route"] == route_label)]

    # --- Project dose for each compound ---
    results = []
    for _, cpd in invitro.iterrows():
        cid = cpd["compound_id"]
        row = {"compound_id": cid}

        # Get IC50 and MW
        ic50 = cpd.get("IC50_nM")
        mw = cpd.get("MW")
        if pd.isna(ic50) or pd.isna(mw):
            row["error"] = "Missing IC50_nM or MW"
            results.append(row)
            continue

        row["IC50_nM"] = ic50
        row["MW"] = mw

        # Get fu for this species
        fu_col = f"fu_plasma_{species}"
        fu = cpd.get(fu_col)
        if pd.isna(fu):
            row["error"] = f"Missing {fu_col}"
            results.append(row)
            continue
        row["fu"] = fu

        # Get CL from IV data
        cpd_iv = pk_iv[pk_iv["compound_id"] == cid]
        if cpd_iv.empty:
            row["error"] = f"No IV PK data for {species}"
            results.append(row)
            continue
        cl = cpd_iv.iloc[0]["CL_mL_min_kg"]
        row["CL_mL_min_kg"] = cl

        # Get F for PO/SC routes
        f_val = None
        if route == "iv":
            f_val = 1.0
            row["F"] = 1.0
        elif pk_extravascular is not None:
            cpd_ev = pk_extravascular[pk_extravascular["compound_id"] == cid]
            if cpd_ev.empty:
                row["error"] = f"No {route.upper()} PK data for {species}"
                results.append(row)
                continue
            f_pct = cpd_ev.iloc[0].get("F_pct")
            if pd.isna(f_pct):
                row["error"] = f"Missing F_pct in {route.upper()} data"
                results.append(row)
                continue
            f_val = f_pct / 100.0
            row["F"] = f_val

        # --- Dose projection ---
        try:
            dose = efficacious_dose_mg_kg(
                ic50_nm=ic50, mw=mw, cl_mL_min_kg=cl,
                fu=fu, tau_h=tau, f=f_val,
                coverage_multiple=coverage, route=route
            )
            row["dose_mg_kg"] = round(dose, 2)

            css = steady_state_css(
                dose, cl, tau, f=f_val, route=route
            )
            row["Css_total_ng_mL"] = round(css, 1)
            row["Cu_ss_ng_mL"] = round(css * fu, 2)
        except Exception as e:
            row["error"] = str(e)
            results.append(row)
            continue

        # --- IVIVE (if human microsomal data available) ---
        t_half_human = cpd.get("microsomal_t_half_min_human")
        fu_human = cpd.get("fu_plasma_human")
        if not pd.isna(t_half_human) and not pd.isna(fu_human):
            try:
                ivive = ivive_workflow(t_half_human, fu_human, species="human")
                row["human_CLh_mL_min_kg"] = round(ivive["cl_hepatic_mL_min_kg"], 2)
                row["human_Fh"] = round(ivive["fh"], 3)
            except Exception:
                pass

        # --- Absorption assessment (PO only) ---
        if route == "po":
            papp = cpd.get("papp_caco2_cm_s")
            sol = cpd.get("solubility_ug_mL")
            if not pd.isna(papp):
                row["permeability"] = classify_permeability(papp)
                row["predicted_Fa"] = round(predict_fa(papp), 2)
            if not pd.isna(sol) and not pd.isna(papp):
                total_dose_mg = dose * 60  # approximate for 60 kg human
                bcs = classify_bcs(sol / 1000, papp, dose_mg=total_dose_mg)
                row["BCS_class"] = bcs["bcs_class"]

        # --- HED/MRSD (if NOAEL provided) ---
        if args.noael is not None:
            try:
                hed = hed_from_noael(args.noael, species)
                row["HED_mg_kg"] = round(hed["hed_mg_kg"], 2)
                row["MRSD_mg_kg"] = round(hed["mrsd_mg_kg"], 3)
                row["MRSD_total_mg"] = round(hed["mrsd_total_mg"], 1)
            except Exception:
                pass

        row["error"] = ""
        results.append(row)

    # --- Display results ---
    df = pd.DataFrame(results)

    # Reorder columns for readability
    priority_cols = [
        "compound_id", "IC50_nM", "MW", "fu", "CL_mL_min_kg", "F",
        "dose_mg_kg", "Css_total_ng_mL", "Cu_ss_ng_mL",
    ]
    ordered = [c for c in priority_cols if c in df.columns]
    remaining = [c for c in df.columns if c not in ordered]
    df = df[ordered + remaining]

    print(f"\n{'=' * 70}")
    print(f"DOSE PROJECTION RESULTS — {route.upper()} route, {species}, "
          f"Q{int(tau)}H, {coverage}x IC50")
    print(f"{'=' * 70}\n")

    # Print each compound
    for _, row in df.iterrows():
        cid = row["compound_id"]
        if row.get("error"):
            print(f"  {cid}: SKIPPED — {row['error']}")
            continue

        print(f"  {cid}:")
        print(f"    IC50 = {row['IC50_nM']} nM, MW = {row['MW']}, "
              f"fu = {row['fu']}, CL = {row['CL_mL_min_kg']} mL/min/kg")
        print(f"    F = {row['F']}, Route = {route.upper()}")
        print(f"    --> Projected dose: {row['dose_mg_kg']} mg/kg")
        print(f"    --> Css,total: {row['Css_total_ng_mL']} ng/mL, "
              f"Cu,ss: {row['Cu_ss_ng_mL']} ng/mL")

        if "human_CLh_mL_min_kg" in row and not pd.isna(row.get("human_CLh_mL_min_kg")):
            print(f"    --> Human CLh (IVIVE): {row['human_CLh_mL_min_kg']} mL/min/kg, "
                  f"Fh: {row['human_Fh']}")

        if route == "po":
            if "permeability" in row and not pd.isna(row.get("permeability")):
                print(f"    --> Permeability: {row['permeability']}, "
                      f"Predicted Fa: {row['predicted_Fa']}")
            if "BCS_class" in row and not pd.isna(row.get("BCS_class")):
                print(f"    --> BCS Class: {row['BCS_class']}")

        if "HED_mg_kg" in row and not pd.isna(row.get("HED_mg_kg")):
            print(f"    --> HED: {row['HED_mg_kg']} mg/kg, "
                  f"MRSD: {row['MRSD_mg_kg']} mg/kg "
                  f"({row['MRSD_total_mg']} mg total)")
        print()

    # --- Save to file ---
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Results saved to: {args.output}")

    # Summary
    successful = df[df["error"] == ""] if "error" in df.columns else df
    skipped = df[df["error"] != ""] if "error" in df.columns else pd.DataFrame()
    print(f"Done: {len(successful)} compound(s) projected, "
          f"{len(skipped)} skipped.")

    return df


def main():
    args = parse_args()
    try:
        run(args)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Data error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
