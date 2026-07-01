#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse


# ARGPARSE


parser = argparse.ArgumentParser(description="Part 6: Merge TMfiltered RBH + AAI + Taxonomy")

parser.add_argument("--tm-csv", required=True,
                    help="proteome_RBH_scored_TMfiltered.csv from Foldseek Part 1")

parser.add_argument("--aai-summary", required=True,
                    help="pairwise_summary.csv from Part 4")

parser.add_argument("--taxonomy", required=True,
                    help="Taxonomy key CSV with columns File,Domain,Kingdom,Phylum,Class,Order,Family,Genus")

parser.add_argument("--outdir", required=True,
                    help="Output directory")

args = parser.parse_args()

TM_CSV = args.tm_csv
AAI_SUMMARY = args.aai_summary
TAX = args.taxonomy
OUTDIR = Path(args.outdir)
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "TMfiltered_with_AAI_and_Taxonomy.csv"


# LOAD FILES


tm_df = pd.read_csv(TM_CSV)
aai_df = pd.read_csv(AAI_SUMMARY)
tax_df = pd.read_csv(TAX)

tm_df.columns = tm_df.columns.str.strip()
aai_df.columns = aai_df.columns.str.strip()
tax_df.columns = tax_df.columns.str.strip()


# NORMALIZE GENOME PAIR ORDERING


def normalize_pair(a, b):
    a = str(a).strip()
    b = str(b).strip()
    return tuple(sorted([a, b]))

tm_df["NormA"], tm_df["NormB"] = zip(*tm_df.apply(
    lambda r: normalize_pair(r["Genome_A"], r["Genome_B"]),
    axis=1
))

aai_df["NormA"], aai_df["NormB"] = zip(*aai_df.apply(
    lambda r: normalize_pair(r["GenomeA"], r["GenomeB"]),
    axis=1
))


# MERGE TM RBHs WITH AAI SUMMARY


merged = tm_df.merge(
    aai_df[["NormA", "NormB", "Mean_AAI", "Total_RBH"]],
    on=["NormA", "NormB"],
    how="left"
)


# RENAME COLUMNS


rename_map = {
    "RBH_count": "ASI_RBHs",
    "Total_RBH": "AAI_RBHs"
}

for old, new in rename_map.items():
    if old in merged.columns:
        merged = merged.rename(columns={old: new})


# CLEAN UP TEMP COLUMNS


merged = merged.drop(columns=["NormA", "NormB"])


# NORMALIZE TAXONOMY KEY


tax_df = tax_df.rename(columns={"File": "Genome"})
tax_df["Genome"] = tax_df["Genome"].astype(str).str.strip()


# MERGE TAXONOMY FOR GENOME_A


taxA = tax_df.copy()
taxA.columns = ["Genome_A"] + [f"{c}_A" for c in tax_df.columns if c != "Genome"]

merged = merged.merge(taxA, on="Genome_A", how="left")


# MERGE TAXONOMY FOR GENOME_B


taxB = tax_df.copy()
taxB.columns = ["Genome_B"] + [f"{c}_B" for c in tax_df.columns if c != "Genome"]

merged = merged.merge(taxB, on="Genome_B", how="left")


# SAVE OUTPUT


merged.to_csv(OUTPUT, index=False)

print(f"Saved final ASI+AAI+Taxonomy table: {OUTPUT}")
