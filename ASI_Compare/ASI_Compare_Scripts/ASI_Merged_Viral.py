#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

# ARGPARSE

parser = argparse.ArgumentParser(
    description="Merge TM+AAI with ICTV taxonomy using species/genus lookup and deepest shared rank."
)

parser.add_argument("--tm-csv", required=True,
                    help="proteome_RBH_scored_TMfiltered.csv")

parser.add_argument("--aai-summary", required=True,
                    help="pairwise_summary.csv")

parser.add_argument("--taxonomy", required=True,
                    help="ICTV taxonomy CSV")

parser.add_argument("--outdir", required=True,
                    help="Output directory")

args = parser.parse_args()

TM_CSV = args.tm_csv
AAI_SUMMARY = args.aai_summary
TAX = args.taxonomy
OUTDIR = Path(args.outdir)
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUTDIR / "TM_AAI_Taxonomy_Merged.csv"

# LOAD INPUTS

tm_df = pd.read_csv(TM_CSV)
aai_df = pd.read_csv(AAI_SUMMARY)
tax_df = pd.read_csv(TAX)

tm_df.columns = tm_df.columns.str.strip()
aai_df.columns = aai_df.columns.str.strip()
tax_df.columns = tax_df.columns.str.strip()

# CLEAN + NORMALIZE TAXONOMY

# Normalize species names
if "Species" in tax_df.columns:
    tax_df["Species"] = tax_df["Species"].astype(str).str.replace(" ", "_")

# Build species/genus dictionaries
species_dict = {}
genus_dict = {}

for _, row in tax_df.iterrows():
    species = row["Species"]
    genus = row["Genus"]

    species_dict[species] = row.to_dict()

    if genus not in genus_dict:
        genus_dict[genus] = row.to_dict()

# Rank order (deepest → highest)
rank_order = [
    "Realm","Subrealm","Kingdom","Subkingdom","Phylum","Subphylum",
    "Class","Subclass","Order","Suborder","Family","Subfamily",
    "Genus","Subgenus","Species"
]

# TAXONOMY LOOKUP FUNCTIONS

def lookup_taxonomy(name):
    """Try species match first, then genus match."""
    if name in species_dict:
        return species_dict[name]

    genus = name.split("_")[0]
    if genus in genus_dict:
        return genus_dict[genus]

    return {}

def deepest_shared_rank(taxA, taxB):
    """Return deepest shared taxon and rank."""
    for rank in reversed(rank_order):  # deepest → highest
        if taxA.get(rank) and taxB.get(rank) and taxA[rank] == taxB[rank]:
            return taxA[rank], rank
    return "None", "None"

# NORMALIZE GENOME PAIRS

def normalize_pair(a, b):
    return tuple(sorted([str(a).strip(), str(b).strip()]))

tm_df["NormA"], tm_df["NormB"] = zip(*tm_df.apply(
    lambda r: normalize_pair(r["Genome_A"], r["Genome_B"]),
    axis=1
))

aai_df["NormA"], aai_df["NormB"] = zip(*aai_df.apply(
    lambda r: normalize_pair(r["GenomeA"], r["GenomeB"]),
    axis=1
))

# MERGE TM + AAI

merged = tm_df.merge(
    aai_df[["NormA", "NormB", "Mean_AAI", "Total_RBH"]],
    on=["NormA", "NormB"],
    how="left"
)

merged = merged.rename(columns={
    "RBH_count": "ASI_RBHs",
    "Total_RBH": "AAI_RBHs"
})

merged = merged.drop(columns=["NormA", "NormB"])

# MERGE TAXONOMY COLUMNS (KEEPING YOUR ORIGINAL FORMAT)

# Build A-side taxonomy columns
taxA = tax_df.copy().rename(columns={"Species": "Species_A"})
taxA = taxA.rename(columns={c: f"{c}_A" for c in tax_df.columns if c != "Species"})
taxA = taxA.rename(columns={"Species_A": "Genome_A"})

# Build B-side taxonomy columns
taxB = tax_df.copy().rename(columns={"Species": "Species_B"})
taxB = taxB.rename(columns={c: f"{c}_B" for c in tax_df.columns if c != "Species"})
taxB = taxB.rename(columns={"Species_B": "Genome_B"})

merged = merged.merge(taxA, on="Genome_A", how="left")
merged = merged.merge(taxB, on="Genome_B", how="left")

# APPLY NEW TAXONOMY LOGIC

closest_taxon = []
closest_rank = []

for _, row in merged.iterrows():
    gA = row["Genome_A"]
    gB = row["Genome_B"]

    taxA = lookup_taxonomy(gA)
    taxB = lookup_taxonomy(gB)

    taxon, rank = deepest_shared_rank(taxA, taxB)

    closest_taxon.append(taxon)
    closest_rank.append(rank)

merged["Closest_Taxon_Match"] = closest_taxon
merged["Closest_Rank"] = closest_rank

# SAVE OUTPUT

merged.to_csv(OUTPUT, index=False)
print(f"Saved merged TM+AAI+taxonomy file: {OUTPUT}")
