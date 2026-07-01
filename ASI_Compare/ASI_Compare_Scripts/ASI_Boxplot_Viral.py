#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse


# ARGPARSE

parser = argparse.ArgumentParser(description="Part 7: Boxplot of Avg_TM and Mean_AAI by taxonomic rank")

parser.add_argument("--input", required=True,
                    help="TMfiltered_with_AAI_and_Taxonomy.csv from Part 6")

parser.add_argument("--outdir", required=True,
                    help="Output directory for figures")

args = parser.parse_args()

INPUT = args.input
OUTDIR = Path(args.outdir)
OUTDIR.mkdir(parents=True, exist_ok=True)


# LOAD DATA


df = pd.read_csv(INPUT)
df.columns = df.columns.str.strip()


# TAXONOMIC RANKS (UPDATED FOR YOUR DATASET)


tax_levels = [
    "Realm", "Subrealm",
    "Kingdom", "Subkingdom",
    "Phylum", "Subphylum",
    "Class", "Subclass",
    "Order", "Suborder",
    "Family", "Subfamily",
    "Genus", "Subgenus",
    "Species"
]


# DETERMINE CLOSEST TAXONOMIC RANK


def closest_rank(row):
    for lvl in reversed(tax_levels):  # Species → Realm
        colA = f"{lvl}_A"
        colB = f"{lvl}_B"
        if colA in row and colB in row:
            if pd.notna(row[colA]) and pd.notna(row[colB]):
                if row[colA] == row[colB]:
                    return lvl
    return "Interrealm"

df["Closest_Taxonomic_Rank"] = df.apply(closest_rank, axis=1)


# BOXPLOT: Avg_TM vs Mean_AAI by taxonomic rank


rank_order = ["Interrealm"] + tax_levels

plot_data_tm = []
plot_data_aai = []
labels = []

for rank in rank_order:
    subset = df[df["Closest_Taxonomic_Rank"] == rank]

    tm_vals = pd.to_numeric(subset["Avg_TM"], errors="coerce").dropna()
    aai_vals = pd.to_numeric(subset["Mean_AAI"], errors="coerce").dropna()

    if len(tm_vals) == 0 or len(aai_vals) == 0:
        continue

    plot_data_tm.append(tm_vals)
    plot_data_aai.append(aai_vals)
    labels.append(f"{rank}\n(n={len(subset)})")


# DRAW FIGURE


plt.figure(figsize=(16, 8))

x = np.arange(len(plot_data_tm)) * 2.0

# TM boxplots
plt.boxplot(
    plot_data_tm,
    positions=x,
    widths=0.6,
    patch_artist=True,
    boxprops=dict(facecolor="steelblue", alpha=0.6)
)

# AAI boxplots
plt.boxplot(
    plot_data_aai,
    positions=x + 0.7,
    widths=0.6,
    patch_artist=True,
    boxprops=dict(facecolor="salmon", alpha=0.6)
)

plt.xticks(x + 0.35, labels, rotation=45)
plt.ylabel("Similarity Score")
plt.title("Avg_TM vs Mean_AAI by Closest Taxonomic Rank")

plt.legend(
    handles=[
        plt.Rectangle((0,0),1,1, color="steelblue", alpha=0.6, label="Avg_TM"),
        plt.Rectangle((0,0),1,1, color="salmon", alpha=0.6, label="Mean_AAI")
    ],
    loc="upper right"
)

plt.tight_layout()

outpath = OUTDIR / "TM_vs_AAI_by_Taxonomic_Rank.png"
plt.savefig(outpath, dpi=300)
plt.close()

print(f"Saved figure: {outpath}")
