import pandas as pd
import numpy as np
import os

# =====================
# User settings
# =====================
INPUT_CSV = "pairwise_RBH_similarity.csv"         # original pairwise similarity
TAXONOMY_CSV = "Automated_Bacteria_Taxonomy.csv"  # taxonomy info
OUTPUT_CSV = "ESM2_with_taxonomy_COMP_TEMP.csv"         # enriched pairwise CSV

# Get base name without extension for dynamic metric files
output_base = os.path.splitext(OUTPUT_CSV)[0]

# Metrics → dynamic output stats files
METRICS = {
    "greedy_1to1_avg": f"{output_base}_taxonomy_stats_greedy_1to1_avg.csv",
    "symmetric_avg": f"{output_base}_taxonomy_stats_symmetric_avg.csv",
    "symmetric_coverage": f"{output_base}_taxonomy_stats_symmetric_coverage.csv",
    "avg_rbh_score": f"{output_base}_taxonomy_stats_avg_rbh_score.csv",
    "rbh_fraction_symmetric": f"{output_base}_taxonomy_stats_rbh_fraction_symmetric.csv"
}

# Example: printing the file names
for metric, file in METRICS.items():
    print(metric, "→", file)



# =====================
# Exclusion rules
# =====================
# Format: {taxonomy_level: [values_to_exclude]}
# Case-insensitive. Either A or B matching is excluded.
EXCLUDE_TAXA = {
    "domain": ["Arch"],   # exclude Domain = Archaea
    # "phylum": ["firmicutes"],  # example for later
}

# =====================
# Load pairwise similarity CSV
# =====================
df = pd.read_csv(INPUT_CSV)

# =====================
# Helper function to extract genus/species
# =====================
def split_genome_name(name):
    parts = str(name).split("_", 2)
    genus = parts[0]
    species = parts[1] if len(parts) > 1 else ""
    return genus, species

# =====================
# Extract genus/species for Genome A and Genome B
# =====================
A_genus, A_species = zip(*df["Genome_A"].map(split_genome_name))
B_genus, B_species = zip(*df["Genome_B"].map(split_genome_name))

# =====================
# Add genus/species columns to dataframe
# =====================
df = pd.concat(
    [
        df.iloc[:, 0:2],  # Genome_A, Genome_B
        pd.DataFrame({"a_genus": A_genus, "a_species": A_species}),
        pd.DataFrame({"b_genus": B_genus, "b_species": B_species}),
        df.iloc[:, 2:]    # rest of original columns
    ],
    axis=1
)

# =====================
# Load taxonomy CSV
# =====================
taxonomy = pd.read_csv(TAXONOMY_CSV)

# Normalize column names
taxonomy.columns = taxonomy.columns.str.strip().str.lower()

# Expected taxonomy columns
expected_cols = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]

# Check which columns actually exist
existing_cols = [c for c in expected_cols if c in taxonomy.columns]
missing_cols = [c for c in expected_cols if c not in taxonomy.columns]

if missing_cols:
    print(f"⚠ Warning: these expected taxonomy columns are missing in the CSV: {missing_cols}")

taxonomy = taxonomy[existing_cols]

# Normalize df column names
df.columns = df.columns.str.strip().str.lower()

# =====================
# Merge taxonomy for Genome A and B
# =====================
df = df.merge(taxonomy.add_prefix("a_"), left_on="a_genus", right_on="a_genus", how="left")
df = df.merge(taxonomy.add_prefix("b_"), left_on="b_genus", right_on="b_genus", how="left")

# Clean up duplicate join keys if any
df = df.drop(columns=["a_genus_y", "b_genus_y"], errors="ignore")
df = df.rename(columns={"a_genus_x": "a_genus", "b_genus_x": "b_genus"})

# =====================
# Reorder columns taxonomically
# =====================
A_tax_cols = ["a_domain","a_kingdom","a_phylum","a_class","a_order","a_family","a_genus","a_species"]
B_tax_cols = ["b_domain","b_kingdom","b_phylum","b_class","b_order","b_family","b_genus","b_species"]
core_cols = ["genome_a","genome_b"]
metric_cols = [c for c in df.columns if c not in core_cols + A_tax_cols + B_tax_cols]

df = df[core_cols + A_tax_cols + B_tax_cols + metric_cols]

# =====================
# Determine closest taxonomy match
# =====================
taxonomy_ranks = [
    ("domain","a_domain","b_domain"),
    ("kingdom","a_kingdom","b_kingdom"),
    ("phylum","a_phylum","b_phylum"),
    ("class","a_class","b_class"),
    ("order","a_order","b_order"),
    ("family","a_family","b_family"),
    ("genus","a_genus","b_genus"),
    ("species","a_species","b_species")
]

def closest_taxonomy_match(row):
    for level, a_col, b_col in reversed(taxonomy_ranks):
        if pd.notna(row[a_col]) and row[a_col] == row[b_col]:
            return level
    return "none"

df["closest_taxonomy"] = df.apply(closest_taxonomy_match, axis=1)

# =====================
# Apply taxonomy exclusions
# =====================
def apply_taxonomy_exclusions(df, rules):
    mask = pd.Series(True, index=df.index)

    for level, values in rules.items():
        a_col = f"a_{level}"
        b_col = f"b_{level}"

        if a_col not in df.columns or b_col not in df.columns:
            print(f"⚠ Warning: exclusion level '{level}' not found, skipping.")
            continue

        values = [v.lower() for v in values]

        level_mask = ~(
            df[a_col].astype(str).str.lower().isin(values) |
            df[b_col].astype(str).str.lower().isin(values)
        )

        mask &= level_mask

    removed = len(df) - mask.sum()
    if removed > 0:
        print(f"✓ Excluded {removed} rows based on taxonomy rules.")

    return df[mask]

df = apply_taxonomy_exclusions(df, EXCLUDE_TAXA)

# =====================
# Compute stats for each metric
# =====================
for metric, out_csv in METRICS.items():
    if metric not in df.columns:
        print(f"⚠ Warning: {metric} not found in dataframe, skipping.")
        continue

    output_rows = []

    for level, _, _ in taxonomy_ranks:
        subset = df[df["closest_taxonomy"] == level]
        if subset.empty:
            continue

        values = subset[metric].dropna()
        if values.empty:
            continue

        stats = {
            "taxonomy_level": level,
            "count": len(values),
            "mean": values.mean(),
            "std": values.std(),
            "min": values.min(),
            "max": values.max(),
            "range": values.max() - values.min()
        }
        output_rows.append(stats)

    stats_df = pd.DataFrame(output_rows)
    stats_df.to_csv(out_csv, index=False)
    print(f"✓ {metric} taxonomy stats saved to: {out_csv}")

# =====================
# Save enriched output
# =====================
df.to_csv(OUTPUT_CSV, index=False)
print(f"✓ Taxonomy successfully added. Saved to: {OUTPUT_CSV}")
