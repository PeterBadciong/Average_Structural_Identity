import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("ESM2_with_taxonomy_COMP_AVG.csv")

# Set the order of taxonomy levels
levels = ["domain","kingdom","class","order","family","genus","species"]

plt.figure(figsize=(15, 9))
ax = sns.boxplot(
    x="closest_taxonomy", 
    y="combined_avg", 
    data=df, 
    order=levels,
    palette="Greys"  # greyscale
)

plt.ylabel("((avg_rbh_score×rbh_fraction_symmetric)+symmetric_avg+greedy_1to1_avg)/3")
plt.xlabel("Closest Taxonomy Level")
plt.title("Similarity by Taxonomy Level")

# Compute global range for positioning annotations
global_range = df["combined_avg"].max() - df["combined_avg"].min()

for i, level in enumerate(levels):
    vals = df.loc[df["closest_taxonomy"] == level, "combined_avg"].dropna()
    n_total = len(vals)
    if n_total == 0:
        continue

    # Annotate number of comparisons below the box
    y_pos = df["combined_avg"].min() - 0.05 * global_range
    ax.text(
        i,
        y_pos,
        f"n={n_total}",
        ha="center",
        va="top",  # place text above the y position
        fontsize=10,
        fontweight="bold",
        color="black"
    )

# Adjust plot limits to make sure the annotation is visible
ax.set_ylim(df["combined_avg"].min() - 0.1 * global_range, df["combined_avg"].max() + 0.05 * global_range)

# Save plot
plt.savefig("similarity_by_taxonomy_combined_avg_COMP_greyscale.png", dpi=300, bbox_inches="tight")
plt.close()
print("Boxplot saved as similarity_by_taxonomy_combined_avg_COMP_greyscale.png")
