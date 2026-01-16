import pandas as pd
from ete3 import NCBITaxa

# Initialize NCBI
ncbi = NCBITaxa()

known_domains = {"Bacteria", "Archaea", "Eukaryota"}

INPUT_CSV = "ASI_Strains_RefSeq.csv"             # your input file
OUTPUT_CSV = "Automated_Bacteria_Taxonomy.csv"  # output file

# Load CSV
df = pd.read_csv(INPUT_CSV)

# Extract genus from species name
df['Genus'] = df['filename'].apply(lambda x: x.split()[0])

# Get unique genera
unique_genera = df['Genus'].unique()

# Cache taxonomy results
taxonomy_cache = {}

for genus in unique_genera:
    try:
        taxid = ncbi.get_name_translator([genus])[genus][0]
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        taxonomy = {}

        # Find domain
        domain_taxid = None
        for tid in lineage:
            if names[tid] in known_domains:
                domain_taxid = tid
                break
        taxonomy["Domain"] = names[domain_taxid] if domain_taxid else ""

        # Other ranks
        for tid in lineage:
            rank = ranks.get(tid, "")
            name = names[tid]
            if rank == "kingdom":
                taxonomy["Kingdom"] = name
            elif rank == "phylum":
                taxonomy["Phylum"] = name
            elif rank == "class":
                taxonomy["Class"] = name
            elif rank == "order":
                taxonomy["Order"] = name
            elif rank == "family":
                taxonomy["Family"] = name
            elif rank == "genus":
                taxonomy["Genus"] = name

        # Fill missing
        full_tax = {k: taxonomy.get(k, "") for k in ["Domain","Kingdom","Phylum","Class","Order","Family","Genus"]}

        taxonomy_cache[genus] = full_tax

    except Exception as e:
        # If genus not found, fill blanks
        taxonomy_cache[genus] = {k: "" for k in ["Domain","Kingdom","Phylum","Class","Order","Family","Genus"]}

# Convert cache to DataFrame
tax_df = pd.DataFrame(taxonomy_cache.values())

# Remove duplicates (just in case)
tax_df = tax_df.drop_duplicates()

# Save CSV
tax_df.to_csv(OUTPUT_CSV, index=False)

print(f"Saved unique taxonomy CSV to {OUTPUT_CSV}")
