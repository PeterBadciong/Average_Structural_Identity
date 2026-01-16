import pandas as pd
import os
import csv
import re
import requests
import zipfile
import io
from tqdm import tqdm

# =========================
# 🔧 User settings
# =========================
INPUT_TSV = "assemblies.tsv"  # your input TSV
OUTPUT_CSV = "ASI_Strains_RefSeq_Staphylococcus.csv"
OUTPUT_DIR = "ASI_Fasta_Files_Staphylococcus2"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================
# 🧼 Helper functions
# =========================
def format_name(row):
    """Generate full organism name and pick best RefSeq accession"""
    organism_name = str(row['Organism Name'])
    parts = organism_name.split()
    genus = parts[0] if len(parts) > 0 else ''
    species = parts[1] if len(parts) > 1 else ''
    
    strain = row.get('Organism Infraspecific Names Strain', '')
    full_name = f"{genus} {species}"
    if pd.notna(strain) and strain != '':
        full_name += f" {strain}"
    
    # Prefer GCF over GCA if paired
    refseq = row['Assembly Accession']
    paired = str(row.get('Assembly Paired Assembly Accession', ''))
    if str(refseq).startswith("GCA_") and paired.startswith("GCF_"):
        refseq = paired
    
    return pd.Series([full_name, refseq])

def sanitize_filename(name):
    """Make filename safe but keep parentheses, spaces, dots, dashes, and underscores"""
    return re.sub(r'[^\w\-\.\(\) ]', '_', name)

def download_assembly_fasta(assembly, out_name):
    """Download genome FASTA from NCBI Datasets"""
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{assembly}/download?include_annotation_type=GENOME_FASTA"
    )
    try:
        r = requests.get(url)
        r.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            for file in z.namelist():
                if file.endswith(".fna"):
                    out_path = os.path.join(OUTPUT_DIR, out_name)
                    with z.open(file) as f, open(out_path, "wb") as out:
                        out.write(f.read())
                    return True
    except requests.HTTPError as e:
        print(f"⚠️ HTTP error for {assembly}: {e}")
    except Exception as e:
        print(f"⚠️ Error downloading {assembly}: {e}")
    return False

# =========================
# Step 1: Process TSV -> CSV
# =========================
df = pd.read_csv(INPUT_TSV, sep="\t")  # adjust sep="\t" if tab-delimited
output_df = df.apply(format_name, axis=1)
output_df.columns = ["Name", "RefSeq"]
output_df = output_df.drop_duplicates()

# Append to existing CSV (or create new)
if os.path.exists(OUTPUT_CSV):
    existing_df = pd.read_csv(OUTPUT_CSV, header=None, names=["Name", "RefSeq"])
    combined_df = pd.concat([existing_df, output_df], ignore_index=True)
    combined_df = combined_df.drop_duplicates()
    combined_df.to_csv(OUTPUT_CSV, index=False, header=False)
else:
    output_df.to_csv(OUTPUT_CSV, index=False, header=False)

print(f"✅ Added {len(output_df)} entries to '{OUTPUT_CSV}'.")

# =========================
# Step 2: Download FASTA files
# =========================
with open(OUTPUT_CSV, newline="") as f:
    reader = csv.reader(f)
    for row in tqdm(list(reader), desc="Downloading genomes"):
        if len(row) < 2:
            continue
        name, accession = row[0].strip(), row[1].strip()
        if not accession:
            print(f"⚠️ Skipping {name}, missing accession")
            continue

        safe_name = sanitize_filename(name)
        out_fasta = f"{safe_name}.fasta"
        out_path = os.path.join(OUTPUT_DIR, out_fasta)

        if os.path.exists(out_path):
            continue  # already downloaded

        success = download_assembly_fasta(accession, out_fasta)
        if not success:
            print(f"⚠️ Failed to find FASTA for {accession}")
