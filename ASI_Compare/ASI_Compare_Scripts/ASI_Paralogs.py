#!/usr/bin/env python3

import os
import csv
import subprocess
import argparse
from pathlib import Path
from collections import defaultdict


# SETTINGS


parser = argparse.ArgumentParser(description="Identify paralogs within each genome using Foldseek TM scores")

parser.add_argument("--input-dir", required=True,
                    help="Directory containing input PDB files")

parser.add_argument("--workdir", required=True,
                    help="Working/output directory")

parser.add_argument("--threads", type=int, default=8,
                    help="Number of CPU threads to use")

parser.add_argument("--sensitivity", type=str, default="9.5",
                    help="Foldseek sensitivity parameter")

parser.add_argument("--tm-threshold", type=float, default=0.80,
                    help="Minimum TM-score to consider proteins paralogs")

args = parser.parse_args()

INPUT_DIR = args.input_dir
WORKDIR = args.workdir
THREADS = args.threads
SENSITIVITY = args.sensitivity
TM_THRESHOLD = args.tm_threshold

THREAD_LIMIT = min(THREADS, os.cpu_count() or THREADS)


# THREAD CONTROL


os.environ["OMP_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["MKL_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["OPENBLAS_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["NUMEXPR_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(THREAD_LIMIT)
os.environ["MMSEQS_NUM_THREADS"] = str(THREAD_LIMIT)

print("THREAD LIMIT ENFORCED:", THREAD_LIMIT)


# PATHS


SPLIT_DIR = f"{WORKDIR}/split"
DB_PATH = f"{WORKDIR}/db"
RESULTS_PATH = f"{WORKDIR}/results"
TMP_PATH = f"{WORKDIR}/tmp"

RAW_TSV = f"{WORKDIR}/self_vs_self.tsv"
MAP_FILE = f"{WORKDIR}/protein_map.tsv"

PARALOG_TSV = f"{WORKDIR}/paralogs.tsv"


# CREATE DIRECTORIES


Path(WORKDIR).mkdir(parents=True, exist_ok=True)
Path(SPLIT_DIR).mkdir(parents=True, exist_ok=True)


# STORAGE


protein_to_genome = {}


# STEP 1: SPLIT PROTEOMES


print("Splitting proteomes...")

for pdb in Path(INPUT_DIR).glob("*.pdb"):

    genome = pdb.stem

    with open(pdb) as f:
        lines = f.readlines()

    current = []
    current_protein_id = None

    for line in lines:

        if line.startswith("REMARK") and "#" in line:
            try:
                gene_id = line.split()[1]
            except Exception:
                continue

            new_protein_id = gene_id

            if current and current_protein_id != new_protein_id:
                out = f"{SPLIT_DIR}/{current_protein_id}.pdb"
                with open(out, "w") as w:
                    w.writelines(current)
                protein_to_genome[current_protein_id] = genome
                current = []

            current_protein_id = new_protein_id

        if line.startswith(("ATOM", "TER", "END", "REMARK")):
            current.append(line)

    if current and current_protein_id:
        out = f"{SPLIT_DIR}/{current_protein_id}.pdb"
        with open(out, "w") as w:
            w.writelines(current)
        protein_to_genome[current_protein_id] = genome

print("Proteins:", len(protein_to_genome))

# Save map
with open(MAP_FILE, "w") as f:
    for p, g in protein_to_genome.items():
        f.write(f"{p}\t{g}\n")


# STEP 2: CREATE DB


print("\nCreating Foldseek database...")

subprocess.run([
    "foldseek", "createdb",
    SPLIT_DIR,
    DB_PATH
], check=True)


# STEP 3: SEARCH (self vs self)


print("Running Foldseek search (self vs self)...")

subprocess.run([
    "foldseek", "search",
    DB_PATH, DB_PATH,
    RESULTS_PATH, TMP_PATH,
    "--threads", str(THREAD_LIMIT),
    "-s", SENSITIVITY,
    "-a"
], check=True)


# STEP 4: CONVERT


print("Converting alignments...")

subprocess.run([
    "foldseek", "convertalis",
    DB_PATH, DB_PATH,
    RESULTS_PATH,
    RAW_TSV,
    "--threads", str(THREAD_LIMIT),
    "--format-output",
    "query,target,qtmscore,ttmscore,alntmscore,evalue"
], check=True)


# STEP 5: IDENTIFY PARALOGS


print("\nIdentifying paralogs (TM ≥", TM_THRESHOLD, ")")

used = set()
count_paralogs = 0

with open(PARALOG_TSV, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Genome", "Protein_A", "Protein_B", "TM_score", "AlnTM_score"])

    with open(RAW_TSV) as f:
        reader = csv.reader(f, delimiter="\t")

        for row in reader:
            try:
                q = row[0]
                t = row[1]
                tm = float(row[2])
                aln = float(row[4])
            except:
                continue

            if q == t:
                continue

            if q not in protein_to_genome or t not in protein_to_genome:
                continue

            gq = protein_to_genome[q]
            gt = protein_to_genome[t]

            # Only same-genome comparisons
            if gq != gt:
                continue

            if tm < TM_THRESHOLD:
                continue

            pair = tuple(sorted((q, t)))
            if pair in used:
                continue
            used.add(pair)

            writer.writerow([gq, pair[0], pair[1], round(tm, 4), round(aln, 4)])
            count_paralogs += 1

print("Paralog pairs found:", count_paralogs)
print("Output:", PARALOG_TSV)

# Cleanup
if os.path.exists(RAW_TSV):
    os.remove(RAW_TSV)

# CLEANUP: Remove all files except the paralog TSV


print("\nCleaning up intermediate files...")

for item in Path(WORKDIR).iterdir():
    # Keep only the paralog TSV
    if item.name == os.path.basename(PARALOG_TSV):
        continue
    # Remove directories recursively
    if item.is_dir():
        subprocess.run(["rm", "-rf", str(item)])
    else:
        item.unlink()

print("✓ Cleanup complete — only paralog file retained.")



print("\nDONE")
