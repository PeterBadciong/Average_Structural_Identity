#!/usr/bin/env python3

import os
import csv
import subprocess
import argparse
from pathlib import Path
from collections import defaultdict


# ============================================================
# SETTINGS
# ============================================================

parser = argparse.ArgumentParser(description="Foldseek RBH + Bidirectional TM Pipeline")

parser.add_argument("--input-dir", required=True,
                    help="Directory containing input PDB files")

parser.add_argument("--workdir", required=True,
                    help="Working/output directory")

parser.add_argument("--threads", type=int, default=8,
                    help="Number of CPU threads to use")

parser.add_argument("--sensitivity", type=str, default="9.5",
                    help="Foldseek sensitivity parameter")

parser.add_argument("--tm-threshold", type=float, default=0.50,
                    help="Minimum TM-score for TM-filtered RBHs and bidirectional TM pairs")

args = parser.parse_args()

INPUT_DIR = args.input_dir
WORKDIR = args.workdir
THREADS = args.threads
THREAD_LIMIT = min(THREADS, os.cpu_count() or THREADS)
SENSITIVITY = args.sensitivity
TM_THRESHOLD = args.tm_threshold

# ============================================================
# THREAD CONTROL
# ============================================================

os.environ["OMP_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["MKL_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["OPENBLAS_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["NUMEXPR_NUM_THREADS"] = str(THREAD_LIMIT)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(THREAD_LIMIT)
os.environ["MMSEQS_NUM_THREADS"] = str(THREAD_LIMIT)

print("THREAD LIMIT ENFORCED:", THREAD_LIMIT)

# ============================================================
# PATHS
# ============================================================

SPLIT_DIR = f"{WORKDIR}/split"
DB_PATH = f"{WORKDIR}/db"
RESULTS_PATH = f"{WORKDIR}/results"
TMP_PATH = f"{WORKDIR}/tmp"

RAW_TSV = f"{WORKDIR}/all_vs_all.tsv"

MAP_FILE = f"{WORKDIR}/protein_map.tsv"

PROTEOME_OUT_ALL = f"{WORKDIR}/proteome_RBH_scored_ALL.csv"
PROTEOME_OUT_FILTERED = f"{WORKDIR}/proteome_RBH_scored_TMfiltered.csv"

PAIRWISE_MATCHES_ALL = f"{WORKDIR}/pairwise_RBH_matches_ALL.tsv"
PAIRWISE_MATCHES_TM = f"{WORKDIR}/pairwise_RBH_matches_TMfiltered.tsv"

# NEW: bidirectional TM-passing pairs (no best-hit requirement)
PAIRWISE_TM_BIDIRECTIONAL = f"{WORKDIR}/pairwise_TM_bidirectional.tsv"

# ============================================================
# CREATE DIRECTORIES
# ============================================================

Path(WORKDIR).mkdir(parents=True, exist_ok=True)
Path(SPLIT_DIR).mkdir(parents=True, exist_ok=True)

# ============================================================
# STORAGE
# ============================================================

protein_to_genome = {}
genome_size = defaultdict(int)

# ============================================================
# STEP 1: SPLIT PROTEOMES
# ============================================================

print("Splitting proteomes...")

protein_count = 0

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
                genome_size[genome] += 1
                protein_count += 1

                current = []

            current_protein_id = new_protein_id

        if line.startswith(("ATOM", "TER", "END", "REMARK")):
            current.append(line)

    if current and current_protein_id:

        out = f"{SPLIT_DIR}/{current_protein_id}.pdb"

        with open(out, "w") as w:
            w.writelines(current)

        protein_to_genome[current_protein_id] = genome
        genome_size[genome] += 1
        protein_count += 1

print("\nGenomes:", len(genome_size))
print("Proteins:", protein_count)

# ============================================================
# SAVE MAP
# ============================================================

print("Saving protein map...")

with open(MAP_FILE, "w") as f:
    for protein, genome in protein_to_genome.items():
        f.write(f"{protein}\t{genome}\n")

# ============================================================
# STEP 2: CREATE DB
# ============================================================

print("\nCreating Foldseek database...")

subprocess.run([
    "foldseek",
    "createdb",
    SPLIT_DIR,
    DB_PATH
], check=True)

# ============================================================
# STEP 3: SEARCH
# ============================================================

print("Running Foldseek search...")

subprocess.run([
    "foldseek",
    "search",
    DB_PATH,
    DB_PATH,
    RESULTS_PATH,
    TMP_PATH,
    "--threads", str(THREAD_LIMIT),
    "-s", SENSITIVITY,
    "-a"
], check=True)

# ============================================================
# STEP 4: CONVERT
# ============================================================

print("Converting alignments...")

subprocess.run([
    "foldseek",
    "convertalis",
    DB_PATH,
    DB_PATH,
    RESULTS_PATH,
    RAW_TSV,
    "--threads", str(THREAD_LIMIT),
    "--format-output",
    "query,target,qtmscore,ttmscore,alntmscore,evalue"
], check=True)

# ============================================================
# STEP 5: BUILD BEST HITS + TM-PASSING MAP
# ============================================================

print("Building best-hit dictionaries and TM-passing map...")

best_hit_all = {}
best_hit_tm = {}

# For bidirectional presence: store all TM-passing directional hits
# key: (query, target), value: (qtmscore, alntmscore)
tm_pass_hits = {}

rows = 0

with open(RAW_TSV) as f:

    reader = csv.reader(f, delimiter="\t")

    for row in reader:

        try:
            q = row[0]
            t = row[1]
            qt = float(row[2])
            aln = float(row[4])
        except Exception:
            continue

        if q == t:
            continue

        if q not in protein_to_genome or t not in protein_to_genome:
            continue

        gq = protein_to_genome[q]
        gt = protein_to_genome[t]

        if gq == gt:
            continue

        key = (q, gt)

        # -----------------------------
        # ALL best hit (no TM filter)
        # -----------------------------
        prev = best_hit_all.get(key)

        if prev is None or qt > prev[1]:
            best_hit_all[key] = (t, qt, aln)

        # -----------------------------
        # TM-filtered best hit
        # -----------------------------
        if qt >= TM_THRESHOLD:

            prev_tm = best_hit_tm.get(key)

            if prev_tm is None or qt > prev_tm[1]:
                best_hit_tm[key] = (t, qt, aln)

            # Store TM-passing directional hit for bidirectional presence
            tm_pass_hits[(q, t)] = (qt, aln)

        rows += 1

        if rows % 10_000_000 == 0:
            print(f"Processed {rows:,} Foldseek rows")

print("\nFinished loading Foldseek hits")
print("ALL best hits:", len(best_hit_all))
print("TM best hits:", len(best_hit_tm))
print("TM-passing directional hits (for bidirectional):", len(tm_pass_hits))

# ============================================================
# RBH FUNCTION (GENOME+PROTEIN NORMALIZATION)
# ============================================================

def compute_rbh(best_hit, pairwise_tsv, summary_csv, label):

    print(f"\nComputing RBHs → {os.path.basename(summary_csv)} [{label}]")

    used = set()

    rbh = defaultdict(lambda: {
        "tm": [],
        "aln": [],
        "count": 0
    })

    processed = 0

    with open(pairwise_tsv, "w", newline="") as out_f:

        writer = csv.writer(out_f, delimiter="\t")

        writer.writerow([
            "Genome_A",
            "Genome_B",
            "Protein_A_full_id",
            "Protein_B_full_id",
            "TM_score",
            "AlnTM_score"
        ])

        for (q, gt), (t, qt, aln) in best_hit.items():

            gq = protein_to_genome[q]

            reverse_key = (t, gq)

            reverse = best_hit.get(reverse_key)

            if reverse is None:
                continue

            reverse_target, reverse_tm, reverse_aln = reverse

            if reverse_target != q:
                continue

            pair_id = tuple(sorted((q, t)))

            if pair_id in used:
                continue

            used.add(pair_id)

            # Normalize genomes AND proteins together
            if gq <= gt:
                genomeA, genomeB = gq, gt
                protA, protB = q, t
            else:
                genomeA, genomeB = gt, gq
                protA, protB = t, q

            # WRITE PAIRWISE OUTPUT
            writer.writerow([
                genomeA,
                genomeB,
                protA,
                protB,
                round(qt, 4),
                round(aln, 4)
            ])

            # SUMMARY
            rbh[(genomeA, genomeB)]["tm"].extend([
                qt,
                reverse_tm
            ])

            rbh[(genomeA, genomeB)]["aln"].extend([
                aln,
                reverse_aln
            ])

            rbh[(genomeA, genomeB)]["count"] += 1

            processed += 1

            if processed % 1_000_000 == 0:
                print(f"RBHs processed [{label}]: {processed:,}")

    print(f"Building genome summary [{label}]...")

    final = []

    for (a, b), vals in rbh.items():

        avg_tm = sum(vals["tm"]) / len(vals["tm"])
        avg_aln = sum(vals["aln"]) / len(vals["aln"])

        count = vals["count"]

        denom = min(genome_size[a], genome_size[b])
        norm = count / denom if denom else 0

        score = avg_tm * norm

        final.append([
            a,
            b,
            avg_tm,
            avg_aln,
            count,
            norm,
            score
        ])

    final.sort(key=lambda x: x[6], reverse=True)

    with open(summary_csv, "w", newline="") as f:

        writer = csv.writer(f)

        writer.writerow([
            "Genome_A",
            "Genome_B",
            "Avg_TM",
            "Avg_alnTM",
            "RBH_count",
            "RBH_norm_min_size",
            "Composite_score"
        ])

        for r in final:

            writer.writerow([
                r[0],
                r[1],
                round(r[2], 4),
                round(r[3], 4),
                r[4],
                round(r[5], 6),
                round(r[6], 6)
            ])

    print("Finished RBH summary:", summary_csv)


# ============================================================
# BIDIRECTIONAL TM-PASSING FUNCTION
# ============================================================

def compute_bidirectional_tm_pairs(tm_pass_hits, out_tsv):

    print(f"\nComputing bidirectional TM-passing pairs → {os.path.basename(out_tsv)}")

    used = set()
    count_pairs = 0

    with open(out_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "Genome_A",
            "Genome_B",
            "Protein_A_full_id",
            "Protein_B_full_id",
            "TM_A_to_B",
            "TM_B_to_A",
            "AlnTM_A_to_B",
            "AlnTM_B_to_A"
        ])

        for (q, t), (qt, aln_qt) in tm_pass_hits.items():

            # Check reverse direction
            reverse = tm_pass_hits.get((t, q))
            if reverse is None:
                continue

            rt, aln_rt = reverse

            # Avoid double-counting
            pair_id = tuple(sorted((q, t)))
            if pair_id in used:
                continue
            used.add(pair_id)

            gq = protein_to_genome[q]
            gt = protein_to_genome[t]

            # Normalize genomes and proteins consistently
            if gq <= gt:
                genomeA, genomeB = gq, gt
                protA, protB = q, t
                tm_A_to_B, tm_B_to_A = qt, rt
                aln_A_to_B, aln_B_to_A = aln_qt, aln_rt
            else:
                genomeA, genomeB = gt, gq
                protA, protB = t, q
                tm_A_to_B, tm_B_to_A = rt, qt
                aln_A_to_B, aln_B_to_A = aln_rt, aln_qt

            writer.writerow([
                genomeA,
                genomeB,
                protA,
                protB,
                round(tm_A_to_B, 4),
                round(tm_B_to_A, 4),
                round(aln_A_to_B, 4),
                round(aln_B_to_A, 4)
            ])

            count_pairs += 1

            if count_pairs % 1_000_000 == 0:
                print(f"Bidirectional TM pairs processed: {count_pairs:,}")

    print("Finished bidirectional TM-passing pairs:", count_pairs)


# ============================================================
# RUN RBH (ALL + TM-FILTERED) AND BIDIRECTIONAL
# ============================================================

compute_rbh(
    best_hit_all,
    PAIRWISE_MATCHES_ALL,
    PROTEOME_OUT_ALL,
    label="ALL"
)

compute_rbh(
    best_hit_tm,
    PAIRWISE_MATCHES_TM,
    PROTEOME_OUT_FILTERED,
    label=f"TM≥{TM_THRESHOLD}"
)

compute_bidirectional_tm_pairs(
    tm_pass_hits,
    PAIRWISE_TM_BIDIRECTIONAL
)

# ============================================================
# CLEANUP
# ============================================================

if os.path.exists(RAW_TSV):
    os.remove(RAW_TSV)

print("\nDONE")
