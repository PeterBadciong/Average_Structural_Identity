#!/usr/bin/env python3
import os
import sys
import subprocess
import itertools
import pandas as pd
import argparse
from multiprocessing import Pool


# ARGPARSE


parser = argparse.ArgumentParser(description="Parallel Pairwise RBH AAI Calculator")

parser.add_argument("--faa-folder", required=True,
                    help="Folder containing .faa files")

parser.add_argument("--outdir", required=True,
                    help="Output directory for pairwise RBH AAI results")

parser.add_argument("--threads", type=int, default=16,
                    help="Number of parallel DIAMOND jobs (each uses 1 thread)")

parser.add_argument("--aai-cutoff", type=float, default=30.0,
                    help="Minimum AAI required to keep an RBH")

parser.add_argument("--cov-cutoff", type=float, default=0.0,
                    help="Minimum average coverage required to keep an RBH")

args = parser.parse_args()

FAADIR = args.faa_folder
OUTDIR = args.outdir
PARALLEL_JOBS = args.threads     # <-- now controls parallelism
MIN_AAI = args.aai_cutoff
MIN_COV = args.cov_cutoff


# UTILS


def run(cmd):
    print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def load_diamond_tsv(path):
    cols = [
        "qseqid","sseqid","pident",
        "mismatch","gapopen",
        "qstart","qend","sstart","send",
        "evalue","bitscore"
    ]
    return pd.read_csv(path, sep="\t", names=cols)

def load_protein_lengths(faa_folder):
    lengths = {}
    for fname in os.listdir(faa_folder):
        if not fname.endswith(".faa"):
            continue
        path = os.path.join(faa_folder, fname)
        with open(path) as f:
            seq_id = None
            seq = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if seq_id is not None:
                        lengths[seq_id] = len("".join(seq))
                    seq_id = line[1:].split()[0]
                    seq = []
                else:
                    seq.append(line)
            if seq_id is not None:
                lengths[seq_id] = len("".join(seq))
    return lengths

protein_lengths = load_protein_lengths(FAADIR)


# RBH COMPUTATION


def compute_rbh(tsv_AB, tsv_BA):
    AB = load_diamond_tsv(tsv_AB)
    BA = load_diamond_tsv(tsv_BA)

    best_AB = AB.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    best_BA = BA.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")

    rbh = best_AB.merge(
        best_BA,
        left_on=["qseqid", "sseqid"],
        right_on=["sseqid", "qseqid"],
        suffixes=("_AB", "_BA")
    )

    rbh["AAI"] = (rbh["pident_AB"] + rbh["pident_BA"]) / 2
    rbh["Alignment_length"] = (rbh["qend_AB"] - rbh["qstart_AB"]).abs() + 1

    def qcov(row):
        qlen = protein_lengths.get(row["qseqid_AB"], 1)
        return row["Alignment_length"] / qlen

    def scov(row):
        slen = protein_lengths.get(row["sseqid_AB"], 1)
        return row["Alignment_length"] / slen

    rbh["Query_coverage"] = rbh.apply(qcov, axis=1)
    rbh["Subject_coverage"] = rbh.apply(scov, axis=1)
    rbh["Avg_coverage"] = (rbh["Query_coverage"] + rbh["Subject_coverage"]) / 2

    return rbh


# PARALLEL JOB FUNCTION


def process_pair(pair):
    A, B = pair
    Aname = os.path.basename(A).replace(".faa", "")
    Bname = os.path.basename(B).replace(".faa", "")

    out_AB = os.path.join(OUTDIR, f"{Aname}_vs_{Bname}.tsv")
    out_BA = os.path.join(OUTDIR, f"{Bname}_vs_{Aname}.tsv")

    # A → B
    run(
        f"diamond blastp -q {A} -d {OUTDIR}/{Bname}.dmnd "
        f"-o {out_AB} -f 6 qseqid sseqid pident mismatch gapopen "
        f"qstart qend sstart send evalue bitscore --threads 1"
    )

    # B → A
    run(
        f"diamond blastp -q {B} -d {OUTDIR}/{Aname}.dmnd "
        f"-o {out_BA} -f 6 qseqid sseqid pident mismatch gapopen "
        f"qstart qend sstart send evalue bitscore --threads 1"
    )

    rbh = compute_rbh(out_AB, out_BA)

    # Apply cutoffs
    rbh = rbh[(rbh["AAI"] >= MIN_AAI) & (rbh["Avg_coverage"] >= MIN_COV)]

    # Cleanup DIAMOND outputs
    os.remove(out_AB)
    os.remove(out_BA)

    summary = {
        "GenomeA": Aname,
        "GenomeB": Bname,
        "Total_RBH": rbh.shape[0],
        "Mean_AAI": rbh["AAI"].mean() if rbh.shape[0] > 0 else 0
    }

    details = []
    for _, row in rbh.iterrows():
        details.append({
            "GenomeA": Aname,
            "GenomeB": Bname,
            "ProteinA": row["qseqid_AB"],
            "ProteinB": row["sseqid_AB"],
            "AAI": row["AAI"],
            "Alignment_length": row["Alignment_length"],
            "Query_coverage": row["Query_coverage"],
            "Subject_coverage": row["Subject_coverage"],
            "Avg_coverage": row["Avg_coverage"]
        })

    return summary, details


# MAIN


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    faa_files = sorted([
        os.path.join(FAADIR, f)
        for f in os.listdir(FAADIR)
        if f.endswith(".faa")
    ])

    print(f"Found {len(faa_files)} genomes")

    # Build DIAMOND DBs
    for faa in faa_files:
        base = os.path.basename(faa).replace(".faa", "")
        db_path = os.path.join(OUTDIR, f"{base}.dmnd")
        if not os.path.exists(db_path):
            run(f"diamond makedb --in {faa} -d {OUTDIR}/{base}")
        else:
            print(f"[SKIP] DB exists for {base}")

    pairs = list(itertools.combinations(faa_files, 2))

    print(f"Running {len(pairs)} pairwise comparisons using {PARALLEL_JOBS} parallel jobs")

    with Pool(PARALLEL_JOBS) as pool:
        results = pool.map(process_pair, pairs)

    summary_rows = []
    all_rbh_rows = []

    for summary, details in results:
        summary_rows.append(summary)
        all_rbh_rows.extend(details)

    pd.DataFrame(summary_rows).to_csv(os.path.join(OUTDIR, "pairwise_summary.csv"), index=False)
    pd.DataFrame(all_rbh_rows).to_csv(os.path.join(OUTDIR, "all_rbh_proteins.csv"), index=False)

    # Cleanup DBs
    for faa in faa_files:
        base = os.path.basename(faa).replace(".faa", "")
        db_path = os.path.join(OUTDIR, f"{base}.dmnd")
        if os.path.exists(db_path):
            os.remove(db_path)
            print(f"[CLEAN] Removed {db_path}")

    print("\n=== DONE ===")
    print(f"Summary table: {OUTDIR}/pairwise_summary.csv")
    print(f"All RBH proteins: {OUTDIR}/all_rbh_proteins.csv")


if __name__ == "__main__":
    main()
