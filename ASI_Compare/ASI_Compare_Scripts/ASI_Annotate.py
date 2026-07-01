import csv
import subprocess
from pathlib import Path
import os
from concurrent.futures import ProcessPoolExecutor
import argparse

# ============================================================
# ARGPARSE (RETROFITTED)
# ============================================================

parser = argparse.ArgumentParser(description="Legacy AAI + Annotation Pipeline (Argparse Retrofitted)")

parser.add_argument("--input-tsv", required=True,
                    help="TSV file from Script 1")

parser.add_argument("--faa-folder", required=True,
                    help="Folder containing .faa files")

parser.add_argument("--diamond", required=True,
                    help="Path to DIAMOND executable")

parser.add_argument("--hmmsearch", required=True,
                    help="Path to hmmsearch executable")

parser.add_argument("--merged-hmm", required=True,
                    help="Merged HMM database")

parser.add_argument("--merged-annot", required=True,
                    help="Merged annotation table")

parser.add_argument("--threads", type=int, default=32,
                    help="CPU threads")

parser.add_argument("--evalue", type=float, default=1e-5,
                    help="E-value cutoff")

parser.add_argument("--run-annotation", action="store_true",
                    help="Enable HMM annotation")

parser.add_argument("--workdir", required=True,
                    help="Working/output directory")

args = parser.parse_args()

# ============================================================
# ASSIGN ARGUMENTS TO ORIGINAL VARIABLE NAMES
# ============================================================

INPUT_TSV = args.input_tsv
FAA_FOLDER = args.faa_folder

DIAMOND = args.diamond
HMMSEARCH = args.hmmsearch

MERGED_HMM = args.merged_hmm
MERGED_ANNOT = args.merged_annot

EVALUE = args.evalue
RUN_ANNOTATION = args.run_annotation

WORK_DIR = args.workdir
TMP_DIR = f"{WORK_DIR}/tmp"

THREADS = args.threads
MAX_WORKERS = 8
CHUNK_SIZE = 250000

os.environ["OMP_NUM_THREADS"] = str(THREADS)

Path(WORK_DIR).mkdir(parents=True, exist_ok=True)
Path(TMP_DIR).mkdir(parents=True, exist_ok=True)

# ============================================================
# NORMALIZATION
# ============================================================

def norm_id(x):
    if x is None:
        return None
    return str(x).strip().split()[0]

# ============================================================
# FASTA UTILITIES
# ============================================================

def load_fasta(path):
    seqs = {}
    with open(path) as f:
        h = None
        s = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if h:
                    seqs[h] = "".join(s)
                h = norm_id(line[1:])
                s = []
            else:
                s.append(line)

        if h:
            seqs[h] = "".join(s)

    return seqs


def write_fasta(seqs, out_file):
    with open(out_file, "w") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n{v}\n")

# ============================================================
# TSV PARSING
# ============================================================

def parse_tsv(tsv_file):
    pairs = []
    proteins = set()

    with open(tsv_file, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")

        for row in r:
            try:
                g1 = row["Genome_A"]
                g2 = row["Genome_B"]

                p1 = norm_id(row["Protein_A_full_id"])
                p2 = norm_id(row["Protein_B_full_id"])

                tm = float(row["TM_score"])
                altm = float(row["AlnTM_score"])
            except:
                continue

            pairs.append((g1, g2, p1, p2, tm, altm))
            proteins.update([p1, p2])

    return pairs, proteins

# ============================================================
# CHUNKING
# ============================================================

def chunk_dict(d, size):
    items = list(d.items())
    for i in range(0, len(items), size):
        yield dict(items[i:i + size])


def write_chunk_fasta(seqs, chunk_id):
    path = Path(TMP_DIR) / f"chunk_{chunk_id}.faa"
    write_fasta(seqs, path)
    return path

# ============================================================
# HMM WORKER
# ============================================================

def run_hmm_chunk(args):
    faa_file, chunk_id = args

    out_tbl = Path(TMP_DIR) / f"hmm_{chunk_id}.tbl"

    cmd = [
        HMMSEARCH,
        "--cpu", str(THREADS),
        "--tblout", str(out_tbl),
        "-E", str(EVALUE),
        str(MERGED_HMM),
        str(faa_file),
    ]

    subprocess.run(cmd, check=True)
    return out_tbl

# ============================================================
# HMM PARSER
# ============================================================

def parse_hmm_tbl(tbl_file):
    best_hits = {}

    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            p = line.split()
            if len(p) < 18:
                continue

            try:
                protein = norm_id(p[0])
                model = norm_id(p[2])
                evalue = float(p[4])
                bitscore = float(p[5])
            except:
                continue

            if evalue > EVALUE:
                continue

            prev = best_hits.get(protein)

            if prev is None or bitscore > prev["bitscore"]:
                best_hits[protein] = {
                    "model": model,
                    "bitscore": bitscore,
                    "evalue": evalue,
                }

    return best_hits

# ============================================================
# ANNOTATION TABLE
# ============================================================

def load_annotation_table(path):
    ann = {}

    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")

        for row in r:
            try:
                model = norm_id(row["model_id"])
                ann[model] = f"{row['source']}:{row['annotation']}"
            except:
                continue

    return ann

# ============================================================
# MAIN
# ============================================================

def main():

    print("Parsing TSV...")
    pairs, proteins = parse_tsv(INPUT_TSV)
    print(f"Proteins: {len(proteins)}")

    # --------------------------------------------------------
    # LOAD SEQUENCES
    # --------------------------------------------------------

    seqs = {}
    for faa in Path(FAA_FOLDER).glob("*.faa"):
        seqs.update(load_fasta(faa))

    seqs = {k: v for k, v in seqs.items() if k in proteins}
    print(f"Loaded sequences: {len(seqs)}")

    # ========================================================
    # OPTIONAL ANNOTATION
    # ========================================================

    if RUN_ANNOTATION:
        print("Running annotation (HMMER + annotation table)...")

        chunks = []
        for i, chunk in enumerate(chunk_dict(seqs, CHUNK_SIZE)):
            fasta_path = write_chunk_fasta(chunk, i)
            chunks.append((fasta_path, i))

        print(f"Running {len(chunks)} HMMER chunks...")

        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
            tbl_files = list(ex.map(run_hmm_chunk, chunks))

        best_hits = {}
        for tbl in tbl_files:
            best_hits.update(parse_hmm_tbl(tbl))

        ann_map = load_annotation_table(MERGED_ANNOT)

        hmm_hits = {}
        hmm_bitscore = {}
        hmm_evalue = {}
        annotations = {}

        for p in proteins:
            hit = best_hits.get(p)

            if not hit:
                hmm_hits[p] = "NA"
                hmm_bitscore[p] = "NA"
                hmm_evalue[p] = "NA"
                annotations[p] = "NA"
            else:
                model = norm_id(hit["model"])

                hmm_hits[p] = model
                hmm_bitscore[p] = hit["bitscore"]
                hmm_evalue[p] = hit["evalue"]
                annotations[p] = ann_map.get(model, model)

        print(f"Annotated proteins: {sum(v != 'NA' for v in annotations.values())}")

    else:
        print("Skipping annotation step.")
        hmm_hits = {p: "NA" for p in proteins}
        hmm_bitscore = {p: "NA" for p in proteins}
        hmm_evalue = {p: "NA" for p in proteins}
        annotations = {p: "NA" for p in proteins}

    # ========================================================
    # FAST DIAMOND AAI (OPTIMIZED)
    # ========================================================

    print("\nRunning optimized DIAMOND AAI...")

    all_faa = Path(TMP_DIR) / "all.faa"
    write_fasta(seqs, all_faa)

    db = Path(TMP_DIR) / "all_db"

    subprocess.run([
        DIAMOND, "makedb",
        "--in", str(all_faa),
        "-d", str(db),
    ], check=True)

    pair_map = {}
    for _, _, p1, p2, _, _ in pairs:
        if p1 in seqs and p2 in seqs:
            pair_map.setdefault(p1, set()).add(p2)

    CHUNK_SIZE_DIAMOND = 20000
    protein_list = list(pair_map.keys())

    aai = {}

    for i in range(0, len(protein_list), CHUNK_SIZE_DIAMOND):

        chunk = protein_list[i:i + CHUNK_SIZE_DIAMOND]

        query_faa = Path(TMP_DIR) / f"aai_chunk_{i}.faa"
        write_fasta({p: seqs[p] for p in chunk}, query_faa)

        out_file = Path(TMP_DIR) / f"aai_chunk_{i}.tsv"

        subprocess.run([
            DIAMOND, "blastp",
            "--sensitive",
            "--comp-based-stats", "1",
            "--masking", "1",
            "-q", str(query_faa),
            "-d", str(db),
            "-o", str(out_file),
            "--outfmt", "6", "qseqid", "sseqid", "pident", "bitscore",
            "--max-target-seqs", "50",
            "--threads", str(THREADS * MAX_WORKERS),
        ], check=True)

        best = {}

        with open(out_file) as f:
            for line in f:
                q, s, pid, bits = line.strip().split("\t")

                q = norm_id(q)
                s = norm_id(s)

                if q not in pair_map or s not in pair_map[q]:
                    continue

                try:
                    pid = float(pid)
                    bits = float(bits)
                except:
                    continue

                key = (q, s)

                if key not in best or bits > best[key][1]:
                    best[key] = (pid, bits)

        for (q, s), (pid, _) in best.items():
            aai[(q, s)] = pid
            aai[(s, q)] = pid

    print(f"AAI computed for {len(aai)//2} protein pairs")

    # ========================================================
    # OUTPUT
    # ========================================================

    out_csv = Path(WORK_DIR) / "ALL_GENOMES_FINAL.csv"

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)

        w.writerow([
            "Genome1", "Genome2",
            "Protein1", "Protein2",
            "TM", "AlnTM",
            "AAI",
            "HMM1", "HMM2",
            "Protein1Bitscore", "Protein2Bitscore",
            "Protein1Evalue", "Protein2Evalue",
            "Annotation1", "Annotation2",
        ])

        for g1, g2, p1, p2, tm, altm in pairs:
            w.writerow([
                g1, g2,
                p1, p2,
                tm, altm,
                aai.get((p1, p2), "NA"),
                hmm_hits.get(p1, "NA"),
                hmm_hits.get(p2, "NA"),
                hmm_bitscore.get(p1, "NA"),
                hmm_bitscore.get(p2, "NA"),
                hmm_evalue.get(p1, "NA"),
                hmm_evalue.get(p2, "NA"),
                annotations.get(p1, "NA"),
                annotations.get(p2, "NA"),
            ])

    print("\nDONE:", out_csv)


if __name__ == "__main__":
    main()
