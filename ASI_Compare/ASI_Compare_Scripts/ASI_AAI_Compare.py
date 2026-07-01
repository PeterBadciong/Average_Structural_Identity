#!/usr/bin/env python3
import csv
from pathlib import Path
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# ============================================================
# ARGPARSE
# ============================================================

parser = argparse.ArgumentParser(
    description="Part 5: Compare Foldseek RBHs and AAI RBHs with PDB diagnostics + pLDDT + bidirectional TM"
)

parser.add_argument("--tm-tsv", required=True,
                    help="TMfiltered Foldseek RBH file")

parser.add_argument("--tm-bidir", required=True,
                    help="pairwise_TM_bidirectional.tsv from Part 1")

parser.add_argument("--aai-csv", required=True,
                    help="all_rbh_proteins.csv from Part 4")

parser.add_argument("--pdb-dir", required=True,
                    help="Directory containing PDB files")

parser.add_argument("--map-file", required=True,
                    help="protein_map.tsv from Part 1")

parser.add_argument("--outdir", required=True)

parser.add_argument("--compute-plddt", action="store_true")

parser.add_argument("--threads", type=int, default=64,
                    help="Number of threads for pLDDT and matching")

args = parser.parse_args()

TM_TSV = args.tm_tsv
TM_BIDIR = args.tm_bidir
AAI_CSV = args.aai_csv
PDB_DIR = Path(args.pdb_dir)
MAP_FILE = Path(args.map_file)
OUTDIR = Path(args.outdir)
OUTDIR.mkdir(parents=True, exist_ok=True)

N_THREADS = max(1, args.threads)

# ============================================================
# HELPERS
# ============================================================

def chunked_iterable(iterable, size=10000):
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

# ============================================================
# OUTPUT PATHS
# ============================================================

FOLDSEEK_ONLY_OUT = OUTDIR / "Foldseek_Exclusive.csv"
AAI_ONLY_OUT = OUTDIR / "AAI_Exclusive.csv"
BOTH_OUT = OUTDIR / "Foldseek_and_AAI.csv"
STATS_OUT = OUTDIR / "Part5_Stats.csv"
GENOME_SUMMARY_OUT = OUTDIR / "Genome_Pair_Summary.csv"

# ============================================================
# LOAD PROTEIN→GENOME MAP
# ============================================================

protein_to_genome = {}
with open(MAP_FILE) as f:
    for line in f:
        p, g = line.strip().split("\t")
        protein_to_genome[p] = g

# ============================================================
# LOAD PARALOG INFORMATION
# ============================================================

PARALOG_FILE = OUTDIR.parent / "paralogs" / "paralogs_TM0.8.tsv"
paralog_members = set()

if PARALOG_FILE.exists():
    with open(PARALOG_FILE) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            _, pA, pB, *_ = parts
            paralog_members.add(pA)
            paralog_members.add(pB)
else:
    print(f"[WARNING] No paralog file found at {PARALOG_FILE}")

# ============================================================
# pLDDT EXTRACTION (PARALLELIZED)
# ============================================================

def compute_plddt(pdb_path: Path):
    if not pdb_path.exists():
        return None
    scores = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    scores.append(float(line[60:66]))
                except:
                    pass
    return sum(scores) / len(scores) if scores else None

def pdb_info_single(protein_id: str):
    pdb_path = PDB_DIR / f"{protein_id}.pdb"
    exists = pdb_path.exists()
    plddt = compute_plddt(pdb_path) if exists and args.compute_plddt else None
    return protein_id, exists, plddt

# ============================================================
# NORMALIZATION
# ============================================================

def normalize_pair(pA, pB):
    gA = protein_to_genome[pA]
    gB = protein_to_genome[pB]

    if gA < gB:
        return gA, gB, pA, pB
    elif gB < gA:
        return gB, gA, pB, pA
    else:
        return gA, gB, min(pA, pB), max(pA, pB)

# ============================================================
# LOAD TM-FILTERED FOLDSEEK RBHs (PARALLEL)
# ============================================================

foldseek_pairs = {}

def process_foldseek_chunk(chunk):
    local = {}
    for row in chunk:
        pA = row["Protein_A_full_id"].strip()
        pB = row["Protein_B_full_id"].strip()

        if pA not in protein_to_genome or pB not in protein_to_genome:
            continue

        gA, gB, pA_norm, pB_norm = normalize_pair(pA, pB)

        try:
            tm = float(row["TM_score"])
            alntm = float(row["AlnTM_score"])
        except:
            continue

        local[(gA, gB, pA_norm, pB_norm)] = (tm, alntm)
    return local

with open(TM_TSV, newline="") as f:
    r = list(csv.DictReader(f, delimiter="\t"))

with ThreadPoolExecutor(max_workers=N_THREADS) as ex:
    futures = [ex.submit(process_foldseek_chunk, chunk)
               for chunk in chunked_iterable(r, 10000)]
    for fut in as_completed(futures):
        foldseek_pairs.update(fut.result())

# ============================================================
# LOAD BIDIRECTIONAL TM (PARALLEL)
# ============================================================

bidir_tm = {}

def process_bidir_chunk(chunk):
    local = {}
    for row in chunk:
        pA = row["Protein_A_full_id"].strip()
        pB = row["Protein_B_full_id"].strip()

        if pA not in protein_to_genome or pB not in protein_to_genome:
            continue

        gA, gB, pA_norm, pB_norm = normalize_pair(pA, pB)

        try:
            tmAB = float(row["TM_A_to_B"])
            tmBA = float(row["TM_B_to_A"])
            avg_tm = (tmAB + tmBA) / 2.0
        except:
            avg_tm = "NA"

        local[(gA, gB, pA_norm, pB_norm)] = avg_tm
    return local

with open(TM_BIDIR, newline="") as f:
    r = list(csv.DictReader(f, delimiter="\t"))

with ThreadPoolExecutor(max_workers=N_THREADS) as ex:
    futures = [ex.submit(process_bidir_chunk, chunk)
               for chunk in chunked_iterable(r, 10000)]
    for fut in as_completed(futures):
        bidir_tm.update(fut.result())

# ============================================================
# LOAD NO-CUTOFF FOLDSEEK MATCHES (ALL.tsv, PARALLEL)
# ============================================================

NO_CUTOFF_FILE = TM_TSV.replace("TMfiltered", "ALL")
no_cutoff_scores = {}

def process_no_cutoff_chunk(chunk):
    local = {}
    for row in chunk:
        pA = row["Protein_A_full_id"].strip()
        pB = row["Protein_B_full_id"].strip()

        if pA not in protein_to_genome or pB not in protein_to_genome:
            continue

        gA, gB, pA_norm, pB_norm = normalize_pair(pA, pB)

        try:
            tm = float(row["TM_score"])
        except:
            tm = 0.0

        local[(gA, gB, pA_norm, pB_norm)] = tm
    return local

with open(NO_CUTOFF_FILE, newline="") as f:
    r = list(csv.DictReader(f, delimiter="\t"))

with ThreadPoolExecutor(max_workers=N_THREADS) as ex:
    futures = [ex.submit(process_no_cutoff_chunk, chunk)
               for chunk in chunked_iterable(r, 10000)]
    for fut in as_completed(futures):
        no_cutoff_scores.update(fut.result())

# ============================================================
# LOAD AAI RBHs (WITH COVERAGE, PARALLEL)
# ============================================================

aai_pairs = {}

def process_aai_chunk(chunk):
    local = {}
    for row in chunk:
        pA = row["ProteinA"].strip()
        pB = row["ProteinB"].strip()

        if pA not in protein_to_genome or pB not in protein_to_genome:
            continue

        gA, gB, pA_norm, pB_norm = normalize_pair(pA, pB)

        try:
            aai = float(row["AAI"])
            aln_len = float(row["Alignment_length"])
            qcov = float(row["Query_coverage"])
            scov = float(row["Subject_coverage"])
            avgcov = float(row["Avg_coverage"])
        except:
            continue

        local[(gA, gB, pA_norm, pB_norm)] = {
            "AAI": aai,
            "Alignment_length": aln_len,
            "Query_coverage": qcov,
            "Subject_coverage": scov,
            "Avg_coverage": avgcov
        }
    return local

with open(AAI_CSV, newline="") as f:
    r = list(csv.DictReader(f))

with ThreadPoolExecutor(max_workers=N_THREADS) as ex:
    futures = [ex.submit(process_aai_chunk, chunk)
               for chunk in chunked_iterable(r, 10000)]
    for fut in as_completed(futures):
        aai_pairs.update(fut.result())

# ============================================================
# PRECOMPUTE PDB INFO IN PARALLEL
# ============================================================

all_proteins = set()
for (_, _, pA, pB) in foldseek_pairs.keys():
    all_proteins.add(pA)
    all_proteins.add(pB)
for (_, _, pA, pB) in aai_pairs.keys():
    all_proteins.add(pA)
    all_proteins.add(pB)

pdb_cache = {}

if args.compute_plddt and all_proteins:
    print(f"Computing pLDDT for {len(all_proteins)} proteins using {N_THREADS} threads...")
else:
    print(f"Caching PDB existence for {len(all_proteins)} proteins using {N_THREADS} threads...")

def worker(protein_id):
    return pdb_info_single(protein_id)

with ThreadPoolExecutor(max_workers=N_THREADS) as ex:
    futures = {ex.submit(worker, pid): pid for pid in all_proteins}
    for fut in as_completed(futures):
        pid, exists, plddt = fut.result()
        pdb_cache[pid] = (exists, plddt)

def pdb_info(protein_id):
    return pdb_cache.get(protein_id, (False, None))

# ============================================================
# SETS AND CATEGORIES
# ============================================================

fs_keys = set(foldseek_pairs.keys())
aai_keys = set(aai_pairs.keys())

both_keys = fs_keys & aai_keys
fs_only_keys = fs_keys - aai_keys
aai_only_keys = aai_keys - fs_keys

# ============================================================
# WRITE CATEGORY CSVs
# ============================================================

def write_category_csv(path, rows, extra_fields=None):
    base_fields = [
        "GenomeA", "GenomeB",
        "ProteinA", "ProteinB",
        "TM", "AlnTM", "AAI",
        "Alignment_length", "Query_coverage", "Subject_coverage", "Avg_coverage",
        "ProteinA_in_PDB", "ProteinB_in_PDB",
        "pLDDT_A", "pLDDT_B",
        "Avg_pLDDT"
    ]
    if extra_fields:
        base_fields += extra_fields

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=base_fields)
        w.writeheader()
        w.writerows(rows)

# Foldseek-exclusive
foldseek_only_rows = []
for gA, gB, pA, pB in fs_only_keys:
    tm, alntm = foldseek_pairs[(gA, gB, pA, pB)]
    A_exists, A_plddt = pdb_info(pA)
    B_exists, B_plddt = pdb_info(pB)

    avg_plddt = None
    if A_plddt is not None and B_plddt is not None:
        avg_plddt = (A_plddt + B_plddt) / 2.0

    foldseek_only_rows.append({
        "GenomeA": gA, "GenomeB": gB,
        "ProteinA": pA, "ProteinB": pB,
        "TM": tm, "AlnTM": alntm, "AAI": "NA",
        "Alignment_length": "NA",
        "Query_coverage": "NA",
        "Subject_coverage": "NA",
        "Avg_coverage": "NA",
        "ProteinA_in_PDB": A_exists, "ProteinB_in_PDB": B_exists,
        "pLDDT_A": A_plddt, "pLDDT_B": B_plddt,
        "Avg_pLDDT": avg_plddt
    })

write_category_csv(FOLDSEEK_ONLY_OUT, foldseek_only_rows)

# AAI-exclusive
aai_only_rows = []
for gA, gB, pA, pB in aai_only_keys:
    info = aai_pairs[(gA, gB, pA, pB)]
    A_exists, A_plddt = pdb_info(pA)
    B_exists, B_plddt = pdb_info(pB)

    key = (gA, gB, pA, pB)
    no_cutoff_tm = no_cutoff_scores.get(key, 0.0)
    avg_bidir_tm = bidir_tm.get(key, "NA")

    avg_plddt = None
    if A_plddt is not None and B_plddt is not None:
        avg_plddt = (A_plddt + B_plddt) / 2.0

    aai_only_rows.append({
        "GenomeA": gA, "GenomeB": gB,
        "ProteinA": pA, "ProteinB": pB,
        "TM": "NA", "AlnTM": "NA", "AAI": info["AAI"],
        "Alignment_length": info["Alignment_length"],
        "Query_coverage": info["Query_coverage"],
        "Subject_coverage": info["Subject_coverage"],
        "Avg_coverage": info["Avg_coverage"],
        "ProteinA_in_PDB": A_exists, "ProteinB_in_PDB": B_exists,
        "pLDDT_A": A_plddt, "pLDDT_B": B_plddt,
        "Avg_pLDDT": avg_plddt,
        "No_TM_Cutoff_Score": no_cutoff_tm,
        "Avg_bidirectional_TM": avg_bidir_tm,
        "ProteinA_has_paralog": pA in paralog_members,
        "ProteinB_has_paralog": pB in paralog_members
    })

write_category_csv(
    AAI_ONLY_OUT,
    aai_only_rows,
    extra_fields=[
        "No_TM_Cutoff_Score",
        "Avg_bidirectional_TM",
        "ProteinA_has_paralog",
        "ProteinB_has_paralog"
    ]
)

# Both
both_rows = []
for gA, gB, pA, pB in both_keys:
    tm, alntm = foldseek_pairs[(gA, gB, pA, pB)]
    info = aai_pairs[(gA, gB, pA, pB)]
    A_exists, A_plddt = pdb_info(pA)
    B_exists, B_plddt = pdb_info(pB)

    avg_plddt = None
    if A_plddt is not None and B_plddt is not None:
        avg_plddt = (A_plddt + B_plddt) / 2.0

    both_rows.append({
        "GenomeA": gA, "GenomeB": gB,
        "ProteinA": pA, "ProteinB": pB,
        "TM": tm, "AlnTM": alntm, "AAI": info["AAI"],
        "Alignment_length": info["Alignment_length"],
        "Query_coverage": info["Query_coverage"],
        "Subject_coverage": info["Subject_coverage"],
        "Avg_coverage": info["Avg_coverage"],
        "ProteinA_in_PDB": A_exists, "ProteinB_in_PDB": B_exists,
        "pLDDT_A": A_plddt, "pLDDT_B": B_plddt,
        "Avg_pLDDT": avg_plddt
    })

write_category_csv(BOTH_OUT, both_rows)

# ============================================================
# pLDDT THRESHOLD STATS
# ============================================================

def compute_avg_plddt(row):
    A = row["pLDDT_A"]
    B = row["pLDDT_B"]
    if A is None or B is None:
        return None
    return (A + B) / 2.0

def count_thresholds(rows):
    stats = {
        ">=0.5": 0,
        ">=0.7": 0,
        ">=0.9": 0
    }
    for r in rows:
        avg = compute_avg_plddt(r)
        if avg is None:
            continue
        if avg >= 0.5: stats[">=0.5"] += 1
        if avg >= 0.7: stats[">=0.7"] += 1
        if avg >= 0.9: stats[">=0.9"] += 1
    return stats

foldseek_plddt_stats = count_thresholds(foldseek_only_rows)
aai_plddt_stats = count_thresholds(aai_only_rows)
both_plddt_stats = count_thresholds(both_rows)


# ============================================================
# STATS CSV
# ============================================================

aai_exclusive_with_paralog = sum(
    1
    for (gA, gB, pA, pB) in aai_only_keys
    if (pA in paralog_members) or (pB in paralog_members)
)

with open(STATS_OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["Category", "Count"])
    w.writeheader()
    w.writerow({"Category": "Foldseek_Exclusive", "Count": len(fs_only_keys)})
    w.writerow({"Category": "AAI_Exclusive", "Count": len(aai_only_keys)})
    w.writerow({"Category": "Both", "Count": len(both_keys)})
    w.writerow({"Category": "Total_Foldseek_RBH", "Count": len(fs_keys)})
    w.writerow({"Category": "Total_AAI_RBH", "Count": len(aai_keys)})
    w.writerow({"Category": "AAI_Exclusive_with_Paralog",
                "Count": aai_exclusive_with_paralog})
    
        # pLDDT threshold rows
    for label, count in foldseek_plddt_stats.items():
        w.writerow({"Category": f"Foldseek_Exclusive_pLDDT_{label}", "Count": count})

    for label, count in aai_plddt_stats.items():
        w.writerow({"Category": f"AAI_Exclusive_pLDDT_{label}", "Count": count})

    for label, count in both_plddt_stats.items():
        w.writerow({"Category": f"Both_pLDDT_{label}", "Count": count})


# ============================================================
# GENOME SUMMARY
# ============================================================

tm_stats = defaultdict(lambda: {"tm_sum": 0.0, "alntm_sum": 0.0, "count": 0})
aai_stats = defaultdict(lambda: {"aai_sum": 0.0, "count": 0})

for (gA, gB, pA, pB), (tm, alntm) in foldseek_pairs.items():
    tm_stats[(gA, gB)]["tm_sum"] += tm
    tm_stats[(gA, gB)]["alntm_sum"] += alntm
    tm_stats[(gA, gB)]["count"] += 1

for (gA, gB, pA, pB), info in aai_pairs.items():
    aai_stats[(gA, gB)]["aai_sum"] += info["AAI"]
    aai_stats[(gA, gB)]["count"] += 1

genome_pairs = set(tm_stats.keys()) | set(aai_stats.keys())

with open(GENOME_SUMMARY_OUT, "w", newline="") as f:
    fieldnames = [
        "GenomeA", "GenomeB",
        "TM_mean", "AlnTM_mean", "TM_RBH_count",
        "AAI_mean", "AAI_RBH_count",
        "Unique_Foldseek_RBH", "Unique_AAI_RBH",
        "Total_below_cutoff_Foldseek_RBH",
        "Count_bidirectional_TM_present"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()

    for gA, gB in sorted(genome_pairs):
        tm_info = tm_stats.get((gA, gB))
        aai_info = aai_stats.get((gA, gB))

        if tm_info and tm_info["count"] > 0:
            tm_mean = tm_info["tm_sum"] / tm_info["count"]
            alntm_mean = tm_info["alntm_sum"] / tm_info["count"]
            tm_count = tm_info["count"]
        else:
            tm_mean = "NA"
            alntm_mean = "NA"
            tm_count = 0

        if aai_info and aai_info["count"] > 0:
            aai_mean = aai_info["aai_sum"] / aai_info["count"]
            aai_count = aai_info["count"]
        else:
            aai_mean = "NA"
            aai_count = 0

        unique_fs = sum(
            1 for (GA, GB, _, _) in fs_only_keys
            if (GA, GB) == (gA, gB)
        )

        unique_aai = sum(
            1 for (GA, GB, _, _) in aai_only_keys
            if (GA, GB) == (gA, gB)
        )

        below_cutoff = sum(
            1
            for (GA, GB, pA, pB) in aai_only_keys
            if (GA, GB) == (gA, gB)
            and no_cutoff_scores.get((GA, GB, pA, pB), 0.0) > 0
        )

        count_bidir_present = sum(
            1
            for (GA, GB, pA, pB) in aai_only_keys
            if (GA, GB) == (gA, gB)
            and bidir_tm.get((GA, GB, pA, pB), "NA") != "NA"
        )

        w.writerow({
            "GenomeA": gA, "GenomeB": gB,
            "TM_mean": tm_mean, "AlnTM_mean": alntm_mean,
            "TM_RBH_count": tm_count,
            "AAI_mean": aai_mean, "AAI_RBH_count": aai_count,
            "Unique_Foldseek_RBH": unique_fs,
            "Unique_AAI_RBH": unique_aai,
            "Total_below_cutoff_Foldseek_RBH": below_cutoff,
            "Count_bidirectional_TM_present": count_bidir_present
        })

print("Part 5 complete.")
