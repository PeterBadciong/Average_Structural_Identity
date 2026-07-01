import csv
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import argparse


# ARGPARSE (RETROFITTED)


parser = argparse.ArgumentParser(description="Annotation Comparison Classifier (Part 3)")

parser.add_argument("--input-csv", required=True,
                    help="ALL_GENOMES_FINAL.csv from Part 2")

parser.add_argument("--output-dir", required=True,
                    help="Directory to write classification outputs")

parser.add_argument("--threads", type=int, default=32,
                    help="Number of worker processes")

parser.add_argument("--chunk-size", type=int, default=200000,
                    help="Rows per chunk")

args = parser.parse_args()


# ASSIGN ARGUMENTS TO ORIGINAL VARIABLE NAMES


INPUT_CSV = args.input_csv
OUTPUT_DIR = Path(args.output_dir)
THREADS = args.threads
CHUNK_SIZE = args.chunk_size

OUTPUT_DIR.mkdir(exist_ok=True)

HYPOTHETICAL_OUT = OUTPUT_DIR / "Hypothetical_Matches.csv"
BOTH_HYPOTHETICAL_OUT = OUTPUT_DIR / "Both_Hypothetical_Matches.csv"
EXACT_OUT = OUTPUT_DIR / "Exact_Matches.csv"
PARTIAL_OUT = OUTPUT_DIR / "Partial_Matches.csv"
UNMATCHED_OUT = OUTPUT_DIR / "Unmatched_Matches.csv"
UNANNOTATED_OUT = OUTPUT_DIR / "Unannotated_Matches.csv"

STATS_OUT = OUTPUT_DIR / "Category_Stats.csv"


# WORD EQUIVALENCE SYSTEM


WORD_EQUIVALENCE = {
    "head": "capsid",
    "anti-termination": "antitermination",
}

def canonicalize(word: str):
    return WORD_EQUIVALENCE.get(word, word)


# CLEANING RULES


IGNORE_WORDS = {
    "protein",
    "family",
    "multispecies",
    "subunit",
    "putative",
    "major",
}

ACCESSION_PATTERN = re.compile(r"\b[A-Z0-9]{3,}\d+\S*\b")

def clean_annotation(annotation: str) -> str:
    annotation = annotation.strip()
    annotation = re.sub(r"^[A-Za-z0-9_]+:", "", annotation)
    return annotation

def extract_words(annotation: str):
    annotation = clean_annotation(annotation)
    annotation = ACCESSION_PATTERN.sub(" ", annotation)

    words = re.findall(r"[A-Za-z\-]+", annotation.lower())

    return {
        canonicalize(w)
        for w in words
        if w not in IGNORE_WORDS and len(w) >= 3
    }

def normalize(annotation: str):
    return frozenset(extract_words(annotation))

def is_unknown(annotation: str):
    ann = annotation.lower()
    return (
        "hypothetical" in ann
        or ann.strip() == "phrog:na"
        or ann.strip() == "na"
    )

def is_phrog_na(annotation: str):
    return "phrog:na" in annotation.lower()

def is_na(annotation: str):
    if annotation is None:
        return True
    ann = annotation.strip().lower()
    return ann in {"na", "n/a", "", "nan"}

def contains_match(words1, words2):
    for w1 in words1:
        for w2 in words2:
            if w1 in w2 or w2 in w1:
                return True
    return False

def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# PARALLEL CLASSIFICATION


def process_chunk(rows):

    hypothetical_rows = []
    both_hypothetical_rows = []
    exact_rows = []
    partial_rows = []
    unmatched_rows = []
    unannotated_rows = []

    for row in rows:

        ann1 = row["Annotation1"].strip()
        ann2 = row["Annotation2"].strip()

        # PRIORITY -1: BOTH UNANNOTATED
        if is_na(ann1) and is_na(ann2):
            unannotated_rows.append(row)
            continue

        # PRIORITY 0: PHROG:NA RULE
        ann1_phrog_na = is_phrog_na(ann1)
        ann2_phrog_na = is_phrog_na(ann2)

        if ann1_phrog_na and ann2_phrog_na:

            if row["HMM1"] == row["HMM2"]:
                exact_rows.append(row)
                continue

            both_hypothetical_rows.append(row)
            continue

        # NORMALIZATION
        norm1 = normalize(ann1)
        norm2 = normalize(ann2)

        ann1_unknown = is_unknown(ann1)
        ann2_unknown = is_unknown(ann2)

        # PRIORITY 1: EXACT MATCHES
        if norm1 == norm2:
            exact_rows.append(row)
            continue

        if len(norm1) == 0 and len(norm2) == 0:
            exact_rows.append(row)
            continue

        # PRIORITY 2: BOTH UNKNOWN
        if ann1_unknown and ann2_unknown:
            both_hypothetical_rows.append(row)
            continue

        # PRIORITY 3: ONE UNKNOWN
        if ann1_unknown != ann2_unknown:
            hypothetical_rows.append(row)
            continue

        # PRIORITY 4: PARTIAL MATCH
        if contains_match(norm1, norm2):
            partial_rows.append(row)
            continue

        # PRIORITY 5: UNMATCHED
        unmatched_rows.append(row)

    return (
        hypothetical_rows,
        both_hypothetical_rows,
        exact_rows,
        partial_rows,
        unmatched_rows,
        unannotated_rows,
    )


# CHUNKING


def chunked(iterable, size):
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


# MAIN EXECUTION


hypothetical_rows = []
both_hypothetical_rows = []
exact_rows = []
partial_rows = []
unmatched_rows = []
unannotated_rows = []

with open(INPUT_CSV, "r", newline="") as infile:

    reader = csv.DictReader(infile)
    fieldnames = reader.fieldnames

    print(f"Running with {THREADS} processes")

    with ProcessPoolExecutor(max_workers=THREADS) as executor:

        for result in executor.map(process_chunk, chunked(reader, CHUNK_SIZE)):

            h, bh, e, p, u, ua = result

            hypothetical_rows.extend(h)
            both_hypothetical_rows.extend(bh)
            exact_rows.extend(e)
            partial_rows.extend(p)
            unmatched_rows.extend(u)
            unannotated_rows.extend(ua)


# SAVE CSVs


write_csv(HYPOTHETICAL_OUT, hypothetical_rows, fieldnames)
write_csv(BOTH_HYPOTHETICAL_OUT, both_hypothetical_rows, fieldnames)
write_csv(EXACT_OUT, exact_rows, fieldnames)
write_csv(PARTIAL_OUT, partial_rows, fieldnames)
write_csv(UNMATCHED_OUT, unmatched_rows, fieldnames)
write_csv(UNANNOTATED_OUT, unannotated_rows, fieldnames)

print("Category CSVs written.")


# STATISTICS


FILES = {
    "Exact": exact_rows,
    "Both_Hypothetical": both_hypothetical_rows,
    "Hypothetical": hypothetical_rows,
    "Partial": partial_rows,
    "Unmatched": unmatched_rows,
    "Unannotated": unannotated_rows,
}

total = sum(len(v) for v in FILES.values())

rows = []

for category, data in FILES.items():

    df = pd.DataFrame(data)

    if df.empty:
        continue

    df["TM"] = pd.to_numeric(df["TM"], errors="coerce")
    df["AAI"] = pd.to_numeric(df["AAI"], errors="coerce") / 100.0

    tm = df["TM"].dropna()
    aai = df["AAI"]

    n = len(df)
    pct_total = (n / total * 100) if total > 0 else 0

    aai_na = aai.isna().sum()
    aai_na_pct = (aai_na / n * 100) if n > 0 else 0

    rows.append({
        "Category": category,
        "N": n,
        "Percent_of_Total": round(pct_total, 2),

        "AAI_NA": int(aai_na),
        "AAI_NA_Percent_of_Category": round(aai_na_pct, 2),

        "TM_mean": tm.mean(),
        "TM_median": tm.median(),
        "TM_min": tm.min(),
        "TM_max": tm.max(),
        "TM_sd": tm.std(),

        "AAI_mean": aai.dropna().mean(),
        "AAI_median": aai.dropna().median(),
        "AAI_min": aai.dropna().min(),
        "AAI_max": aai.dropna().max(),
        "AAI_sd": aai.dropna().std(),
    })

stats_df = pd.DataFrame(rows)
stats_df.to_csv(STATS_OUT, index=False)

print(f"Stats written to: {STATS_OUT}")


# SUMMARY


print("\nDONE")
print(f"Total rows: {total}")
print(f"Exact: {len(exact_rows)}")
print(f"Both Hypothetical: {len(both_hypothetical_rows)}")
print(f"Hypothetical (one-sided): {len(hypothetical_rows)}")
print(f"Partial: {len(partial_rows)}")
print(f"Unmatched: {len(unmatched_rows)}")
print(f"Unannotated: {len(unannotated_rows)}")
