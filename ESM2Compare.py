from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
from itertools import combinations
from tqdm import tqdm

# =========================
# User settings
# =========================
input_folder = "ESM2_NPZ"
similarity_threshold = 0.0000
output_csv = "pairwise_RBH_similarity.csv"
num_threads = 20  # adjust to your CPU cores
save_every = 100   # save CSV every N completed comparisons

# =========================
# Functions
# =========================
def load_npz_embeddings(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    embeddings = []
    protein_ids = []
    for key in data.files:
        emb = data[key]
        if emb.ndim != 1:
            continue
        embeddings.append(emb)
        protein_ids.append(key)
    if len(embeddings) == 0:
        raise ValueError(f"No 1D embeddings found in {npz_path}")
    dim_set = {e.shape[0] for e in embeddings}
    if len(dim_set) != 1:
        raise ValueError(f"Inconsistent embedding dimensions in {npz_path}: {dim_set}")
    return np.vstack(embeddings), protein_ids

def greedy_one_to_one_threshold(sim_matrix, threshold):
    sim = sim_matrix.copy()
    matches = []
    while True:
        i, j = np.unravel_index(sim.argmax(), sim.shape)
        score = sim[i, j]
        if score < threshold:
            break
        matches.append((i, j, score))
        sim[i, :] = -1
        sim[:, j] = -1
    return matches

def compute_pair(pair):
    file_a, file_b = pair
    path_a = os.path.join(input_folder, file_a)
    path_b = os.path.join(input_folder, file_b)

    # --- Load embeddings ---
    embeddings_a, proteins_a = load_npz_embeddings(path_a)
    embeddings_b, proteins_b = load_npz_embeddings(path_b)

    # --- Cosine similarity matrix ---
    sim_matrix = cosine_similarity(embeddings_a, embeddings_b)

    # --- Greedy 1-to-1 matching above threshold ---
    matches = greedy_one_to_one_threshold(sim_matrix, similarity_threshold)
    avg_greedy = np.mean([m[2] for m in matches]) if matches else 0.0

    # --- Symmetric average of best hits ---
    avg_a_to_b = sim_matrix.max(axis=1).mean()
    avg_b_to_a = sim_matrix.max(axis=0).mean()
    avg_symmetric = (avg_a_to_b + avg_b_to_a) / 2

    # --- Reciprocal Best Hits (RBH) ---
    best_in_B_for_A = sim_matrix.argmax(axis=1)  # for each protein in A, best B
    best_in_A_for_B = sim_matrix.argmax(axis=0)  # for each protein in B, best A

    rbh_pairs = [(i, j) for i, j in enumerate(best_in_B_for_A)
                 if best_in_A_for_B[j] == i]

    num_rbh = len(rbh_pairs)
    avg_rbh_score = np.mean([sim_matrix[i, j] for i, j in rbh_pairs]) if rbh_pairs else 0.0

    # --- RBH fractions ---
    frac_rbh_a = num_rbh / embeddings_a.shape[0] if embeddings_a.shape[0] > 0 else 0
    frac_rbh_b = num_rbh / embeddings_b.shape[0] if embeddings_b.shape[0] > 0 else 0
    frac_rbh_sym = (frac_rbh_a + frac_rbh_b) / 2

    return {
        "Genome_A": file_a,
        "Genome_B": file_b,
        "Greedy_1to1_avg": avg_greedy,
        "Num_matches_above_threshold": len(matches),
        "Avg_max_A_to_B": avg_a_to_b,
        "Avg_max_B_to_A": avg_b_to_a,
        "Symmetric_avg": avg_symmetric,
        "Num_reciprocal_best_hits": num_rbh,
        "Avg_rbh_score": avg_rbh_score,
        "RBH_fraction_A": frac_rbh_a,
        "RBH_fraction_B": frac_rbh_b,
        "RBH_fraction_symmetric": frac_rbh_sym
    }


# =========================
# Load existing CSV if any
# =========================
if os.path.exists(output_csv):
    df_existing = pd.read_csv(output_csv)
    done_pairs = set(tuple(sorted([row["Genome_A"], row["Genome_B"]])) for _, row in df_existing.iterrows())
    results = df_existing.to_dict("records")
    print(f"Resuming: {len(done_pairs)} comparisons already done.")
else:
    done_pairs = set()
    results = []

# =========================
# Find NPZ files and remaining pairs
# =========================
npz_files = sorted([f for f in os.listdir(input_folder) if f.endswith(".npz")])
all_pairs = list(combinations(npz_files, 2))
remaining_pairs = [pair for pair in all_pairs if tuple(sorted(pair)) not in done_pairs]
print(f"{len(remaining_pairs)} comparisons remaining.")

# =========================
# Run in parallel with incremental saving
# =========================
completed_count = 0
with ProcessPoolExecutor(max_workers=num_threads) as executor:
    futures = {executor.submit(compute_pair, pair): pair for pair in remaining_pairs}
    for f in tqdm(as_completed(futures), total=len(futures), desc="Comparing NPZ pairs"):
        try:
            result = f.result()
            results.append(result)
            completed_count += 1

            # Save every N comparisons
            if completed_count % save_every == 0:
                pd.DataFrame(results).to_csv(output_csv, index=False)
        except Exception as e:
            pair = futures[f]
            print(f"Error processing pair {pair}: {e}")

# Final save
pd.DataFrame(results).to_csv(output_csv, index=False)
print(f"All comparisons saved to {output_csv}")
