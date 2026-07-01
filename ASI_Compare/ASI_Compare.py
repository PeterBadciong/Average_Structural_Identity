import os
import subprocess
from pathlib import Path
import sys

# ============================================================
# Choose Bacterial/Archaeal or Viral
# ============================================================
MODE = "Bacterial"   # options: "viral" or "bacterial"

# ============================================================
# STEP TOGGLES — TURN STEPS ON/OFF HERE
# ============================================================

RUN_STEP1 = True   # Foldseek RBH
RUN_STEP1B = True      # Paralog Identification
RUN_STEP2 = True       # Annotation Pipeline
RUN_STEP3 = True       # Annotation Comparison
RUN_STEP4 = True       # Pairwise RBH AAI
RUN_STEP5 = True       # Foldseek vs AAI comparison (Part 5)
RUN_STEP6 = True       # Make merged File with AAI and ASI
RUN_STEP7 = True
RUN_PLDDT = True
RUN_ANNOTATION = True
# ============================================================
# USER CONFIGURATION (ALL SETTINGS LIVE HERE)
# ============================================================

# ---- Script paths ----
SCRIPT1 = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Foldseek.py"
SCRIPT1B = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Paralogs.py"
SCRIPT2 = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Annotate.py"
SCRIPT3 = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Quality.py"
SCRIPT4 = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_AAI.py"
SCRIPT5 = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_AAI_Compare.py"
SCRIPT6 = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Merged_Bacterial.py"
SCRIPT6B = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Merged_Viral.py"
SCRIPT7 = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Boxplot_Bacterial.py"
SCRIPT7B = "/storage2/scratch/pbadciong/ASI/FoldseekerSuite/FoldseekerAllinOne/ASI_Boxplot_Viral.py"

# ---- Unified parameters ----
WORKDIR = "/storage2/scratch/pbadciong/ASI/TestingFolder/376_Bacteria"
THREADS = 64

# Subdirectories automatically created
WORKDIR_1 = f"{WORKDIR}/foldseek"
WORKDIR_1B = f"{WORKDIR}/paralogs"
WORKDIR_2 = f"{WORKDIR}/sequence"
WORKDIR_3 = f"{WORKDIR}/annotation_comparison"
WORKDIR_4 = f"{WORKDIR}/pairwise_rbh_aai"
WORKDIR_5 = f"{WORKDIR}/foldseek_vs_aai_comparison"

# ---- Script 1 parameters ----
INPUT_DIR = "/storage2/scratch/pbadciong/Average_Structural_ID/ESMFoldPredictedPDBFiles"
SENSITIVITY = "9.5"
TM_THRESHOLD = 0.50
PARALOG_TM_THRESHOLD = 0.8

# ---- Script 2 parameters ----
FAA_FOLDER = "/storage2/scratch/pbadciong/Average_Structural_ID/ASI_Prodigal_Proteins"

DIAMOND = "/usr/bin/diamond"
HMMSEARCH = "/usr/bin/hmmsearch"
# Viral vs Bacterial taxonomy files
VIRAL_TAXONOMY = "/storage2/scratch/pbadciong/Average_Structural_ID/ICTV_Taxonomy.csv"
BACTERIAL_TAXONOMY = "/storage2/scratch/pbadciong/Average_Structural_ID/Automated_Bacteria_Taxonomy.csv"
# Select taxonomy file based on MODE
if MODE == "viral":
    ACTIVE_TAXONOMY = VIRAL_TAXONOMY
else:
    ACTIVE_TAXONOMY = BACTERIAL_TAXONOMY


# Viral vs Bacterial HMM + Annotation files
VIRAL_MERGED_HMM = "/storage2/scratch/pbadciong/Average_Structural_ID/PhrogAndVog/PhrogAndVog.hmm"
VIRAL_MERGED_ANNOT = "/storage2/scratch/pbadciong/Average_Structural_ID/PhrogAndVog/MergedPhrogVog.tsv"

BACTERIAL_MERGED_HMM = "/storage2/scratch/pbadciong/Average_Structural_ID/BacterialHMMs/profiles.hmm"
BACTERIAL_MERGED_ANNOT = "/storage2/scratch/pbadciong/Average_Structural_ID/BacterialHMMs/KEGG_Annotations.tsv"

# Select HMM + annotation files based on MODE
if MODE == "viral":
    ACTIVE_MERGED_HMM = VIRAL_MERGED_HMM
    ACTIVE_MERGED_ANNOT = VIRAL_MERGED_ANNOT
else:
    ACTIVE_MERGED_HMM = BACTERIAL_MERGED_HMM
    ACTIVE_MERGED_ANNOT = BACTERIAL_MERGED_ANNOT


EVALUE = 1e-5

#Part 4 AAI_cutoff
AAI_CUTOFF = 30.0
COVERAGE_CUTOFF = 0.50



# ============================================================
# UTILITIES
# ============================================================

def run(cmd):
    print(f"\n[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def check_exists(path):
    if not os.path.exists(path):
        print(f"[ERROR] Required file not found: {path}")
        sys.exit(1)

# ============================================================
# MAIN PIPELINE
# ============================================================

def main():

    # Ensure workdirs exist
    Path(WORKDIR_1).mkdir(parents=True, exist_ok=True)
    Path(WORKDIR_1B).mkdir(parents=True, exist_ok=True)
    Path(WORKDIR_2).mkdir(parents=True, exist_ok=True)
    Path(WORKDIR_3).mkdir(parents=True, exist_ok=True)
    Path(WORKDIR_4).mkdir(parents=True, exist_ok=True)
    Path(WORKDIR_5).mkdir(parents=True, exist_ok=True)


    # --------------------------------------------------------
    # STEP 1
    # --------------------------------------------------------
    if RUN_STEP1:
        print("\n==============================")
        print(" STEP 1: Running Foldseek RBH ")
        print("==============================")

        cmd1 = (
            f"python3 {SCRIPT1} "
            f"--input-dir {INPUT_DIR} "
            f"--workdir {WORKDIR_1} "
            f"--threads {THREADS} "
            f"--sensitivity {SENSITIVITY} "
            f"--tm-threshold {TM_THRESHOLD}"
        )
        run(cmd1)

        tm_tsv = f"{WORKDIR_1}/pairwise_RBH_matches_TMfiltered.tsv"
        all_tsv = f"{WORKDIR_1}/pairwise_RBH_matches_ALL.tsv"

        print("\nChecking Script 1 outputs...")
        check_exists(tm_tsv)
        check_exists(all_tsv)

        print("✓ Found TM-filtered RBH file")
        print("✓ Found ALL RBH file")

    else:
        tm_tsv = f"{WORKDIR_1}/pairwise_RBH_matches_TMfiltered.tsv"
        check_exists(tm_tsv)
        print("\n[SKIP] Step 1 disabled — using existing RBH file")


        # --------------------------------------------------------
        # STEP 1B — Paralog Detection (within-genome)
        # --------------------------------------------------------
        if RUN_STEP1B:
            print("\n==============================")
            print(" STEP 1B: Detecting Paralogs  ")
            print("==============================")

            cmd1b = (
                f"python3 {SCRIPT1B} "
                f"--input-dir {INPUT_DIR} "
                f"--workdir {WORKDIR_1B} "
                f"--threads {THREADS} "
                f"--sensitivity {SENSITIVITY} "
                f"--tm-threshold {PARALOG_TM_THRESHOLD}"
            )

            run(cmd1b)

            paralog_tsv = f"{WORKDIR_1B}/paralogs_TM0.8.tsv"
            check_exists(paralog_tsv)
            print("✓ Paralog file generated:", paralog_tsv)

        else:
            print("\n[SKIP] Step 1B disabled")


    # --------------------------------------------------------
    # STEP 2
    # --------------------------------------------------------
    if RUN_STEP2:
        print("\n==============================")
        print(" STEP 2: Running AAI Pipeline ")
        print("==============================")


        cmd2 = (
            f"python3 {SCRIPT2} "
            f"--input-tsv {tm_tsv} "
            f"--faa-folder {FAA_FOLDER} "
            f"--diamond {DIAMOND} "
            f"--hmmsearch {HMMSEARCH} "
            f"--merged-hmm {ACTIVE_MERGED_HMM} "
            f"--merged-annot {ACTIVE_MERGED_ANNOT} "
            f"--threads {THREADS} "
            f"--evalue {EVALUE} "
            f"--workdir {WORKDIR_2} "
        )

        if RUN_ANNOTATION:
            cmd2 += " --run-annotation"

        run(cmd2)

    else:
        print("\n[SKIP] Step 2 disabled — using existing AAI output")

    final_csv = f"{WORKDIR_2}/ALL_GENOMES_FINAL.csv"
    check_exists(final_csv)

    # --------------------------------------------------------
    # STEP 3
    # --------------------------------------------------------
    if RUN_STEP3:
        print("\n==============================")
        print(" STEP 3: Annotation Comparison ")
        print("==============================")

        cmd3 = (
            f"python3 {SCRIPT3} "
            f"--input-csv {final_csv} "
            f"--output-dir {WORKDIR_3} "
            f"--threads {THREADS}"
        )
        run(cmd3)

    else:
        print("\n[SKIP] Step 3 disabled")

    # --------------------------------------------------------
    # STEP 4
    # --------------------------------------------------------
    if RUN_STEP4:
        print("\n==============================")
        print(" STEP 4: Pairwise RBH AAI     ")
        print("==============================")

        cmd4 = (
            f"python3 {SCRIPT4} "
            f"--faa-folder {FAA_FOLDER} "
            f"--outdir {WORKDIR_4} "
            f"--threads {THREADS} "
            f"--aai-cutoff {AAI_CUTOFF} "
            f"--cov-cutoff {COVERAGE_CUTOFF}"
        )
        run(cmd4)

    else:
        print("\n[SKIP] Step 4 disabled")

    # --------------------------------------------------------
    # STEP 5
    # --------------------------------------------------------
    if RUN_STEP5:
        print("\n==============================")
        print(" STEP 5: Foldseek vs AAI RBH  ")
        print("==============================")

        tm_tsv = f"{WORKDIR_1}/pairwise_RBH_matches_TMfiltered.tsv"
        aai_csv = f"{WORKDIR_4}/all_rbh_proteins.csv"
        pdb_dir = f"{WORKDIR_1}/split"

        check_exists(tm_tsv)
        check_exists(aai_csv)
        check_exists(pdb_dir)

        map_file = f"{WORKDIR_1}/protein_map.tsv"
        check_exists(map_file)

        cmd5 = (
            f"python3 {SCRIPT5} "
            f"--tm-tsv {WORKDIR_1}/pairwise_RBH_matches_TMfiltered.tsv "
            f"--tm-bidir {WORKDIR_1}/pairwise_TM_bidirectional.tsv "
            f"--aai-csv {WORKDIR_4}/all_rbh_proteins.csv "
            f"--pdb-dir {WORKDIR_1}/split "
            f"--map-file {WORKDIR_1}/protein_map.tsv "
            f"--outdir {WORKDIR_5} "
            f"{'--compute-plddt' if RUN_PLDDT else ''}"
        )


        run(cmd5)

    if RUN_STEP6:
        print("\n==============================")
        print(" STEP 6: Merge TM RBH with AAI")
        print("==============================")

        script6_to_run = SCRIPT6B if MODE == "viral" else SCRIPT6

        cmd6 = (
            f"python3 {script6_to_run} "
            f"--tm-csv {WORKDIR_1}/proteome_RBH_scored_TMfiltered.csv "
            f"--aai-summary {WORKDIR_4}/pairwise_summary.csv "
            f"--taxonomy {ACTIVE_TAXONOMY} "
            f"--outdir {WORKDIR}/tm_with_aai"
        )

        run(cmd6)

        # --------------------------------------------------------
        # STEP 7
        # --------------------------------------------------------
        if RUN_STEP7:
            print("\n==============================")
            print(" STEP 7: Taxonomic TM vs AAI Figures ")
            print("==============================")

            merged_input = f"{WORKDIR}/tm_with_aai/TM_AAI_Taxonomy_Merged.csv"
            check_exists(merged_input)

            fig_outdir = f"{WORKDIR}/figures"
            Path(fig_outdir).mkdir(parents=True, exist_ok=True)

            script7_to_run = SCRIPT7B if MODE == "viral" else SCRIPT7

            cmd7 = (
                f"python3 {script7_to_run} "
                f"--input {merged_input} "
                f"--outdir {fig_outdir}"
        )

            run(cmd7)

        else:
            print("\n[SKIP] Step 7 disabled")



    else:
        print("\n[SKIP] Step 5 disabled")


    print("\n==============================")
    print(" FULL PIPELINE COMPLETE ")
    print("==============================")

if __name__ == "__main__":
    main()
