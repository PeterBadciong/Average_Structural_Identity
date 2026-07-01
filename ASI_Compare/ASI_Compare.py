import os
import subprocess
from pathlib import Path
import sys
import argparse

# Command-line argument parsing

def parse_args():
    parser = argparse.ArgumentParser(
        description="ASI Pipeline — Foldseek RBH, AAI, Annotation, TM/AAI Merge"
    )

    # REQUIRED INPUTS
    parser.add_argument(
        "-a", "--aminoacids",
        required=True,
        help="Folder containing original protein FASTA files"
    )
    parser.add_argument(
        "-s", "--structures",
        required=True,
        help="Folder containing predicted protein structures (PDB files)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory"
    )

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "-b", "--bacterial",
        action="store_true",
        help="Run pipeline in Bacterial/Archaeal mode"
    )
    mode_group.add_argument(
        "-v", "--viral",
        action="store_true",
        help="Run pipeline in Viral mode"
    )

    parser.add_argument(
        "-t", "--taxonomy",
        required=True,
        help="Taxonomy CSV file"
    )

    # OPTIONAL INPUTS
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--annotate", type=str, choices=["true", "false"], default="true")
    parser.add_argument("--TM_Threshold", type=float, default=0.5)
    parser.add_argument("--Paralog_TM_Threshold", type=float, default=0.8)
    parser.add_argument("--sensitivity", type=float, default=9.5)

    return parser.parse_args()



# Utility functions


def run(cmd):
    print(f"\n[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def check_exists(path):
    if not os.path.exists(path):
        print(f"[ERROR] Required file not found: {path}")
        sys.exit(1)



# MAIN PIPELINE


def main():

    args = parse_args()

    # Required inputs
    FAA_FOLDER = args.aminoacids
    INPUT_DIR = args.structures
    WORKDIR = args.output
    ACTIVE_TAXONOMY = args.taxonomy

    # Mode selection
    MODE = "viral" if args.viral else "bacterial"

    # Optional inputs
    THREADS = args.threads
    RUN_ANNOTATION = (args.annotate.lower() == "true")
    TM_THRESHOLD = args.TM_Threshold
    PARALOG_TM_THRESHOLD = args.Paralog_TM_Threshold
    SENSITIVITY = args.sensitivity

    # Step toggles (same as your original script)
    RUN_STEP1 = True
    RUN_STEP1B = True
    RUN_STEP2 = True
    RUN_STEP3 = True
    RUN_STEP4 = True
    RUN_STEP5 = True
    RUN_STEP6 = True
    RUN_STEP7 = True
    RUN_PLDDT = True

    # Script paths (unchanged)
    SCRIPT1 = "ASI_Compare_Scripts/ASI_Foldseek.py"
    SCRIPT1B = "ASI_Compare_Scripts/ASI_Paralogs.py"
    SCRIPT2 = "ASI_Compare_Scripts/ASI_Annotate.py"
    SCRIPT3 = "ASI_Compare_Scripts/ASI_Quality.py"
    SCRIPT4 = "ASI_Compare_Scripts/ASI_AAI.py"
    SCRIPT5 = "ASI_Compare_Scripts/ASI_AAI_Compare.py"
    SCRIPT6 = "ASI_Compare_Scripts/ASI_Merged_Bacterial.py"
    SCRIPT6B = "ASI_Compare_Scripts/ASI_Merged_Viral.py"
    SCRIPT7 = "ASI_Compare_Scripts/ASI_Boxplot_Bacterial.py"
    SCRIPT7B = "ASI_Compare_Scripts/ASI_Boxplot_Viral.py"

    # Workdir subfolders
    WORKDIR_1 = f"{WORKDIR}/foldseek"
    WORKDIR_1B = f"{WORKDIR}/paralogs"
    WORKDIR_2 = f"{WORKDIR}/sequence"
    WORKDIR_3 = f"{WORKDIR}/annotation_comparison"
    WORKDIR_4 = f"{WORKDIR}/pairwise_rbh_aai"
    WORKDIR_5 = f"{WORKDIR}/foldseek_vs_aai_comparison"

    # AAI parameters
    DIAMOND = "/usr/bin/diamond"
    HMMSEARCH = "/usr/bin/hmmsearch"
    EVALUE = 1e-5
    AAI_CUTOFF = 30.0
    COVERAGE_CUTOFF = 0.50

    # Mode-specific HMMs
    if MODE == "viral":
        ACTIVE_MERGED_HMM = "HMMs/PhrogAndVog.hmm"
        ACTIVE_MERGED_ANNOT = "HMMs/MergedPhrogVog.tsv"
    else:
        ACTIVE_MERGED_HMM = "HMMs/profiles.hmm"
        ACTIVE_MERGED_ANNOT = "HMMs/KEGG_Annotations.tsv"

    # Ensure workdirs exist
    for d in [WORKDIR_1, WORKDIR_1B, WORKDIR_2, WORKDIR_3, WORKDIR_4, WORKDIR_5]:
        Path(d).mkdir(parents=True, exist_ok=True)

    
    # STEP 1 — Foldseek RBH
   
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

        check_exists(tm_tsv)
        check_exists(all_tsv)

    else:
        tm_tsv = f"{WORKDIR_1}/pairwise_RBH_matches_TMfiltered.tsv"
        check_exists(tm_tsv)

    
    # STEP 1B — Paralog detection
   
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

    
    # STEP 2 — AAI Pipeline
   
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

    final_csv = f"{WORKDIR_2}/ALL_GENOMES_FINAL.csv"
    check_exists(final_csv)

    # STEP 3 — Annotation Comparison
    
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

    
    # STEP 4 — Pairwise RBH AAI

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

    
    # STEP 5 — Foldseek vs AAI RBH
   
    if RUN_STEP5:
        print("\n==============================")
        print(" STEP 5: Foldseek vs AAI RBH  ")
        print("==============================")

        tm_tsv = f"{WORKDIR_1}/pairwise_RBH_matches_TMfiltered.tsv"
        aai_csv = f"{WORKDIR_4}/all_rbh_proteins.csv"
        pdb_dir = f"{WORKDIR_1}/split"
        map_file = f"{WORKDIR_1}/protein_map.tsv"

        for f in [tm_tsv, aai_csv, pdb_dir, map_file]:
            check_exists(f)

        cmd5 = (
            f"python3 {SCRIPT5} "
            f"--tm-tsv {tm_tsv} "
            f"--tm-bidir {WORKDIR_1}/pairwise_TM_bidirectional.tsv "
            f"--aai-csv {aai_csv} "
            f"--pdb-dir {pdb_dir} "
            f"--map-file {map_file} "
            f"--outdir {WORKDIR_5} "
            f"{'--compute-plddt' if RUN_PLDDT else ''}"
        )
        run(cmd5)

    
    # STEP 6 — Merge TM RBH + AAI
   
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

        
        # STEP 7 — Figures
       
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

    print("\n==============================")
    print(" FULL PIPELINE COMPLETE ")
    print("==============================")


if __name__ == "__main__":
    main()
