# Determining Average Structural Identity

## 1. Obtaining Whole Bacterial/Archaeal Genomes Using NCBI Taxonomy and Genome Browser

1. **To quickly obtain genomes:**  
   - Download a TSV file containing the genomes you want.  
   - Run `NCBIDownloader.py` and use the TSV file as input.

2. **Taxonomy of Genomes:**
   - Run `TaxonomyAdder.py` using the `.csv` file given from `NCBIDownloader.py`.  
   - This will give a taxonomy from domain to genus for every genome added.  
   - You will need to have **Ete3 installed** or run in a **conda environment**.
  
## 1.1 Obtaining Viral Genomes using ICTV taxonomy
   - Download the viral Genomes you wish to compare, and use the ICTV taxonomy file already provided

---

## 2. 3D Structural Predictions

1. **Requirements:**  
   - ESMFold is GPU dependant, access to a cluster is advised for large scale predictions.  
   - Alternatively, you can use the Colab notebook provided for single whole genome predictions.
   - A local machine version is currently in the works for batch prediction of whole genome protein structures

---

## 3. Running ASI_Compare

1. Using the predicted ESMFold structures and their original protein fasta files.
   - Inputs Required
     * -a/--aminoacids Folder containing the original protein fasta files
     * -s/--structures Folder containing the predicted protein structures
     * -o/--output Output directory
     * -b for bacterial/archaeal -v for viral
     * -t/--taxonomy Taxonomy.csv file
   - Optional Inputs
     * --threads (default 8)
     * --annotate true/false (default True)
     * --TM_Threshold (default 0.5)
     * --Paralog_TM_Threshold (default 0.8)
     * --sensitivity (default 9.5)





​​
