# Determining Average Structural Identity

## 1. Obtaining Whole Genomes Using NCBI Taxonomy and Genome Browser

1. **To quickly obtain genomes:**  
   - Download a TSV file containing the genomes you want.  
   - Run `NCBIDownloader.py` and use the TSV file as input.

2. **Taxonomy of Genomes:**
   - Run `TaxonomyAdder.py` using the `.csv` file given from `NCBIDownloader.py`.  
   - This will give a taxonomy from domain to genus for every genome added.  
   - You will need to have **Ete3 installed** or run in a **conda environment**.

---

## 2. Embedding with the ESM2 Model

1. **Requirements:**  
   - ESM2 embeddings are GPU-dependent, so you will need access to a GPU.  
   - Alternatively, you can use the Colab notebook: `ColabESM2Embedder.ipynb`.

2. **Model:**  
   - For the current edition of this project, use the **ESM2-T6-M8** model.

---

## 3. Running a Pairwise Comparison

1. Using the ESM2 Embedded files, run `ESM2Compare.py`.
   - Select count for threads used; depending on the number of genomes, the time will increase greatly under the formula:

$$
\text{Total Comparisons} = \frac{n \cdot (n - 1)}{2}
$$

2. Run `StatisticalTaxonomyTest.py` using your pairwise comparisons from `ESM2Compare.py` and your taxonomy file from `TaxonomyAdder.py`.

3.  This will provide an output based on different scores
      - Greedy-1to1-avg
      - symmetric-avg
      - Avg-rbh-score
      - rbh-fraction-symmetric
  
4.  This values provide the final score using the formula
   
$$
\text{Average score} = \frac{(\texttt{avg-rbh-score} \times \texttt{rbh-fraction-symmetric}) + \texttt{symmetric-avg} + \texttt{greedy-1to1-avg}}{3}
$$





​​
