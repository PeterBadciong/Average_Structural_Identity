# Determining Average Structural Identity

## 1. Obtaining Whole Genomes Using NCBI Taxonomy and Genome Browser

1. To quickly obtain genomes:  
   - Download a TSV file containing the genomes you want.  
   - Run `NCBI_Downloader.py` and use the TSV file as input.

## 2. Embedding with the ESM2 Model

1. **Requirements:**  
   - ESM2 embeddings are GPU-dependent, so you will need access to a GPU.  
   - Alternatively, you can use the Colab notebook: `Colab_ESM2_Embedder.ipynb`.  

2. **Model:**  
   - For the current edition of this project, use the **ESM2-T6-M8** model.

3. **Preparation:**  
   - Move the folder of `.npz` files to the directory where you will run step 3.

## 3. Pairwise Comparison

- Perform the pairwise comparison of embeddings to determine average structural identity.
