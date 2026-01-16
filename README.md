Determining Average Structural Identity
1. Obtaining whole genomes using NCBI Taxonomy and Genome browser
1a. To obtain these genomes quick, download a tsv file from the genomes you want, and run NCBI_Downloader.py, use the tsv file as your input
2. Embedding with the ESM2 Model
2a. ESM2 Embeddings are GPU-Dependant, so you will need access to a GPU, alternatively, you can use the colab notebook Colab_ESM2_Embedder.ipynb , for the current edition of this project, use the ESM2-T6-M8 model
2b. Move the folder of .npz files to the directory you wish to run step 3 in
3. Pairwise Comparison
