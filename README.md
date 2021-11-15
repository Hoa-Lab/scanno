# SCanno Tutorial

### 1. Requirements:
- Python (v3.8.2)
- Scanpy (v1.81)
- Pandas (V1.3.4)
- Scikit_learn (v1.0.1)
- Tensorflow (v2.6.0)


### 2. Install
`pip install scanno`


### 3. Prepare Input Files
- Normalized (log1p) datasets in h5ad format (usually preprocessed using Scanpy)
- Marker gene information in csv format (the example file could be found in the example folder), which contains three columns:  
	1) "cell": This column contains expected cell names.  
	2) "marker": This column contains known marker genes for the given cell type. If there are more than one marker genes, insert "," between each gene symbol. Recomend using **two** markers for each cell type.
	3) "blocklist": This column contains the genes that used as markers for other cell types, but may also expressed in the given cell type.  
	

### 4. Annotation
