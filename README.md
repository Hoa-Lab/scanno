# SCanno Tutorial

### 1. Requirements:
- Python (v3.8.2)
- Scanpy (v1.81)
- Pandas (v1.3.4)
- Numpy (v1.20)
- Scikit_learn (v1.0.1)
- Tensorflow (v2.6.0)


### 2. Install
`pip install scanno`


### 3. Prepare Input Files
- Normalized (log1p) datasets in h5ad format (usually preprocessed using Scanpy)
- Marker gene information in csv format (the example file could be found in the example folder), which contains three columns:  
	1) "cell": This column contains expected cell names.  
	2) "marker": This column contains known marker genes for the given cell type. If there are more than one marker genes, insert "," between each gene symbol. Recommend using **two** markers for each cell type.
	3) "blocklist": This column contains the genes that used as markers for other cell types, but may also expressed in the given cell type.  
	

### 4. Prelabel Individual Dataset
- import packages and set variables  
```
import pandas  as pd
import scanpy as sc
from pathlib import Path
import scanno.scanno as sca
import tensorflow.keras as keras
f_marker='file_of_marker_gene_info'
fd_in='folder_of_input_h5ad_files'
fd_prelbl='folder_to_save_prelabeled_h5ad'
```

- load input files
```
df_marker=pd.read_csv(f_marker, index_col=0)
l_fname=list(Path(fd_in).glob('*.h5ad'))
```

- pre-label each individual dataset
```
for fname in l_fname:
	name=Path(fname).stem
	ada=sc.read(fname)
	ada=sca.prelabel(ada, df_marker)
	ada.write(f'{fd_prelbl}/{name}.h5ad')
```

### 5. Prepare Training Set for Annotation
- set variables
```
fd_train='folder_to_save_training_data'
```

- load prelabeled dataset
```
l_fname=list(Path(fd_prelbl).glob('*.h5ad'))
l_ada=[sc.read(i) for i in l_fname]
```

- prepare training data
```
df_train=sca.anno_traindata(l_ada)
df_train.to_csv(f'{fd_train}/data.csv')
```

- save features (genes)
```
Path(f'{fd_train}/gene.txt').write_text('\n'.join(df_train.columns.tolist()[0:-1]))
```
	*Note:
	Last columns is "prelabel", which should not be used as features.

### 6. Train the Annotation Model
- set variables
```
fd_mod='folder_to_save_model'
```

- load training data
```
df_train=pd.read_csv(f'{fd_train}/data.csv', index_col=0)
```

- Build Model
```
n_out=len(df_train['prelabel'].unique().tolist())
n_in=df_train.shape[1]-1
mod=sca.build_mod(n_in, n_out)
```
	*Note:  
	1) Customized model architecture could be used instead of calling "build_mod" function;  
	2) Default "build_mod" using dense model contains 2 hidden layers, with 15000 nodes each layer. These parameters could be adjusted by setting "layer" and "node" arguments.

- Train the Model
```
mod, lbl_enc=sca.train_model(df_train, mod, fd_mod)
```
	*Notes:  
	Default "train_model" using 2 epochs. This can be adjusted by setting "n_epoch" argument.

### 7. Annotation
- set variables
```
fd_anno='folder_to_save_annotated_dataset'
```

- load prelabeled dataset and gene list
```
l_fname=list(Path(fd_prelbl).glob('*.h5ad'))
l_gene=Path(f'{fd_train}/gene.txt').read_text().split('\n')
```

- Annotation Each Dataset
```
for fname in l_fname:
	name=Path(fname).stem
	ada=sc.read(fname)
	ada, df=sca.annotate(ada, mod, l_gene, lbl_enc)
	#save
	ada.write(f'{fd_anno}/{name}.h5ad')
	df.to_csv(f'{fd_anno}/{name}.csv')
```
	*Note:  
	1) Function "annotate" returns the annotated anndata (with cell labeling in the 'pred' column), and the dataframe contains probablities for all cell types;
	2) By default, "annotate" function set cells with all probablities less than 0.7 as "Unknown". This threshold can be ajusted by set "t" argument.
	
### 8. Identify Marker Genes
- set variables
```
cell='cell_of_interest'
fd_mg='folder_to_save_result'
```

- load annotated anndata and probablity dataframe
```
l_fname=list(Path(fd_anno).glob('*.h5ad'))
l_ada=[sc.read(i) for i in l_fname]
l_fname=list(Path(fd_anno).glob('*.csv'))
l_df=[pd.read_csv(i, index_col=0) for i in l_fname]
```

- Identify Marker Genes
```
df_mg=sca.get_marker(cell, l_ada, l_df)
df_mg.to_csv(f'{fd_mg}/{cell}.csv')
```

### 9. Prepare Inputs for Trajectory Analysis
- set variables
```
fd_cell='folder_to_save_combined_anndata'
```

- extract cells from each dataset, and combined into one anndata
```
l_fname=list(Path(fd_anno).glob('*.h5ad'))
l_ada=[]
for fname in l_fname:
	ada=sc.read(fname)
	ada=ada[ada.obs['pred']==cell, :].copy()
	l_ada.append(ada)
ada=ad.concat(l_ada)
ada.write(f'{fd_cell}/{cell}.h5ad')
```

### 10. Prepare Training Data for Trajectory Analysis
- set variables
```
fd_trajtrain='folder_to_save_trajectory_training_data'
col='column_of_true_timestamp'
begin='start_group'
end='end_group'
```
	*Note:
	1) "col" is the column name in ada.obs, which contains the info of start/end cells
	2) "begin"/"end" is the label value in "col", which indicates the start/end cells.

- prepare training data
```
df_trajtrain=sca.traj_traindata(ada, col, begin, end)
df_trajtrain.to_csv(f'{fd_trajtrain}/data.csv')
Path(f'{fd_trajtrain}/gene.txt').write_text('\n'.join(df_trajtrain.columns.tolist()[0:-1]))
```

### 11. Train the Trajectory Model
- set variables
```
fd_trajmod='folder_to_save_trajectory_model'
```

- build model
```
n_in=df_trajtrain.shape[1]-1
mod=sca.build_mod(n_in, n_out=2)
```

- train the model
```
mod, lbl_enc=sca.train_model(df_trajtrain, mod, fd_trajmod, col_lbl='true_label')
```

### 12. Calculate Trajectory
- set variables
```
fd_trajgene='folder_to_save_trajectory_genes'
```

- load inputs
```
ada=sc.read(f'{fd_cell}/{cell}.h5ad')
l_gene=Path(f'{fd_trajtrain}/gene.txt').read_text().split('\n')
```

- calculate probability score (pseudotime)
```
ada=sca.pseudotime(ada, mod, l_gene, lbl_enc, end, raw=False)
```
	*Note: scores are saved in the "prob_score" column in ada.obs

- identify up/down regulated genes during trajectory
```
df_trajgene=sca.traj_gene(ada, raw=False)
df_trajgene.to_csv(f'{fd_trajgene}/{cell}.csv')
```


