import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import zscore
import tensorflow.keras as keras
from sklearn import preprocessing
import pickle
from collections import Counter

#------------------------------------------------------
def prelabel(ada, df_meta, raw=True, z=1.2):
    #load to df
    if raw:
        ada_tmp=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
        df=ada_tmp.to_df()
    else:
        df=ada.to_df()
    #clean meta
    df_meta=df_meta.fillna('')
    l_all=df_meta['marker'].tolist()
    l_all=[i.strip().split(',') for i in l_all]
    l_all=sum(l_all, [])
    l_all=[i.strip() for i in l_all]
    l_all=[i for i in l_all if i in df.columns]
    #zscore
    df_tmp=df.loc[:, l_all].copy()
    df_tmp=df_tmp.apply(zscore)
    #label
    ada.obs['prelabel']='TBD'
    for cell in df_meta.index:
        #get genes
        l_marker=df_meta['marker'][cell].strip().split(',')
        l_marker=[i.strip() for i in l_marker]
        l_block=df_meta['blocklist'][cell].strip().split(',')
        l_block=[i.strip() for i in l_block]
        #check marker
        if len([i for i in l_marker if i not in df.columns])>0:
            print(f'{cell} markers not expressed...')
            continue
        #label
        l_other=[i for i in l_all if i not in l_marker+l_block]
        dfi=df_tmp.copy()
        for gene in l_marker:
            dfi=dfi.loc[dfi[gene]>z, :]
        for gene in l_other:
            dfi=dfi.loc[dfi[gene]<z, :]
        if dfi.shape[0]>0:
            ada.obs.loc[dfi.index, ['prelabel']]=cell
    return ada


#-----------------------------------------------------------------
def anno_traindata(l_ada, raw=True, n_min=5):
    #concat
    l_df=[]
    for ada in l_ada:
        if raw:
            ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
        ada=ada[ada.obs['prelabel']!='TBD', :].copy()
        df=ada.to_df()
        df['prelabel']=ada.obs['prelabel']
        l_df.append(df)
    df=pd.concat(l_df).dropna(axis=1)
    #filter rare cells
    d_cnt=Counter(df['prelabel'].tolist())
    l_cell=[i for i in df['prelabel'].unique().tolist() if d_cnt[i]>n_min]
    df=df.loc[df['prelabel'].isin(l_cell), :]
    return df


#--------------------------------------------------------------------
def build_mod(n_in, n_out, layer=2, node=15000, r=0.5, lr=0.0001):
 	#build model
	model=keras.models.Sequential()
	model.add(keras.layers.Input(n_in,))
	for i in range(layer):
		model.add(keras.layers.Dropout(r))
		model.add(keras.layers.Dense(node, activation='swish', kernel_initializer='he_normal'))
	model.add(keras.layers.Dense(n_out, activation='softmax'))
	print(model.summary())
    #compile
	opt=keras.optimizers.Adam(learning_rate=lr)
	model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['categorical_accuracy'])
	return model


def train_model(df, model, fd_out, n_epoch=2, batch=32, vs=0.1, col_lbl='prelabel'):
    #input
    le=preprocessing.LabelEncoder()
    y=le.fit_transform(df[col_lbl].values)
    y=keras.utils.to_categorical(y)
    X=df.drop(col_lbl, axis=1).values
    #save le
    with open(f'{fd_out}/label_encode.pkl', 'wb') as f:
        pickle.dump(le, f)
    #train
    cb_checkpoint=keras.callbacks.ModelCheckpoint(f'{fd_out}/model', save_best_only=True)
    cb_tensorboard=keras.callbacks.TensorBoard(f'{fd_out}/log')
    history=model.fit(X, y, epochs=n_epoch, batch_size=batch, callbacks=[cb_checkpoint, cb_tensorboard], validation_split=vs)
    return model, le
   
#-----------------------------------------------------------------------
def annotate(ada, model, l_gene, le, raw=True, t=0.7):
    #get df
    if raw:
        ada_tmp=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
        df=ada_tmp.to_df()
    else:
        df=ada.to_df()
    df=df.loc[:, l_gene].fillna(0)
    #pred
    df_pred=pd.DataFrame(model.predict(df.values), index=df.index)
    df_pred.columns=le.inverse_transform(df_pred.columns.tolist())
    #label
    df_pred['pred']=df_pred.idxmax(axis=1)
    df_pred.loc[df_pred.max(axis=1)<t, ['pred']]='Unknown'
    #add anno to adata
    ada.obs=ada.obs.merge(df_pred.loc[:, ['pred']], left_index=True, right_index=True)
    ada.obs.loc[ada.obs['prelabel']!='TBD', ['pred']]=ada.obs['prelabel']
    return ada, df_pred

#-----------------------------------------------------------------------
def get_marker(cell, l_ada, l_df, raw=True):
    #concat df
    df=pd.concat(l_df)
    #concat ada
    l_cnt=[]
    for ada in l_ada:
        if raw:
            ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
        dfi=ada.to_df()
        l_cnt.append(dfi)
    df_cnt=pd.concat(l_cnt).fillna(0)
    #corr
    s=df_cnt.corrwith(df[cell])
    df_corr=pd.DataFrame(s, columns=['score']).dropna()
    df_corr=df_corr.sort_values('score', ascending=False)
    return df_corr

#------------------------------------------------------------------------
def traj_traindata(ada, col, begin, end, raw=False):
    if raw:
        ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
    df=ada.to_df()
    df['true_label']=ada.obs[col]
    #filter
    df=df.loc[df['true_label'].isin([begin, end]), :].copy()
    return df


#-------------------------------------------------------------------------
def pseudotime(ada, model, l_gene, le, end, raw=True):
    #load
    if raw:
        ada_tmp=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
        df=ada_tmp.to_df()
    else:
        df=ada.to_df()
    df=df.loc[:, l_gene].fillna(0)
    #predict
    df_pred=pd.DataFrame(model.predict(df.values), index=df.index)
    df_pred.columns=le.inverse_transform(df_pred.columns.tolist())
    #clean
    df_pred=df_pred.loc[:, [end]]
    df_pred.columns=['prob_score']
    ada.obs=ada.obs.merge(df_pred, left_index=True, right_index=True)
    return ada

#-------------------------------------------------------------------------
def traj_gene(ada, col='prob_score', raw=True):
    #load
    if raw:
        ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
    df=ada.to_df()
    #corr
    s=df.corrwith(ada.obs[col])
    df_corr=pd.DataFrame(s, columns=['score']).dropna()
    df_corr=df_corr.sort_values('score', ascending=False)
    return df_corr
