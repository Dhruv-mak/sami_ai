import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import re

import logging
logger = logging.getLogger(__name__)

def csv2h5ad(compound_type, df, data_saving_path, split=True):
    """
    Convert a DataFrame containing omics data into .h5ad format for use in Scanpy.

    This function processes omics data by:
    - Filtering out rows where all feature values are zero.
    - Extracting spatial coordinates (x, y).
    - Assigning metadata including regions and omics type.
    - Optionally, splitting and saving separate .h5ad files for each unique region.

    Parameters
    ----------
    compound_type : str
        The type of omics data (e.g., 'metabolomics', 'lipidomics', etc.).
    df : pandas.DataFrame
        DataFrame containing the omics data.
    data_saving_path : str
        Path to the directory where the .h5ad files will be saved.
    split : bool, optional
        If True, saves separate .h5ad files for each unique region. Default is True.

    Returns
    -------
    None
        Writes .h5ad files to the specified directory and logs the processing status.
    """

    df =df.rename(columns={'tissue_id':'region'})
    df['region'] = df['region'].astype(str)
    feat_cols = df.columns[3:]
    df = df[df[feat_cols].sum(axis=1)!=0]
    nodefeats = df[feat_cols]
    nodefeats = np.array(nodefeats)
    df['obsindex'] = df['x'].astype(str)+"x"+df['y'].astype(str)
    xs=df["x"]
    ys=df["y"]
    xs=np.array(xs)
    ys=np.array(ys)
    pos=np.stack((xs,ys), axis=-1)
    adata=sc.AnnData(X=nodefeats,dtype='float32')
    obs=df.set_index('obsindex')[['region']]
    adata.obs=obs
    adata.obsm['spatial'] = pos
    var = pd.DataFrame({'omics':compound_type},index=df[feat_cols].columns)
    adata.var = var
    adata.raw = adata
    
    # adata.write(os.path.join(data_saving_path,f'{split_file[0]}_{split_file[1]}.h5ad'))
    if split==True:
        unique_id=df['region'].unique().tolist()
        for i in range(len(unique_id)):
            adata_sub=adata[adata.obs['region']==unique_id[i]].copy()
            adata_sub.write(os.path.join(data_saving_path,f'{unique_id[i]}_{compound_type}.h5ad'))
    logger.info(f'csv2h5ad completed')


                    
#merge omics data by coordinates and return .h5ad file
def pooldata(dict_compound_type, data_reading_path, data_saving_path, split=True):
    error_message  =""
    adata_dict={}
    for file, compound_type in dict_compound_type.items():
        logger.info(f'pool step processing {file}')          
        # print("converting", os.path.join(data_path, file))
        data_temp = pd.read_parquet(os.path.join(data_reading_path, file))
        data_temp =data_temp.rename(columns={'tissue_id':'region'})
        data_temp['region'] = data_temp['region'].astype(str)
        feat_cols = data_temp.columns[3:]
        data_temp = data_temp[data_temp[feat_cols].sum(axis=1)!=0]
        nodefeats = data_temp[feat_cols]
        nodefeats = np.array(nodefeats)
        data_temp['obsindex'] = data_temp['x'].astype(str)+"x"+data_temp['y'].astype(str)

        adata_temp=sc.AnnData(X=nodefeats,dtype='float32')
        obsindex = pd.DataFrame(index=data_temp['obsindex'])
        adata_temp.obs = obsindex
        varindex = pd.DataFrame(index=data_temp[feat_cols].columns)
        adata_temp.var = varindex
        adata_temp.raw = adata_temp
        # adding all the lipids, metabolites, glycans to a dictionary
        adata_dict[compound_type]=adata_temp

    #remove duplicate columns between metabolomics and lipidomics from metabolomics
    # common_elements = list(set(adata_dict['metabolomics'].var.index.to_list()).intersection(adata_dict['lipidomics'].var.index.to_list()))
    # adata_dict['metabolomics']=adata_dict['metabolomics'][:,~adata_dict['metabolomics'].var.index.isin(common_elements)]

    # adata_dict
    # Pickle the adata_dict
    # with open(os.path.join(data_saving_path, r'C:\Users\ghari\Documents\OPS\SAMI\adata_dict.pkl'), 'wb') as f:
    #     pickle.dump(adata_dict, f)

    # manually overriding pixel mis alignment
    # adata_dict['lipids'].obs = adata_dict['sm'].obs 

    # combining all the regions of same type of omics data
    adata=ad.concat(adata_dict,join='inner',axis=1,label='omics')

    # the concatenation usually happens on the basis of the index, so if the index is not same then the data will not be concatenated
    # so we need to check if the data is concatenated or not
    # if its empty and the pixes are same, then we can do a blind merge and give a warning to user
    same_pixles = True
    if len(adata.obs.index) == 0:
        error_message = "No common pixels found in the data, So performing BLIND merge, please check the data post clustering"
        logger.warning("mismatch, no common pixels found in the data")
        logger.warning("Attempting to perform BLIND merge")
        # resetting all omics indexes to one common index
        given_omics = list(adata_dict.keys())
        logger.info(f"given omics: {given_omics}")
        for i in range(1, len(given_omics)):
            logger.info(f"checking if {given_omics[i]} has same pixels as {given_omics[0]}")
            logger.info(f"pixel count for {given_omics[i]} and {given_omics[0]} \
                            is {len(adata_dict[given_omics[i]].obs.index)} and {len(adata_dict[given_omics[0]].obs.index)}")
            if len(adata_dict[given_omics[i]].obs.index) == len(adata_dict[given_omics[0]].obs.index):

                adata_dict[given_omics[i]].obs == adata_dict[given_omics[0]].obs
            else:
                same_pixles = False
                break
        
    if not same_pixles:
        logger.warning("pixel count is different in both the datasets, BLIND merge is not possible")
        error_message = "Pixel count is different in both the datasets, BLIND merge is not possible"
        return error_message

    adata=ad.concat(adata_dict,join='inner',axis=1,label='omics')


    obs=data_temp.set_index('obsindex')[['region']]
    obs=obs[obs.index.isin(adata.obs.index)]
    adata.obs=obs
    adata.obsm['spatial']=adata.obs.reset_index()['obsindex'].str.split('x',expand=True).astype(np.float64).to_numpy()
    
    if split==True:
        unique_id=data_temp['region'].unique().tolist()
        for i in range(len(unique_id)):
            adata_sub=adata[adata.obs['region']==unique_id[i]].copy()
            
            # # writing only if it have more than 1 omics, because single omics data is already written in csv2h5ad step
            if len(adata_sub.var['omics'].unique())>1:
                # logger.info(f'writing {adata_sub.var["omics"].unique()}')
                adata_sub.write(os.path.join(data_saving_path,f'{unique_id[i]}_pool.h5ad'))
    
    # this contains all tissues and all omics data
    adata.write(os.path.join(data_saving_path,f'All_tissues_pool.h5ad'))
    logger.info(f'pooling completed')

    return error_message
    


