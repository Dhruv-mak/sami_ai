import os
import scanpy as sc
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
import numpy as np
import logging
import harmonypy as hm
import time
import uuid
logger = logging.getLogger(__name__)

class Clusters:
    
    def __init__(self, adata, results_dir, resolution, harmony_flag=True):
        self.adata = adata
        self.result_dir = results_dir
        self.resolution = resolution
        self.harmony_flag = harmony_flag
    
def clustering(self):

    logger.info("Clustering started")

    try:
        logger.info("CLustering:  Harmony flag: {}".format(self.harmony_flag))
        sc.pp.pca(self.adata)
        if (len(np.unique(self.adata.obs['region']))>1) and self.harmony_flag == 1:
            sc.external.pp.harmony_integrate(self.adata,key='region',max_iter_harmony=20)
            all_pcs = self.adata.obsm['X_pca_harmony'].shape[1]
            n_pcs = min(20, all_pcs)
            sc.pp.neighbors(self.adata,use_rep = 'X_pca_harmony', n_pcs=n_pcs)
        else:
            all_pcs = self.adata.obsm['X_pca'].shape[1]
            n_pcs = min(20, all_pcs) 
            sc.pp.neighbors(self.adata,n_pcs=n_pcs)
            
        sc.tl.umap(self.adata, random_state= 7)
        sc.tl.leiden(self.adata,resolution=self.resolution, random_state= 7)
    except Exception as e:
        logger.error(f"Error in clustering core SAMI: {str(e)}")
        return e

    self.filename = f'{uuid.uuid4().hex[:8]}.h5ad'
    self.adata.write(os.path.join(self.result_dir,self.filename))
    logger.info("Clustering completed and written the output file")

    return self.filename


def plot_umap_cluster(self,size=50,show=False):
    color='leiden'
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist=self.adata.obs[color].cat.categories.tolist()
    n_cluster=str(len(colorlist))

    fig, axs = plt.subplots(1, 2, figsize=(30, 10))

    sc.pl.umap(
        self.adata, 
        legend_loc='on data',
        color=[color], 
        palette=[v for k, v in clusters_colors.items() if k in colorlist],
        ax=axs[0],
        show=False,
    )

    sc.pl.spatial(
        self.adata,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color=color,
        title=f'UMAP Plot',
        size=size,
        palette=[v for k, v in clusters_colors.items() if k in colorlist],
        ax=axs[1],
        show=False
    )

    samples = np.unique(self.adata.obs['region']).tolist()
    for sample in samples:
        adata_sub = self.adata[self.adata.obs['region']==sample]
        x_label = np.min(adata_sub.obsm['spatial'][:,0])
        y_label = np.min(adata_sub.obsm['spatial'][:,1])
        axs[1].text(x_label,y_label, sample,fontsize=10, color='black')
    
    filename = f'{uuid.uuid4().hex[:8]}_umap.png'
    plt.savefig(os.path.join(self.result_dir, filename),dpi=400)
    
    if show==False:
        plt.close()
    else:
        plt.show()

    # return os.path.join(self.result_dir, f'{self.filename}_umap.png'), len(colorlist)
    return os.path.join(self.result_dir, filename)

def plot_cluster(self,size=50,show=False):
    color='leiden'
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist=self.adata.obs[color].cat.categories.tolist()
    n_cluster=str(len(colorlist))
    fig, axs = plt.subplots(1, 1, figsize=(15, 10))
    sc.pl.spatial(
        self.adata,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color=color,
        title=f'Clustering Plot',
        size=size,
        palette=[v for k, v in clusters_colors.items() if k in colorlist],
        show=False
    )

    samples = np.unique(self.adata.obs['region']).tolist()
    for sample in samples:
        adata_sub = self.adata[self.adata.obs['region']==sample]
        x_label = np.min(adata_sub.obsm['spatial'][:,0])
        y_label = np.min(adata_sub.obsm['spatial'][:,1])
        axs.text(x_label,y_label, sample,fontsize=10, color='black')

    filename = f'{uuid.uuid4().hex[:8]}.png'
    plt.savefig(os.path.join(self.result_dir,filename),dpi=400)
    if show==False:
        plt.close()
    else:
        plt.show()
    return os.path.join(self.result_dir,filename)

def plot_select_cluster(self,cluster,size=50,show=False):
    adata = sc.read(os.path.join(self.file_path,f'{self.region}_{self.modality}_{self.resolution}.h5ad'))
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    
    if isinstance(cluster, int):
        cluster = [cluster]
        
    colorlist=adata.obs['leiden'].cat.categories.tolist()
    palette=[v if k in map(str,cluster) else '#A0A0A0' for k, v in clusters_colors.items()]
    fig, axs = plt.subplots(1, 1, figsize=(15, 10))

    sc.pl.spatial(
    adata,
    img_key=None,
    library_id=None,
    spot_size=1, #important to add to avoid error
    color='leiden',
    legend_loc='none',
        title=f'{self.region}_{self.modality}_{self.resolution}_{cluster}',
    size=50,
    palette=palette,
    ax=axs,
    show=False
    )

    samples = np.unique(adata.obs['region']).tolist()
    for sample in samples:
        adata_sub = adata[adata.obs['region']==sample]
        x_label = np.min(adata_sub.obsm['spatial'][:,0])
        y_label = np.min(adata_sub.obsm['spatial'][:,1])
        axs.text(x_label,y_label, sample,fontsize=10, color='black')
    plt.savefig(os.path.join(self.file_path,f'{self.region}_{self.modality}_{self.resolution}_{cluster}.png'),dpi=400) 
    if show==False:
        plt.close()
    else:
        plt.show()

    return os.path.join(self.file_path,f'{self.region}_{self.modality}_{self.resolution}_{cluster}.png')
    
class Cluster_Integration:
    def __init__(self,adata,adata_ref, region1, region2, working_folder,clusterInt_specific_working_folder):

        self.adata1 = adata
        self.adata2 = adata_ref
        self.working_folder = working_folder
        self.region1 = region1
        self.region2 = region2
        # this is the folder where the integrated data and results will be saved
        self.clusterInt_specific_working_folder = clusterInt_specific_working_folder
        
def integrate(self):

    logger.info("INtegration started")
    var_names = self.adata1.var_names.intersection(self.adata2.var_names)
    self.adata1.obs['sample'] = self.region1
    self.adata2.obs['sample'] = self.region2

    try:
        sc.tl.ingest(self.adata1,self.adata2,obs='leiden')
    except Exception as e:
        logger.error(f"Error in integration core SAMI: {str(e)}")
        return
    
    file_path1 = os.path.join(self.clusterInt_specific_working_folder, f'{time.strftime("%m%d%y_%H%M%S")}_integrated.h5ad')
    self.adata1.write(file_path1)
    time.sleep(1)
    file_path2 = os.path.join(self.clusterInt_specific_working_folder, f'{time.strftime("%m%d%y_%H%M%S")}_integrated.h5ad')
    self.adata2.write(file_path2)
    logger.info("Integration completed and written the output file")
    return f"Integration completed and written the output file: {file_path1} and {file_path2}"
    
def splitname(self,adata):
    # name,_ = os.path.splitext(adata)
    # region, modality, resolution = name.split('_')
    # return region, modality, resolution
    # changed for GUI braoder handling
    name = adata.replace('.h5ad', "")
    #region, modality, resolution = name.split('_')
    resolution = name.split('_')[-1]
    modality = name.split('_')[-2]
    region = '_'.join(name.split('_')[:-2])

    return region, modality, resolution

def plot_umap_cluster_int(self,adata, size=50,show=False):
    color='leiden'
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist=adata.obs[color].cat.categories.tolist()
    n_cluster=str(len(colorlist))
    fig, axs = plt.subplots(1, 2, figsize=(30, 10))

    sc.pl.umap(
        adata, 
        legend_loc='on data',
        color=[color], 
        title="UMAP Plot",
        palette=[v for k, v in clusters_colors.items() if k in colorlist],
        ax=axs[0],
        show=False,
    )

    sc.pl.spatial(
        adata,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color=color,
        title="UMAP Plot",
        size=size,
        palette=[v for k, v in clusters_colors.items() if k in colorlist],
        ax=axs[1],
        show=False
    )
    samples = np.unique(adata.obs['region']).tolist()
    for sample in samples:
        adata_sub = adata[adata.obs['region']==sample]
        x_label = np.min(adata_sub.obsm['spatial'][:,0])
        y_label = np.min(adata_sub.obsm['spatial'][:,1])
        axs[1].text(x_label,y_label, sample,fontsize=10, color='black')

    filename = f'{uuid.uuid4().hex[:8]}_umap_integrated.png'
    fig_path = os.path.join(filename)
    plt.savefig(fig_path,dpi=400)

    if show==False:
        plt.close()
    else:
        plt.show()

    return fig_path
        
def plot_overlap_umap_int(self, adata1, adata2, int_folder, show=False):
    adata_concat = adata2.concatenate(adata1, batch_categories=[self.region2,self.region1])
    fig,ax = plt.subplots(figsize=(10,8))
    plt.subplots_adjust(left=0.1,right=0.7)
    sc.pl.umap(adata_concat,color='batch',title=f'Overlay UMAP for integration',size=5,legend_fontsize=20,ax=ax,show=False)
    
    # fig_path = os.path.join(int_folder, f'{time.strftime("%m%d%y_%H%M%S")}_umap_overlap.png')
    fig_path = os.path.join(int_folder, f'{uuid.uuid4().hex[:8]}_umap_overlap.png')
    plt.savefig(fig_path,dpi=200)
    if show==False:
        plt.close()
    else:
        plt.show()
    return fig_path
        
        
def plot_select_cluster_int(self,region,modality, res, cluster,size=50,show=False):
    adata = sc.read(os.path.join(self.clusterInt_specific_working_folder,f'{region}-{modality}-{res}_integrated.h5ad'))
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    
    if isinstance(cluster, int):
        cluster = [cluster]
        
    colorlist=adata.obs['leiden'].cat.categories.tolist()
    palette=[v if k in map(str,cluster) else '#A0A0A0' for k, v in clusters_colors.items()]
    fig, axs = plt.subplots(1, 1, figsize=(15, 10))
    sc.pl.spatial(
        adata,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color='leiden',
        legend_loc='none',
        title=f'{region}_integrated_{cluster}',
        size=size,
        palette=palette,
        show=False,
        ax = axs
    )

    samples = np.unique(adata.obs['region']).tolist()
    for sample in samples:
        adata_sub = adata[adata.obs['region']==sample]
        x_label = np.min(adata_sub.obsm['spatial'][:,0])
        y_label = np.min(adata_sub.obsm['spatial'][:,1])
        axs.text(x_label,y_label, sample,fontsize=10, color='black')

    plt.savefig(os.path.join(self.clusterInt_specific_working_folder,f'{region}-{modality}-{res}_{cluster}.png'),dpi=400) 
    if show==False:
        plt.close()
    else:
        plt.show()


    return os.path.join(self.clusterInt_specific_working_folder,f'{region}-{modality}-{res}_{cluster}.png')