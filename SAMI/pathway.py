import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import networkx as nx
import logging
logger = logging.getLogger(__name__)
try:
    # PyInstaller creates a temp folder and stores path in _MEIPASS
    R_HOME = os.path.join(sys._MEIPASS, "R_HOME")
    path = os.path.join(R_HOME, "bin")
    os.environ["R_HOME"] = R_HOME
    os.environ["PATH"] += path + ";" + os.environ["R_HOME"]
    logger.info("running in executable mode")
except Exception as e:
    os.environ["PATH"] += os.pathsep + os.path.join(os.environ["R_HOME"], "bin")
    
    logger.info("running in dev mode")

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


class Pathway:
    def __init__(self,region,modality, pathways_folder, makers_folder=None):
        
        self.pathways_folder = pathways_folder
        self.markers_folder = makers_folder
        if not os.path.exists(pathways_folder):
            os.mkdir(pathways_folder)
            
        self.region = region
        if modality=="metabolomics_glycomics":
            self.modality = 'metabolomics_glycomics'
        elif modality=='metabolomics' or modality=='glycomics':
            print('Metabolomics and Glycomics are combined for pathway enrichment analysis.')
            self.modality = 'metabolomics_glycomics'
        else:
            self.modality = 'lipidomics'
       
def findpathway(self):

    result = []
    logger.info("Pathway enrichment core started")
    try:
        base_path = sys._MEIPASS
        libraries_path = os.path.dirname(base_path)
        libraries_path = os.path.join(libraries_path, r"R-4.3.3/win-library/4.3")
        libraries_path = str(libraries_path).replace("\\", "/")
        logger.info("running in executable mode")
        logger.info(f"libraries path: {libraries_path}")
    except:
        # libraries_path = "C:/Program Files/R/R-4.3.3/win-library/4.3"
        # libraries_path = os.path.join(os.path.expanduser("~"), "R", "x86_64-redhat-linux-gnu-library", "4.4")
        # libraries_path = os.environ.get("R_LIBS_USER", "")
        logger.info("running in dev mode")
        # logger.info(f"libraries path: {libraries_path}")

    
    # ro.r(f'custom_lib_path <- "{libraries_path}"')
    # ro.r(f'Sys.setenv(R_LIBS_USER = custom_lib_path)')
    # ro.r('.libPaths(c(custom_lib_path, .libPaths()))')

    r_source = ro.r['source']
    logger.info("current working directory: {}".format(os.getcwd()))
    # logger.info("Rlocation: {}".format(R_location))
    logger.info("Pathway enrichment started")
    try:
    # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = os.path.join(sys._MEIPASS , 'config')
    except Exception:
        base_path = os.path.join(os.getcwd(), 'SAMI')
    
    logger.info("base path of R file: {}".format(base_path))
            # finding lib folder
    try:
    # PyInstaller creates a temp folder and stores path in _MEIPASS
        lib_folder = os.path.join(sys._MEIPASS , 'lib')
    except Exception as e:
        lib_folder = os.path.join(os.getcwd(), 'lib')
    finally:
        logger.info(f"lib folder: {lib_folder}")

    r_file =  os.path.join(base_path, 'MetabAnalystR.R')
    r_source(r_file) 
    
    markers = pd.read_csv(os.path.join(self.markers_folder,f'{self.region}_marker.csv'))
    clusters = markers['cluster'].unique().tolist()
    
    if self.modality=='metabolomics_glycomics':
        markers = markers.loc[markers['omics'].isin(['metabolomics','glycomics'])]
        
        result_pathR = ro.StrVector([self.pathways_folder])
        logical_arg = ro.BoolVector([True,False])
        modalityR = ro.StrVector([self.modality])
        clustersR = ro.vectors.IntVector(clusters)
        regionR = ro.StrVector([self.region])
        with localconverter(ro.default_converter + pandas2ri.converter):
            markersR = ro.conversion.py2rpy(markers)

        PathwayEnrichment = ro.r['PathwayEnrichment']
        
        # PathwayEnrichment(result_path=result_pathR,modality=modalityR,region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[1])
        PathwayEnrichment(region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[1],merge=logical_arg[1],
                            results_folder=self.pathways_folder , lib_folder = lib_folder)
        
        # (region,markers,clusters,lipid,merge, results_folder, lib_folder)
        
        logger.info("pathway enrichment analysis is done.")
        #### 2024-08-31 add direction ###
        pathways = pd.read_csv(os.path.join(self.pathways_folder,f'ora_{self.region}_{self.modality}.csv'))
        
        pathways['pathway_avg_log2FC'] = None
        pathways['regulated'] = None
        for idx, row in pathways.iterrows():
            cluster = row['cluster']
            compounds = row['compounds'].split(',')
            markers_fc = markers.loc[(markers['cluster']==cluster)&(markers['feature'].isin(compounds))]
            pathway_fc = round(markers_fc['avg_log2FC'].mean(),4)
            pathways.at[idx, 'pathway_avg_log2FC'] = pathway_fc
            pathways.at[idx, 'regulated'] = 'up' if pathway_fc > 0 else 'down'
        csv_file = f'ora_{self.region}_{self.modality}.csv'
        csv_file_path =  os.path.join(self.pathways_folder,csv_file)
        pathways.to_csv(csv_file_path,index=False)
        return csv_file
                
    else:
        markers = markers.loc[markers['omics']=='lipidomics']
    
        result_pathR = ro.StrVector([self.pathways_folder])
        logical_arg = ro.BoolVector([True,False])
        modalityR = ro.StrVector([self.modality])
        clustersR = ro.vectors.IntVector(clusters)
        regionR = ro.StrVector([self.region])
        with localconverter(ro.default_converter + pandas2ri.converter):
            markersR = ro.conversion.py2rpy(markers)

        PathwayEnrichment = ro.r['PathwayEnrichment']
    
        # PathwayEnrichment(result_path=result_pathR,modality=modalityR,region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[0])
        
        PathwayEnrichment(region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[0],merge=logical_arg[1], 
                            results_folder=self.pathways_folder, lib_folder = lib_folder)
        #### 2024-08-31 add direction ###
        pathways = pd.read_csv(os.path.join(self.pathways_folder,f'ora_{self.region}_{self.modality}.csv'))
        
        pathways['pathway_avg_log2FC'] = None
        pathways['regulated'] = None
        for idx, row in pathways.iterrows():
            cluster = row['cluster']
            compounds = row['compounds'].split(',')
            markers_fc = markers.loc[(markers['cluster']==cluster)&(markers['feature'].isin(compounds))]
            pathway_fc = round(markers_fc['avg_log2FC'].mean(),4)
            pathways.at[idx, 'pathway_avg_log2FC'] = pathway_fc
            pathways.at[idx, 'regulated'] = 'up' if pathway_fc > 0 else 'down'
        
        csv_file = f'ora_{self.region}_{self.modality}.csv'
        csv_file_path =  os.path.join(self.pathways_folder,csv_file)
        pathways.to_csv(csv_file_path,index=False)
        return csv_file
    
def anndict(self):
    anndict = {}
    annoations_directory = os.path.join(*[ os.getcwd(), 'visionApp', 'annotation'])
    for file in os.listdir(annoations_directory):
        if re.match(r'^\w+annotation\.csv',file):
            modality = file.split('_')[0]
            data_temp = pd.read_csv(os.path.join(annoations_directory, file),delimiter=',')
            feat_cols = data_temp['Name'].tolist()
            anndict[modality] = feat_cols
    return anndict


def plot_dot(self,cluster,scale,height,top=None,show=False):
    logger.info("pathways dot plot called")
    logger.info("reading_csv {}".format(os.path.join(self.pathways_folder ,f'ora_{self.region}_{self.modality}.csv')))
    ora_data = pd.read_csv(os.path.join(self.pathways_folder ,f'ora_{self.region}_{self.modality}.csv'))
    ora_data['enrichment_ratio'] = ora_data['hits']/ora_data['expected']
    
    if top is None:
        data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)]
    else:
        data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)].head(top)

    if data.shape[0]!=0:

        fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(10,10))
        fig.subplots_adjust(left=0.45,wspace=0.35)
        #### p-value  #####   
        pathway = data['pathway'].tolist()
        log_trans = -np.log10(data['Raw.p'])

        #axis
        y_pos = np.arange(len(pathway))*5
        x_pos = np.arange(0,np.ceil(max(log_trans))*1.1,5 if max(log_trans)>5 else 1)
        ax[0].set_xticks(x_pos)
        ax[0].tick_params(axis='x', labelsize=15)
        ax[0].set_xlim(-np.ceil(max(log_trans))*0.05,np.ceil(max(log_trans))*1.05)
        ax[0].set_yticks(y_pos,pathway,fontsize=18,fontweight='bold')
        ax[0].invert_yaxis()

        size = data['enrichment_ratio']*scale
        x = log_trans.tolist()
        color = data['Raw.p'].tolist()

        scatter = ax[0].scatter(x,y_pos,s=size,c=color,cmap='PuRd_r',zorder=2)
        ax[0].grid(alpha=0.5,zorder=1)
        ax[0].set_xlabel('-log10(p_value)',fontsize=18,fontweight='bold')

        cbar = fig.colorbar(scatter,shrink=0.5)
        cbar.ax.tick_params(labelsize=15)
        cbar.ax.set_title('p-value',fontsize=15,fontweight='bold')
        loc = cbar.ax.get_position()
        new_loc =[loc.x0+0.02,loc.y0+0.15,loc.width,loc.height]
        cbar.ax.set_position(new_loc)

        fig.set_size_inches(15, 8)

        #### FDR  #####   
        pathway = data['pathway'].tolist()
        log_trans = -np.log10(data['FDR'])

        #axis
        y_pos = np.arange(len(pathway))
        x_pos = np.arange(0,max(0.001,np.ceil(max(log_trans))*1.1),5 if max(log_trans)>5 else 1)
        ax[1].set_xticks(x_pos)
        ax[1].tick_params(axis='x', labelsize=15)
        ax[1].set_xlim(-max(np.ceil(max(log_trans))*0.05,0.05),max(1,np.ceil(max(log_trans))*1.1))
        ax[1].set_yticks(y_pos,'',fontsize=10)
        ax[1].invert_yaxis()

        size = data['enrichment_ratio']*scale
        x = log_trans.tolist()
        color = data['FDR'].tolist()

        scatter = ax[1].scatter(x,y_pos,s=size,c=color,cmap='PuRd_r',zorder=2)
        ax[1].grid(alpha=0.5,zorder=1)
        ax[1].set_xlabel('-log10(FDR)',fontsize=18,fontweight='bold')

        handles, labels = scatter.legend_elements(prop="sizes", alpha=1,num=4,func=lambda x: x/scale)
        legend = ax[1].legend(handles, labels, bbox_to_anchor=(1.8,0.43), ncol=1, title="Enrichment\nRatio",fontsize=15,labelspacing=1,borderpad=0.5)
        legend.get_title().set_fontsize(15)
        legend.get_title().set_weight('bold')

        cbar = fig.colorbar(scatter,shrink=0.5)
        scatter.set_clim(vmax=1)
        cbar.ax.tick_params(labelsize=15)
        cbar.ax.set_title('FDR',fontsize=15,fontweight='bold')
        loc = cbar.ax.get_position()
        new_loc =[loc.x0+0.01,loc.y0+0.15,loc.width,loc.height]
        cbar.ax.set_position(new_loc)
        fig.set_size_inches(20,height)

        plt.savefig(os.path.join(self.pathways_folder,f'{self.region}_{self.modality}_{cluster}.png'))
        if show == False:
            plt.close()
        else:
            plt.show()
        return os.path.join(self.pathways_folder,f'{self.region}_{self.modality}_{cluster}.png')
    else:
        return None
    
def plot_bar(self,cluster,top=None,show=False):
    ora_data = pd.read_csv(os.path.join(self.pathways_folder,f'ora_{self.region}_{self.modality}.csv'))
    if top is None:
        up_data = ora_data[(ora_data['regulated'] == 'up')&(ora_data['cluster']==cluster)]
        down_data = ora_data[(ora_data['regulated'] == 'down')&(ora_data['cluster']==cluster)]
    else:
        up_data = ora_data[(ora_data['regulated'] == 'up')&(ora_data['cluster']==cluster)].head(top)
        down_data = ora_data[(ora_data['regulated'] == 'down')&(ora_data['cluster']==cluster)].head(top)

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    if up_data.shape[0]!=0:

        ax[0].barh(up_data['pathway'], -np.log10(up_data['Raw.p']), color='red')
        ax[0].text(1.05, 0.5, 'Up-Regulated', fontsize=15, fontweight='bold', transform=ax[0].transAxes, ha='right', va='center', rotation=270)
        #ax[0].set_title('Up-Regulated',fontsize=15,fontweight='bold')
        ax[0].set_ylabel('Pathway',fontsize=15,fontweight='bold')
        ax[0].invert_yaxis()
        ax[0].set_yticks(range(len(up_data['pathway'])))
        ax[0].set_yticklabels(up_data['pathway'], fontsize=15, weight='bold')
        ax[0].grid(True, which='both', linestyle='--', linewidth=0.5)
        for label in ax[0].get_xticklabels():
            label.set_fontsize(12)
            label.set_fontweight('bold')

    else:
        fig.delaxes(ax[0])

    if down_data.shape[0]!=0:
        ax[1].barh(down_data['pathway'], -np.log10(down_data['Raw.p']), color='blue')
        ax[1].text(1.05, 0.5, 'Down-Regulated', fontsize=15, fontweight='bold', transform=ax[1].transAxes, ha='right', va='center', rotation=270)
        #ax[1].set_title('Down-Regulated',fontsize=15,fontweight='bold')
        ax[1].set_xlabel('-log10(p-value)',fontsize=15,fontweight='bold')
        ax[1].set_ylabel('Pathway',fontsize=15,fontweight='bold')
        ax[1].invert_yaxis()
        ax[1].set_yticks(range(len(down_data['pathway'])))
        ax[1].set_yticklabels(down_data['pathway'], fontsize=15, weight='bold')
        ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)
        for label in ax[1].get_xticklabels():
            label.set_fontsize(12)
            label.set_fontweight('bold')
    else:
        fig.delaxes(ax[1])

    plt.tight_layout()
    plt.savefig(os.path.join(self.pathways_folder,f'{self.region}_{self.modality}_{cluster}_barplot.png'))
    if show == False:
        plt.close()
    else:
        plt.show()
    return os.path.join(self.pathways_folder,f'{self.region}_{self.modality}_{cluster}_barplot.png')


def pathway_network(self,cluster,top=None,size =800, show=False):
    logger.info("pathways network plot called")
    logger.info("reading_csv {}".format(os.path.join(self.pathways_folder ,f'ora_{self.region}_{self.modality}.csv')))
    ora_data = pd.read_csv(os.path.join(self.pathways_folder ,f'ora_{self.region}_{self.modality}.csv'))
    ora_data['enrichment_ratio'] = ora_data['hits']/ora_data['expected']
    
    if self.modality=='metabolomics_glycomics':
        libname = 'smpdb_pathway'
    else:
        libname = 'sub_class'
    
    if top is None:
        data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)]
    else:
        data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)].head(top)

    if data.shape[0]!=0:
        
        # lib_folder = os.path.join(*[ os.getcwd(), 'visionApp', 'lib'])

        #handling lib folder basing on whether its running in pyinstaller or developement environment
        try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
            lib_folder = os.path.join(sys._MEIPASS , 'lib')
        except Exception as e:
            lib_folder = os.path.join(os.getcwd(), 'lib')

        lib = pd.read_csv( os.path.join(lib_folder, f'{libname}.csv'),index_col=0)
        lib['member']=lib['member'].str.split(';')
        lib = lib.explode('member').reset_index(drop=True)

        pathway_list = data['pathway'].tolist()
        lib_temp = lib.loc[lib['name'].isin(pathway_list)]
        pathway_num = len(pathway_list)

        prop = np.zeros((pathway_num,pathway_num))
        for i in range(pathway_num):
            set1 = set(lib.loc[(lib['name']==pathway_list[i]),'member'])
            for j in range(pathway_num):
                set2 = set(lib.loc[(lib['name']==pathway_list[j]),'member'])
                intersection_size = len(set1.intersection(set2))
                union_size = len(set1.union(set2))
                prop[i][j]=intersection_size/union_size
        prop_df=pd.DataFrame(data=prop,index=pathway_list,columns=pathway_list)

        #network
        fig, ax = plt.subplots(figsize=(50,40))
        plt.subplots_adjust(left=0,right=1,bottom=0,top=1)
        G = nx.Graph()

        for i in range(prop_df.shape[0]):
            G.add_node(prop_df.columns[i])
            for j in range(i+1, prop_df.shape[1]):
                if prop_df.iloc[i,j] > 0.25: # set a correlation threshold here
                    G.add_edge(prop_df.columns[i],prop_df.columns[j],weight=prop_df.iloc[i,j]*50)

        pos = nx.circular_layout(G)
        
        ymax = max([coord[1] for coord in pos.values()])
        ymin = min([coord[1] for coord in pos.values()])
        crit = (ymax-ymin)*0.05
        
        for pathway1, coord1 in pos.items():
            if (ymax>coord1[1]>ymax-crit) | (ymin<coord1[1]<ymin+crit):
                coord1[0] *= 1.5

        plt.xlim(min([coord[0] for coord in pos.values()])*1.5,max([coord[0] for coord in pos.values()])*1.5)

        label_pathways = data['pathway'].tolist()
        node_size = []
        node_color = []
        labels = {}
        for node in G.nodes():
            node_size.append(data.loc[data['pathway']==node,'enrichment_ratio']*size)
            node_color.append(data.loc[data['pathway']==node,'Raw.p'])
            if node in label_pathways:
                labels[node] = node

        label_pos = {}
        for key,cord in pos.items():
            if key in label_pathways:
                label_pos[key] = (cord[0],cord[1]-0.05)
                
        weights = [G[u][v]['weight'] for u,v in G.edges()]

        nx.draw_networkx_labels(G, label_pos,labels, font_size=45,font_weight='bold', verticalalignment='top')
        nx.draw(G, pos, node_size=node_size, node_color=node_color,cmap='PuRd_r',edge_color='lightblue',width=weights,edgecolors='black',linewidths=5)

        plt.savefig(os.path.join(self.pathways_folder,f'network_{self.region}_{self.modality}_{cluster}.png'))
        if show == False:
            plt.close()
        else:
            plt.show()

        return os.path.join(self.pathways_folder,f'network_{self.region}_{self.modality}_{cluster}.png')