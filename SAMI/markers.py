import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import ranksums
from .utils import *
import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import seaborn as sns
import logging

logger = logging.getLogger(__name__)
class Markers:
    def __init__(self, region,):
        self.region = region
        self.list_of_clusters = None

    
    def findmarkers(self,dict_clustermap, adata,adj_pval_cutoff,top,adata2=None, workdir=None):  
        initial_id = 70
        logger.info("Markers core started")
        logger.info(f'Finding markers for {self.region}')
        cluster_param='leiden'
        
        pvalue_result=pd.DataFrame()
        marker_result=pd.DataFrame()
        cluster_list=adata.obs[cluster_param].cat.categories.tolist()
        self.list_of_clusters = cluster_list.copy()
        cluster_list = sorted([int(s) for s in cluster_list])
        logger.info("Total clusters: "+str(len(cluster_list)))

        left_a = dict_clustermap['left']
        right_b = dict_clustermap['right']

        previos_clusters = []
        for left_a, right_b in zip(left_a, right_b): #[0]:
            logger.info(f'Finding markers for cluster {left_a} vs {right_b}')

            # if its a general single 
            if len(left_a)==1 and left_a[0] not in previos_clusters:
                cluster_id = left_a[0]
                previos_clusters.append(cluster_id)
            else:
                cluster_id = initial_id
                initial_id += 1

            
            pvalue_temp=pd.DataFrame(columns=['feature','omics','pvalue','avg_log2FC','abs_avg_log2FC','pct1','pct2','adj_pvalue','cluster','rank', 'left_a', 'right_b'])
            for feature in adata.var.index:
                if adata2 is not None and feature in adata2.var.index:
                    # logger.info("dual")
                    a=adata[adata.obs[cluster_param].isin(left_a)][:,feature].X.copy()
                    b=adata2[adata2.obs[cluster_param].isin(right_b)][:,feature].X.copy()

                    # logger.debug(f"left shape {a.shape} right shape {b.shape}")
                    # logger.debug(f"left: {list(adata[adata.obs[cluster_param].isin(left_a)].obs['region'].unique())} right: {list(adata2[adata2.obs[cluster_param].isin(right_b)].obs['region'].unique())}")

                else:
                    # logger.info("single")
                    a=adata[adata.obs[cluster_param].isin(left_a)][:,feature].X.copy()
                    b=adata[adata.obs[cluster_param].isin(right_b)][:,feature].X.copy()
                    # b=adata[adata.obs[cluster_param]==str(9)][:,feature].X.copy()
                    # logger.debug(f"left: {list(adata[adata.obs[cluster_param].isin(left_a)].obs['region'].unique())} right: {list(adata[adata.obs[cluster_param].isin(right_b)].obs['region'].unique())}")
                    # logger.debug(f"left shape {a.shape} right shape {b.shape}")
                
                # logger.info(f"left shape {a.shape} right shape {b.shape}")

                # casting the array fro 64 to 32 to save sapce
                a = a.astype(np.float32)
                b = b.astype(np.float32)
                    
                if len(a)!=0 and len(b)!=0:
                    pct1=sum(a!=0)/len(a)
                    pct2=sum(b!=0)/len(b)
                    if (pct1>0.1)&(pct2>0.1):
                        stat,pvalue=ranksums(a,b)
                        pvalue=pvalue[0]
                        mean_a=np.mean(a)
                        mean_b=np.mean(b)
                        avg_log2FC=np.log2(mean_a/mean_b)
                        abs_avg_log2FC=abs(avg_log2FC)
                        omics = adata.var.loc[feature,'omics']
                        pvalue_feature=pd.DataFrame({'feature':feature,'omics':omics,'pvalue':[pvalue],'avg_log2FC':avg_log2FC,
                                                     'abs_avg_log2FC':abs_avg_log2FC,'pct1':pct1,'pct2':pct2,'cluster':cluster_id, 
                                                     'left_a': ','.join(left_a), 'right_b': ','.join(right_b)})
                        # logger.info(f"cluster_id {cluster_id}")
                        pvalue_temp=pd.concat([pvalue_temp,pvalue_feature],ignore_index=True)
                    else:
                        # logger.info("continue hit")
                        continue
                        

            #results of all the features
            pvalue_result=pd.concat([pvalue_result,pvalue_temp],ignore_index=True)
            # logger.info(f"pvalue_result shape: {pvalue_result['cluster'].unique()}")

            #results of markers
            result_temp=pvalue_temp.sort_values(by=['pvalue'],ascending=True).reset_index(drop=True)
            result_temp=adj_pvalue(result_temp).sort_values(by=['adj_pvalue'],ascending=True).loc[result_temp['adj_pvalue']<=adj_pval_cutoff,:]
            result_temp=result_temp.sort_values(by=['abs_avg_log2FC'],ascending=False).reset_index(drop=True).head(top)
            marker_result=pd.concat([marker_result,result_temp],ignore_index=True)

        pvalue_result=pvalue_result.drop(columns=['rank','pct1','pct2','adj_pvalue'])
        marker_result=marker_result.drop(columns=['rank'])

        if not os.path.exists(workdir):
            os.mkdir(workdir)

        

        # crude updates
        marker_result['omics'] = marker_result['omics'].replace({'lip':'lipidomics','sm':'metabolomics'})
        marker_result['omics'] = marker_result['omics'].apply(lambda x: 'lipidomics' if x.startswith('lip') else x)
        marker_result['omics'] = marker_result['omics'].apply(lambda x: 'glycomics' if x.startswith('gly') else x)
        marker_result['omics'] = marker_result['omics'].apply(lambda x: 'metabolomics' if x.startswith('met') else x)
        
        pvalue_result.to_csv(os.path.join(workdir,f'{self.region}_pvalue.csv'),index=False)
        marker_result.to_csv(os.path.join(workdir,f'{self.region}_marker.csv'),index=False)
        logger.info("Markers core finished")
        
    def circular_tree(self,clusters=None,top_n=5, show=True, workdir=None):
        markers = pd.read_csv(os.path.join(workdir,f'{self.region}_marker.csv'))
        if clusters != None:
            markers = markers.loc[markers['cluster'].isin(clusters)]
        else:
            clusters = 'all'
        
        top = markers.groupby(['cluster'],as_index=False).head(top_n)

        top['c|f'] = top['cluster'].astype(str)+'|'+top['feature']
        top['center'] = 'center'


        renames = [['c|f','cluster'],['cluster','center']]
        data1 = pd.DataFrame()
        for rename in renames:
            temp = top.copy()
            temp['start'], temp['end'] = temp[rename[0]].copy(), temp[rename[1]].copy()
            temp = temp.drop_duplicates(subset=['start','end'],keep='first').astype(str)
            data1 = pd.concat([data1,temp],axis=0)
        data2 = data1.copy()
        data2 = data2.rename(columns={'start':'end','end':'start'})
        data = pd.concat([data1,data2],axis=0)

        data['label'] = data.apply(lambda x: '' if x['start']=='center' else x['feature'] if x['start'].count('|') == 1 else x['start'] , axis=1)
        data['size'] = data.apply(lambda x: x['abs_avg_log2FC'] if x['start'].count('|') == 1 else 5, axis=1)
        color_list = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        data['color'] = data.apply(lambda x: str(80) if 'center' in x['start'] else str(int(x['cluster'])) , axis=1)
        data['color'] = data['color'].map(color_list)

        plt.figure(figsize=(15, 15))
        G = nx.Graph()

        nodes = pd.unique(data[['start', 'end']].values.ravel('K'))
        G.add_nodes_from(nodes)

        edges = data[['start', 'end']].values
        G.add_edges_from(edges)

        size_dict = dict(zip(data['start'], data['size']))
        nx.set_node_attributes(G, size_dict, 'size')

        color_dict = dict(zip(data['start'], data['color']))
        nx.set_node_attributes(G, color_dict, 'color')


        label_dict = dict(zip(data['start'], data['label']))
        nx.set_node_attributes(G, label_dict, 'label')

        pos = nx.nx_agraph.graphviz_layout(G, prog="twopi", args="")

        center = pos['center']
        label_pos = {}
        for key,cord in pos.items():
            label_pos[key] = ((cord[0]-center[0])*1.05+center[0],(cord[1]-center[1])*1.03+center[1])

        labels = nx.draw_networkx_labels(G, label_pos, label_dict, font_size=15, font_weight='bold')
        theta = {k: np.arctan2(v[1]-center[1], v[0]-center[0]) * 180/np.pi for k, v in label_pos.items() }
        for key,t in labels.items():
            if key.count('|')==1:
                if 90 < theta[key] or theta[key] < -90 :
                    angle = 180 + theta[key]
                    t.set_ha('right')
                else:
                    angle = theta[key]
                    t.set_ha('left')
                t.set_va('center')
                t.set_rotation(angle)
                t.set_rotation_mode('anchor')
            else:
                if 90 < theta[key] or theta[key] < -90 :
                    angle = 180 + theta[key]
                    t.set_ha('right')
                else:
                    angle = theta[key]
                    t.set_ha('left')
                t.set_va('center')
                t.set_rotation(angle)
                t.set_rotation_mode('anchor')

        nx.draw(G, pos, node_color=[d['color'] for n, d in G.nodes(data=True)],node_size=[float(d['size'])*40 for n, d in G.nodes(data=True)],edge_color='#D3D3D3')
        x_values, y_values = zip(*pos.values())
        x_max = max(x_values)
        x_min = min(x_values)
        y_max = max(y_values)
        y_min = min(y_values)
        x_margin = (x_max - x_min) * 0.5
        y_margin = (y_max - y_min) * 0.5
        plt.xlim(x_min - x_margin, x_max + x_margin)
        plt.ylim(y_min - y_margin, y_max + y_margin)

        cluster_str = str(clusters).replace("'","")
        plt.savefig(os.path.join(workdir,f'circulartree_{self.region}_{cluster_str}_top_{top_n}.png'))
        if show==False:
            plt.close()
        else:
            plt.show()
        return os.path.join(workdir,f'circulartree_{self.region}_{cluster_str}_top_{top_n}.png')
    
    def get_clusters(self, work_dir):
        pvalues = pd.read_csv(os.path.join(work_dir,f'{self.region}_pvalue.csv'))
        return list(pvalues['cluster'].unique())
    
    def volcano_plot(self,cluster,show=True, work_dir=None):
        pvalues = pd.read_csv(os.path.join(work_dir,f'{self.region}_pvalue.csv'))

        pval = pvalues.loc[pvalues['pvalue']>0]
        pval = pval.loc[(pval['cluster']==cluster)&(pval['pvalue']>=np.percentile(pval['pvalue'],5))&(pval['abs_avg_log2FC']<=np.percentile(pval['abs_avg_log2FC'],95))]

        x = pval['avg_log2FC']
        y = -np.log10(pval['pvalue'])

        colors = np.where(y <= -np.log10(0.05),'darkgray',np.where(x >= 0, 'red', 'blue'))

        fig, ax = plt.subplots(figsize=(10,10))
        fig.subplots_adjust(left=0.2,bottom=0.2)
        plt.scatter(x,y,c=colors,s=10)

        crit1=pval[(pval['avg_log2FC']<=0.5)&(pval['avg_log2FC']>=0)].sort_values('pvalue').head(2).index
        crit2=pval[(pval['avg_log2FC']>=-0.5)&(pval['avg_log2FC']<=0)].sort_values('pvalue').head(2).index
        ann_df = pval.loc[pval.index.isin(crit1)|pval.index.isin(crit2)]
        ann_dict = {}
        for index, row in ann_df.iterrows():
            x = row['avg_log2FC']
            y = -np.log10(row['pvalue'])
            annotation = row['feature']
            coord = (x, y)
            ann_dict[coord] = annotation

        texts=[]
        for key, values in ann_dict.items():
            texts.append(plt.text(key[0], key[1], values, size=23, weight='bold')) 
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1))

        x_min = np.min(pval['avg_log2FC'])
        x_max = np.max(pval['avg_log2FC'])
        x_abs_max = max(abs(x_min), abs(x_max))*1.1
        plt.xlim(-x_abs_max, x_abs_max)

        plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--',alpha=0.5)

        plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=2))
        plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=1))
        plt.xticks(fontsize=40,fontweight='bold')
        plt.yticks(fontsize=40,fontweight='bold')

        plt.xlabel('log2(FC)',fontsize=40,fontweight='bold')
        plt.ylabel('-log10(pvalue)',fontsize=40,fontweight='bold')

        cluster = str(cluster).replace("'","")
        plt.savefig(os.path.join(work_dir,f'vol_{self.region}_{cluster}.png'))
        if show == False:
            plt.close()
        else:
            plt.show()

        return os.path.join(work_dir,f'vol_{self.region}_{cluster}.png')
            
    def density_plot(self,adata,cluster,marker,show=False):
        adata= adata[:,marker]
        data = pd.concat([adata.to_df(),adata.obs],axis=1)
        data['cluster']=np.where(data['leiden']==str(cluster),data['leiden'],'other')

        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        clusters_colors['other']='#A0A0A0'

        fig,ax = plt.subplots(figsize=(12,5))
        sns.kdeplot(data=data,x=marker,hue='cluster',palette=clusters_colors,fill=True,common_norm=False,alpha=0.5,linewidth=0)
        plt.setp(ax.get_legend().get_texts(), fontsize='15') 
        plt.setp(ax.get_legend().get_title(), fontsize='15',fontweight='bold')

        ax.set_xlabel(marker,fontsize=15,fontweight='bold')
        ax.set_ylabel('Density',fontsize=15,fontweight='bold')
        ax.xaxis.set_major_locator(ticker.MaxNLocator(3))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
        plt.savefig(os.path.join('../results/markers/',f'density_{self.region}_{cluster}_{marker}.png'))
        if show == False:
            plt.close()
        else:
            plt.show()
            
            
