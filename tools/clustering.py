from SAMI.clustering import Clusters, clustering, plot_cluster, plot_umap_cluster
from SAMI.utils import adata_concat
from pydantic import BaseModel, FilePath, Field
from chainlit import make_async
from enum import StrEnum
import os
from tools.types import ToolResultType, ToolResult
import logging
import scanpy as sc

logger = logging.getLogger(__name__)

class OmicsModality(StrEnum):
    """Enumeration for omics modalities."""

    pool = "pool"
    metabolomics = "metabolomics"
    lipidomics = "lipidomics"
    glycomics = "glycomics"

class ClusteringMeta(BaseModel):
    file: FilePath = Field(
        description="Path to the .h5ad file containing the data.",
    )
    regions: list[str] = Field(
        description="List of regions to be clustered.",
    )
    modality: list[OmicsModality] = Field(
        description="List of omics data (e.g., pool, metabolomics, lipidomics, glycomics) for the corresponding region in regions based on index. If not explicitly mentioned then it'd be pool.",
    )
    resolution: float = Field(
        description="Resolution parameter for Leiden clustering.",
        default=1.5,
    )
    harmony: bool = Field(
        description="Whether Harmony was used for batch effect correction.",
        default=False,
    )


class ClusteringInput(BaseModel):
    region_list: list[str] = Field(
        description="List of regions to be clustered.",
    )
    modality: list[OmicsModality] = Field(
        description="list of omics data (e.g., pool, metabolomics, lipidomics, glycomics). for the corresponding region in region_list based on index. If not explicitely mentioned then it'd be pool.",
    )
    resolution: float = Field(
        description="Resolution parameter for Leiden clustering.",
        default=1.5,
    )
    harmony: bool = Field(
        description="Whether to use Harmony for batch effect correction.",
        default=False,
    )


async def clustering_tool(input: ClusteringInput, session_id: str) -> list[ToolResult]:
    """
    Perform clustering on the given input data.

    Args:
        input (ClusteringInput): Input data for clustering.
        session_id (str): Session ID for tracking.

    Returns:
        str: Result of the clustering operation.
    """
    results = []
    work_dir = os.path.join(".files", session_id)
    pool_dir = os.path.join(work_dir, "pooled_data")
    results_dir = os.path.join(work_dir, "results")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    if len(input.region_list) != len(input.modality):
        logger.error("Length of region_list and modality must be the same.")
        return [
            ToolResult(
                type=ToolResultType.text,
                content="Length of region_list and modality must be the same.",
                error=True
            )
        ]
    elif len(input.region_list) == len(input.modality) == 1:
        if input.region_list[0] == "All_tissues" and input.modality[0] != "pool":
            logger.error("If region_list is 'All_tissues', modality must be 'pool'.")
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content="If region_list is 'All_tissues', modality must be 'pool'.",
                    error=True
                )
            ]
        print(f"Reading file: {input.region_list[0]}_{input.modality[0]}.h5ad")
        try:
            read_file = sc.read(os.path.join(pool_dir, f"{input.region_list[0]}_{input.modality[0]}.h5ad"))
            regions = input.region_list[0]
            modality = input.modality[0]
        except FileNotFoundError:
            logger.error(f"File {input.region_list[0]}_{input.modality[0]}.h5ad does not exist.")
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content=f"File {input.region_list[0]}_{input.modality[0]}.h5ad does not exist.",
                    error=True
                )
            ]
    else:
        h5ad_files = []
        for i in range(len(input.region_list)):
            if input.region_list[i] == "All_tissues" and input.modality[i] != "pool":
                logger.error("If region_list is 'All_tissues', modality must be 'pool'.")
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        content="If region_list is 'All_tissues', modality must be 'pool'.",
                        error=True
                    )
                ]
            file_name = f"{input.region_list[i]}_{input.modality[i]}.h5ad"
            if not os.path.exists(os.path.join(pool_dir, file_name)):
                logger.error(f"File {file_name} does not exist.")
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        content=f"File {file_name} does not exist.",
                        error=True
                    )
                ]
            h5ad_files.append(sc.read(os.path.join(pool_dir, file_name)))
        try:
            read_file = adata_concat(h5ad_files)
        except Exception as e:
            logger.error(f"Error while concatenating files: {str(e)}")
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content=f"Error while concatenating files: {str(e)}",
                    error=True
                )
            ]
    cluster = Clusters(
        read_file,
        os.path.join(work_dir, "results"),
        resolution=input.resolution,
        harmony_flag=input.harmony,
    )
    try:
        async_clustering = make_async(clustering)
        cluster_file = await async_clustering(cluster)
        
        results.append(
            ToolResult(
                type=ToolResultType.text,
                content=f"Clustering completed successfully. Results saved to {cluster_file}.",
                error=False
            )
        )
    except Exception as e:
        logger.error(f"Error while clustering: {str(e)}")
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error while clustering: {str(e)}",
                error=True
            )
        ]
    try:
        async_plot_umap_cluster = make_async(plot_umap_cluster)
        umap_file = await async_plot_umap_cluster(cluster)
        results.append(
            ToolResult(
                type=ToolResultType.image,
                content=umap_file,
                error=False,
                desc="UMAP plot of the clustered data."
            )
        )
    except Exception as e:
        logger.error(f"Error while plotting UMAP: {str(e)}")
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error while plotting UMAP: {str(e)}",
                error=True
            )
        ]
    try:
        async_plot_cluster = make_async(plot_cluster)
        cluster_image_file = await async_plot_cluster(cluster)
        results.append(
            ToolResult(
                type=ToolResultType.image,
                content=cluster_image_file,
                error=False,
                desc="Cluster plot of the clustered data."
            )
        )            
    except Exception as e:
        logger.error(f"Error while plotting cluster: {str(e)}")
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error while plotting cluster: {str(e)}",
                error=True
            )
        ]
    return results
