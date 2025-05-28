import os
from pydantic import BaseModel, FilePath, Field
from tools.types import ToolResultType, ToolResult
from SAMI.clustering import Cluster_Integration, integrate, plot_overlap_umap_int, plot_umap_cluster_int
from chainlit import make_async
import scanpy as sc

class IntegrationInput(BaseModel):
    """
    Input model for the integration tool.
    """

    file_path1: str = Field(
        description="Name of the first .h5ad file to integrate.",
    )
    file_path2: str = Field(
        description="Name of the second .h5ad file to integrate.",
    )

async def integration_tool(input: IntegrationInput, session_id: str) -> list[ToolResult]:
    """
    Integrate two .h5ad files using SAMI's integration method.

    Args:
        input (IntegrationInput): Input data for integration.
        session_id (str): Unique identifier for the session.

    Returns:
        list[ToolResult]: List of results from the integration operation.
    """
    if not input.file_path1.endswith(".h5ad"):
        input.file_path1 += ".h5ad"
    if not input.file_path2.endswith(".h5ad"):
        input.file_path2 += ".h5ad"
    results = []
    results_dir = os.path.join(".files", session_id, "results")

    if not os.path.exists(os.path.join(results_dir, input.file_path1)):
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error: File {input.file_path1} does not exist.",
                error=True,
            )
        ]
    if not os.path.exists(os.path.join(results_dir, input.file_path2)):
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error: File {input.file_path2} does not exist.",
                error=True,
            )
        ]
    try:
        integration_folder = os.path.join(".files", session_id, "results", "integration")
        adata1 = sc.read(os.path.join(results_dir, input.file_path1))
        adata2 = sc.read(os.path.join(results_dir, input.file_path2))
        cluster_int = Cluster_Integration(
            adata1,
            adata2,
            input.file_path1.split("/")[-1].split(".")[0],
            input.file_path2.split("/")[-1].split(".")[0],
            os.path.join(".files", session_id, "results", "clustering"),
            integration_folder,
        )
        async_integration = make_async(integrate)
        result = await async_integration(cluster_int)
        results.append(
            ToolResult(
                type=ToolResultType.text,
                content=f"Integration completed successfully. Result: {result}",
            )
        )
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error during integration: {str(e)}",
                error=True,
            )
        ]
        
    try:
        async_plot_umap = make_async(plot_umap_cluster_int)
        umap_file1 = await async_plot_umap(cluster_int, adata1)
        umap_file2 = await async_plot_umap(cluster_int, adata2)
        results.append([
            ToolResult(
                type=ToolResultType.image,
                content=umap_file1,
                error=False,
                desc="UMAP plot of the first integrated dataset."
            )
        ])
        results.append([
            ToolResult(
                type=ToolResultType.image,
                content=umap_file2,
                error=False,
                desc="UMAP plot of the second integrated dataset."
            )
        ])
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error while plotting UMAP: {str(e)}",
                error=True,
            )
        ]
    try:
        async_plot_overlap = make_async(plot_overlap_umap_int)
        overlap_file = await async_plot_overlap(cluster_int, adata1, adata2, integration_folder)
        results.append(
            ToolResult(
                type=ToolResultType.image,
                content=overlap_file,
                error=False,
                desc="UMAP plot of the integrated datasets with overlap.",
            )
        )
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error while plotting overlap UMAP: {str(e)}",
                error=True,
            )
        ]
    return results
    