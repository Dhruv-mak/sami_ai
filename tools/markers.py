from pydantic import BaseModel, Field
from tools.types import ToolResultType, ToolResult
import scanpy as sc
from db import get_record
import logging
import numpy as np
from SAMI.markers import Markers, findmarkers, circular_tree
import uuid
from chainlit import make_async
from pandas.core.common import flatten
import os

logger = logging.getLogger(__name__)


def add_preva_abund(adata):
    abundance = np.sum(adata.X != 0, axis=1) / adata.shape[1]
    prevalence = np.sum(adata.X != 0, axis=0) / adata.shape[0]
    adata.obs["abundance"] = abundance
    adata.var["prevalence"] = prevalence
    return adata


def adata_filter(adata, abundance=0.1, prevalence=0.1):
    adata_temp = add_preva_abund(adata)
    adata_filter = adata_temp[
        adata_temp.obs["abundance"] >= abundance,
        adata_temp.var["prevalence"] >= prevalence,
    ]
    return adata_filter


def validate_entry(entry, highest):
    try:
        if entry != "":
            entry = int(entry)
            if entry <= highest:
                return True
            else:
                return False

    except:
        return False


def get_dict_clustermap(mode, lst_clusters, highest_cluster):
    matrix = []
    if mode == "single":
        for cluster in lst_clusters:
            rest_clusters = [str(c) for c in lst_clusters if c != cluster]
            rest_clusters = ",".join(rest_clusters)
            matrix.append([str(cluster), rest_clusters])
    else:
        for cluster in lst_clusters:
            matrix.append([str(cluster), str(cluster)])
    dict_clustermap = {}
    left_list = []
    right_list = []
    for item in matrix:

        logger.info(f"item: {item}")
        left = item[0].split(",")
        right = item[1].split(",")

        left = [item for item in left if validate_entry(item, highest_cluster)]
        right = [item for item in right if validate_entry(item, highest_cluster)]

        left = list(flatten(left))
        right = list(flatten(right))

        # if either left or right is empty we are skipping that row
        if min(len(left), len(right)) == 0:
            continue
        left_list.append(left)
        right_list.append(right)

    dict_clustermap["left"] = left_list
    dict_clustermap["right"] = right_list
    return dict_clustermap


class MarkersInput(BaseModel):
    """
    Input model for the markers tool.
    """

    file1: str = Field(
        description="Name of the .h5ad file to find markers for (case/primary dataset).",
    )
    file2: str | None = Field(
        default=None,
        description="Name of the second .h5ad file for comparison (control/reference dataset). If None, performs single-dataset marker analysis.",
    )
    abundance: float = Field(
        description="Minimum abundance threshold for filtering features. Features must have this minimum expression level to be considered.",
    )
    adj_pval_cutoff: float = Field(
        default=0.05,
        description="Adjusted p-value cutoff for statistical significance. Only markers with adjusted p-value <= this threshold will be retained.",
    )
    prevalance: float = Field(
        default=0.1,
        description="Minimum prevalence threshold (fraction of cells expressing the feature). Features must be expressed in at least this fraction of cells (pct1 > 0.1 and pct2 > 0.1).",
    )
    top_n_molecules: int = Field(
        default=50,
        description="Maximum number of top marker molecules to return per cluster comparison, ranked by absolute average log2 fold change.",
    )


async def markers_tool(input: MarkersInput, session_id: str) -> list[ToolResult]:
    """
    Find marker molecules for single-cell analysis using SAMI's marker detection method.

    This tool performs statistical analysis to identify differentially expressed features
    between clusters or datasets. It can operate in two modes:
    - Single dataset: Compare clusters within one dataset
    - Dual dataset: Compare corresponding clusters between two datasets (case vs control)

    The analysis includes:
    - Data filtering based on abundance and prevalence thresholds
    - Statistical testing using rank-sum tests
    - Multiple testing correction for adjusted p-values
    - Ranking by absolute log2 fold change
    - Export of results as CSV files

    Args:
        input (MarkersInput): Input parameters containing:
            - file1: Primary .h5ad file for analysis
            - file2: Optional second .h5ad file for comparison
            - abundance: Minimum expression threshold for feature filtering
            - adj_pval_cutoff: Adjusted p-value threshold for significance
            - prevalance: Minimum fraction of cells expressing the feature
            - top_n_molecules: Maximum number of top markers to return per comparison
        session_id (str): Unique identifier for the session to organize results

    Returns:
        list[ToolResult]: List of results containing:
            - Success/error messages for the marker analysis
            - File paths to generated CSV files with marker results
            - Statistical summary of identified markers

    The function generates two output files:
        - {reference_name}_pvalue.csv: Complete statistical results for all features
        - {reference_name}_marker.csv: Filtered top markers meeting significance criteria
    """
    logger.info(f"Starting marker analysis for session {session_id}")
    result = []
    try:
        adata1 = sc.read(f".files/{session_id}/results/{input.file1}")
        logger.info(f"Loaded file1: {input.file1}")
        lst1_regions = get_record(session_id, input.file1)["regions"]
        logger.info(f"Expected regions for file1: {lst1_regions}")
        if list(adata1.obs["region"].unique()) != lst1_regions:
            logger.error(f"Regions in {input.file1} do not match the expected regions.")
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content=f"Error: Regions in {input.file1} do not match the expected regions.",
                    error=True,
                )
            ]
        lst1_regions = adata1.obs["region"].unique().tolist()
        adata1 = adata1[adata1.obs["region"].isin(lst1_regions)]
        highest_cluster = np.sort(
            [int(item) for item in adata1.obs["leiden"].unique()]
        )[-1]
        logger.info(f"Regions in {input.file1}: {lst1_regions}")
        adata_filtered1 = adata_filter(adata1, input.abundance, input.prevalance)
        mode = "single" if input.file2 is None or input.file2 == "" else "dual"
        if input.file2 is not None and input.file2 != "":
            adata2 = sc.read(f".files/{session_id}/results/{input.file2}")
            lst2_regions = get_record(session_id, input.file2)["regions"]
            if list(adata2.obs["region"].unique()) != lst2_regions:
                logger.error(
                    f"Regions in {input.file2} do not match the expected regions."
                )
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        content=f"Error: Regions in {input.file2} do not match the expected regions.",
                        error=True,
                    )
                ]
            adata2 = adata2[adata2.obs["region"].isin(lst2_regions)]
            highest_cluster = np.sort(
                [int(item) for item in adata1.obs["leiden"].unique()]
            )[-1]
            logger.info(f"Regions in {input.file2}: {lst2_regions}")
            adata_filtered2 = adata_filter(adata2, input.abundance, input.prevalance)
        else:
            adata_filtered2 = None
            logger.info("Single dataset mode: No second file provided for comparison.")

        dict_clustermap = get_dict_clustermap(
            mode, list(adata1.obs["leiden"].unique()), highest_cluster
        )
        file_name = f"{uuid.uuid4().hex[:8]}"
        markers = Markers(file_name)

        try:
            async_findmarkers = make_async(findmarkers)
            pval_file, markers_file = await async_findmarkers(
                markers,
                dict_clustermap,
                adata_filtered1,
                input.adj_pval_cutoff,
                input.top_n_molecules,
                adata_filtered2,
                os.path.join(".files", session_id, "markers"),
            )
            result.append(
                ToolResult(
                    type=ToolResultType.text,
                    content=f"Marker analysis completed successfully. pvalues file {pval_file} and markers file {markers_file} generated.",
                    error=False,
                )
            )
        except Exception as e:
            logger.error(f"Error in findmarkers: {e}")
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content=f"Error: during marker analysis {e}",
                    error=True,
                )
            ]
        try:
            async_circular_tree = make_async(circular_tree)
            circular_tree_file = await async_circular_tree(
                markers,
                clusters=None,
                top_n=input.top_n_molecules,
                show=False,
                workdir=os.path.join(".files", session_id, "markers"),
            )
            result.append(
                ToolResult(
                    type=ToolResultType.image,
                    content=circular_tree_file,
                    error=False,
                    desc="Circular tree plot of marker molecules.",
                )
            )
        except Exception as e:
            logger.error(f"Error in circular_tree: {e}")
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content=f"Error: during circular tree generation {e}",
                    error=True,
                )
            ]

    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error: during marker analysis {e}",
                error=True,
            )
        ]
    return result
