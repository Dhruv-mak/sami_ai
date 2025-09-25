from pydantic import BaseModel, Field
from enum import StrEnum
from tools.types import ToolResultType, ToolResult
import os
from SAMI.pathway import Pathway, findpathway
from chainlit import make_async

class PathwayModality(StrEnum):
    """Enumeration for pathway analysis modalities."""

    metabolomics = "metabolomics"
    lipidomics = "lipidomics"
    glycomics = "glycomics"

class PathwayInput(BaseModel):
    """Input model for pathway analysis."""

    marker_file: str = Field(
        description="filename of csv marker file containing markker data.",
    )
    modality: PathwayModality = Field(
        description="Modality of the pathway analysis."
    )


async def pathway_tool(input: PathwayInput, session_id: str) -> list[ToolResult]:
    """
    Perform pathway analysis based on the provided input.

    Args:
        input (PathwayInput): Input data for pathway analysis.

    Returns:
        str: Result of the pathway analysis.
    """
    markers_folder = os.path.join(".files", session_id, "markers")
    pathway_folder = os.path.join(".files", session_id, "pathways")
    if not os.path.exists(pathway_folder):
        os.makedirs(pathway_folder)
    if not os.path.exists(os.path.join(markers_folder, input.marker_file)):
        return [
            ToolResult(
                type=ToolResultType.text,
                message=f"Markers file does not exist, it should have been created after : {markers_folder}",
            )
        ]
    pathway = Pathway(
        region=input.marker_file.replace("_marker.csv", ""),
        modality=input.modality.value,
        pathways_folder=pathway_folder,
        makers_folder=markers_folder
    )
    try:
        findpathway_async = make_async(findpathway)
        pathway_file = await findpathway_async(pathway)
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Pathway analysis completed successfully. Pathways saved to {pathway_file}",
                error=False,
            )
        ]
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error during pathway analysis: {e}",
                error=True,
            )
        ]