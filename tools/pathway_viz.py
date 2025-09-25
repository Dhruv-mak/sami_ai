from pydantic import BaseModel, Field
from tools.types import ToolResultType, ToolResult
from chainlit import make_async
import os
from SAMI.pathway import Pathway, findpathway, pathway_network, plot_bar, plot_dot


class PathwayVizInput(BaseModel):
    """
    Input model for pathway visualization.
    """

    cluster: int = Field(
        ddescription="Cluster number for pathway visualization.",
    )
    pathaway_file: str = Field(
        description="Pathway csv path for visualization.",
    )
    top_n_dot: int = Field(
        default=20,
        description="Number of top pathways to visualize in the dot plot.",
    )
    top_n_network: int = Field(
        default=10,
        description="Number of top pathways to visualize in the network plot.",
    )
    top_n_bar: int = Field(
        default=5,
        description="Number of top pathways to visualize in the bar plot.",
    )


async def pathway_viz_tool(input: PathwayVizInput, session_id: str) -> list[ToolResult]:
    pathway_folder = os.path.join(".files", session_id, "pathways")
    if not os.path.exists(pathway_folder):
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Pathway folder does not exist: {pathway_folder}",
                error=True,
            )
        ]
    if not os.path.exists(os.path.join(pathway_folder,input.pathaway_file)):
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Pathway file does not exist: {input.pathaway_file}",
                error=True,
            )
        ]
    result = []
    sep_file_list = input.pathaway_file.split("_")
    region = sep_file_list[1]
    modality = "_".join(sep_file_list[2:]).replace(".csv", "")
    pathway = Pathway(
        region=region,
        modality=modality,
        pathways_folder=pathway_folder,
        makers_folder=os.path.join(".files", session_id, "markers"),
    )
    try:
        pathway_network_async = make_async(pathway_network)
        network_image = await pathway_network_async(
            pathway,
            cluster=input.cluster,
            top=input.top_n_network,
            size=800,
            show=False,
        )
        result.append(
            ToolResult(
                type=ToolResultType.image,
                content=network_image,
                description="Network visualization of pathways.",
            )
        )
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error during pathway visualization: {str(e)}",
                error=True,
            )
        ]
    try:
        plot_bar_async = make_async(plot_bar)
        bar_plot_image = await plot_bar_async(
            pathway,
            cluster=input.cluster,
            top=input.top_n_bar
        )
        result.append(
            ToolResult(
                type=ToolResultType.image,
                content=bar_plot_image,
                description="Bar plot of pathaway directions.",
            )
        )
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Pathaway Direction plot: {str(e)}",
                error=True,
            )
        ]
    try:
        dot_plot_async = make_async(plot_dot)
        dot_plot_image = await dot_plot_async(
            pathway,
            cluster=input.cluster,
            scale=30,
            height=10,
            top=input.top_n_dot,
        )
        result.append(
            ToolResult(
                type=ToolResultType.image,
                content=dot_plot_image,
                description="Dot plot of pathways.",
            )
        )
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error during Dot plot: {str(e)}",
                error=True,
            )
        ]
    return result