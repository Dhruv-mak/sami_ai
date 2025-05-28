import logging

import os
from SAMI.norm import Normalization
import pandas as pd
from pydantic import BaseModel, FilePath, Field
from enum import StrEnum
import time
from tools.types import ToolResultType, ToolResult
import re
from SAMI.preprocessing import csv2h5ad, pooldata
from chainlit import make_async

logger = logging.getLogger(__name__)


class RowNormMethods(StrEnum):
    """Enumeration for row normalization methods."""

    SumNorm = "SumNorm"
    MedianNorm = "MedianNorm"
    CompNorm = "CompNorm"


class TransformMethods(StrEnum):
    LogTrans2 = "LogTrans2"
    SquareRootTrans = "SquareRootTrans"
    CubeRootTrans = "CubeRootTrans"


class ScaleNormMethods(StrEnum):
    MeanCenter = "MeanCenter"
    AutoNorm = "AutoNorm"
    ParetoNorm = "ParetoNorm"
    RangeNorm = "RangeNorm"


class CompoundTypes(StrEnum):
    metabolomics = "metabolomics"
    lipidomics = "lipidomics"
    glycomics = "glycomics"
    other = "other"


class NormalizationInput(BaseModel):
    row_norm: RowNormMethods = Field(
        default=RowNormMethods.SumNorm,
        description="Row normalization method.",
    )
    trans_norm: TransformMethods = Field(
        default=TransformMethods.LogTrans2,
        description="Transformation method.",
    )
    ref_compound: str | None = Field(
        default=None,
        description="Reference compound for normalization. Only needed for CompNorm.",
    )


def validate_and_process_csv(csv_file: FilePath) -> pd.DataFrame:
    """
    Validate and process a CSV file for analysis.

    This function:
    1. Reads the CSV file
    2. Removes unnecessary columns
    3. Standardizes region column naming
    4. Validates that region data is present

    Args:
        csv_file: Path to the CSV file

    Returns:
        Processed DataFrame or None if validation fails
    """
    print(f"Validating and processing CSV file: {csv_file}")
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        return None

    # Remove unnecessary columns
    columns_to_drop = ["Unnamed: 0", "spotId", "raster", "z", "Date"]
    df = df.drop([col for col in columns_to_drop if col in df.columns], axis=1)

    # Standardize region column naming
    if "tissue_id" in df.columns:
        df = df.rename(columns={"tissue_id": "region"})
    elif "roi_named" in df.columns:
        df = df.rename(columns={"roi_named": "region"})

    # Validate region column exists
    if "region" not in df.columns:
        logger.error("Missing 'region' column in the CSV file.")
        return None

    # Remove rows with missing region data
    df = df.dropna(subset=["region"])

    return df


async def get_compound_type(
    df: pd.DataFrame,
) -> CompoundTypes:
    """Finds out the compound type based on the columns in the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the data.

    Returns:
        CompoundTypes: The compound type based on the columns in the DataFrame.
    """
    colums = df.columns[3:]
    lipidomics_pattern = re.compile(r"^[A-Z]{1,4}\(.*:.*\)$")
    if all(col.startswith("G") or col.startswith("g") for col in colums):
        return CompoundTypes.glycomics
    elif all(re.search(lipidomics_pattern, col) for col in colums):
        return CompoundTypes.lipidomics
    else:
        return CompoundTypes.metabolomics


async def normalisation_tool(
    normalisation_input: NormalizationInput,
    session_id: str,
) -> list[ToolResult]:
    """Load, validate, and normalize data from a CSV file.

    Args:
        row_norm: Row normalization method ('SumNorm', 'MedianNorm', etc.).
        trans_norm: Transformation method ('LogTrans2', etc.).
        ref_compound: Reference compound for normalization (if needed).

    Returns:
        Status message indicating success or failure.
    """
    output = ""
    if len(os.listdir(".files")) == 0:
        logger.error("No files found in the directory. Please upload a CSV files.")
        return [
            ToolResult(
                type=ToolResultType.text,
                content="Error: No files found in the directory. Please upload a CSV files.",
                error=True,
            )
        ]
    work_dir = os.path.join(".files", session_id)
    files = [x for x in os.listdir(work_dir) if x.endswith(".csv")]
    print(f"Files found: {files}")
    regions = set()
    dict_compound_type = {}
    for ind, file in enumerate(files):
        valid_start = time.time()
        df = validate_and_process_csv(os.path.join(work_dir, file))
        if normalisation_input.row_norm == RowNormMethods.CompNorm:
            if normalisation_input.ref_compound == "":
                logger.error(
                    "Reference compound is required for CompNorm normalization."
                )
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        content="Error: Reference compound is required for CompNorm normalization.",
                        error=True,
                    )
                ]
            else:
                # handling the 0s in the peak column (Copied directly from metavision3d)
                df = df[df[normalisation_input.ref_compound] != 0].reset_index(
                    drop=True
                )
        valid_end = time.time()
        logger.info(
            f"Validation and loading completed in {valid_end - valid_start:.2f} seconds for file {ind}."
        )
        if df is None:
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content="Error: Invalid CSV file format.",
                    error=True,
                )
            ]
        try:
            start = time.time()
            normalized_df = Normalization(
                df,
                first_compound_idx=3,
                rowNorm=normalisation_input.row_norm,
                transNorm=normalisation_input.trans_norm,
                c=1,
                log_base=2,
                ref_compound=normalisation_input.ref_compound,
            )

            end = time.time()
            logger.info(
                f"Normalization completed in {end - start:.2f} seconds for file {ind}."
            )

            # Save the normalized data to a new csv file
            normalized_file = f"normalized_data_{ind}.parquet"
            save_start = time.time()
            if not os.path.exists(os.path.join(work_dir, "normalized_data")):
                os.mkdir(os.path.join(work_dir, "normalized_data"))
            normalized_df.to_parquet(
                os.path.join(work_dir, "normalized_data", normalized_file), index=False
            )

            save_end = time.time()
            logger.info(
                f"Data saved to {normalized_file} in {save_end - save_start:.2f} seconds for file {ind}."
            )

            available_regions = set(df["region"].unique())
            if len(regions) == 0:
                regions.update(available_regions)
            else:
                # Check if the regions are consistent across files
                if len(regions - available_regions) != 0:
                    logger.error(
                        f"Regions in {ind} do not match previously loaded files."
                    )
                    return [
                        ToolResult(
                            type=ToolResultType.text,
                            content="Error: Regions do not match across files.",
                            error=True,
                        )
                    ]
            compound_type = await get_compound_type(df)
            dict_compound_type[normalized_file] = str(compound_type)
            data_saving_path = os.path.join(work_dir, "pooled_data")
            if not os.path.exists(data_saving_path):
                os.mkdir(data_saving_path)
            csv2h5ad(compound_type, df.copy(), data_saving_path, split=True)

        except Exception as e:
            logger.error(f"Error normalizing data: {e}")
            return [
                ToolResult(
                    type=ToolResultType.text,
                    content=f"Error: {e} while processing file {ind}.",
                    error=True,
                )
            ]
    try:
        pooldata(
            dict_compound_type,
            os.path.join(work_dir, "normalized_data"),
            os.path.join(work_dir, "pooled_data"),
            split=True,
        )
    except Exception as e:
        logger.error(f"Error pooling data: {e}")
        return [
            ToolResult(
                type=ToolResultType.text,
                content=f"Error: {e} while pooling data. Probably the pixel count is not the same",
                error=True,
            )
        ]

    return [
        ToolResult(
            type=ToolResultType.text,
            content=f"Normalization completed successfully. Files are successfully saved. The regions are: {', '.join(regions)}",
            error=False,
        )
    ]
