from pydantic import BaseModel
import os
import json


class ClusterRecord(BaseModel):
    cluster_file: str
    modality: str
    regions: list[str]

    def __str__(self):
        return f"ClusterToModality(session_id={self.session_id}, cluster_file={self.cluster_file}, modality_file={self.modality_file})"


def initialize_json(session_id: str):
    """
    Initialize the JSON file for storing cluster to modality mappings.
    """
    os.mkdir(os.path.join(".files", session_id))
    file_path = os.path.join(".files", session_id, "session_data.json")
    if not os.path.exists(file_path):
        with open(file_path, "w") as f:
            json.dump({}, f)

def insert_regions(session_id: str, regions: list[str]):
    """
    Insert regions into the JSON file.

    Args:
        session_id (str): Unique identifier for the session.
        regions (list[str]): List of regions to insert.
    """
    file_path = os.path.join(".files", session_id, "session_data.json")
    with open(file_path, "r+") as f:
        data = json.load(f)
        data["regions"] = regions
        f.seek(0)
        json.dump(data, f, indent=4)

def get_regions(session_id: str) -> list[str]:
    """
    Get regions from the JSON file.

    Args:
        session_id (str): Unique identifier for the session.

    Returns:
        list[str]: List of regions.
    """
    file_path = os.path.join(".files", session_id, "session_data.json")
    with open(file_path, "r") as f:
        data = json.load(f)
        return data.get("regions", [])

def insert_cluster_record(
    session_id: str, cluster_file: str, modality: str, regions: list[str]
):
    """
    Insert a record into the JSON file.

    Args:
        cluster_file (str): name of the cluster file.
        modality str: Modality of the cluster file.
        regions (list[str]): List of regions in the cluster file.
    """
    file_path = os.path.join(".files", session_id, "session_data.json")
    record = ClusterRecord(
        cluster_file=cluster_file, modality=modality, regions=regions
    )

    with open(file_path, "r+") as f:
        data = json.load(f)
        if "cluster_records" not in data:
            data["cluster_records"] = []
        data["cluster_records"].append(record.model_dump())
        f.seek(0)
        json.dump(data, f, indent=4)


def get_record(session_id: str, cluster_file: str):
    """
    Retrieve modality associated with a given session ID and cluster file.

    Args:
        session_id (str): Unique identifier for the session.
        cluster_file (str): Path to the cluster file.

    Returns:
        record (ClusterRecord | None): The found record or None if not found.
    """
    file_path = os.path.join(".files", session_id, "session_data.json")

    if not os.path.exists(file_path):
        return []

    with open(file_path, "r") as f:
        data = json.load(f)

    for record in data["cluster_records"]:
        if record["cluster_file"] == cluster_file:
            return record
    return None
