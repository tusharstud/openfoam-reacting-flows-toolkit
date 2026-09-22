#!/usr/bin/env python3
"""
extract_flame_profiles.py
-------------------------
Author: Tushar Sagar (Aerospace Engineering)
Description:
    Extracts and parses centerline, radial, and time-series probe data
    from OpenFOAM postProcessing directories (reactingFoam / rhoReactingFoam).
"""

import os
import glob
import numpy as np
import pandas as pd


def find_latest_time_dir(sample_base_path):
    """Finds the directory corresponding to the latest OpenFOAM simulation time."""
    if not os.path.exists(sample_base_path):
        return None
    subdirs = [d for d in os.listdir(sample_base_path) if os.path.isdir(os.path.join(sample_base_path, d))]
    try:
        subdirs_float = [(d, float(d)) for d in subdirs]
        subdirs_float.sort(key=lambda x: x[1])
        return os.path.join(sample_base_path, subdirs_float[-1][0])
    except ValueError:
        return None


def parse_openfoam_sampled_file(file_path):
    """
    Parses OpenFOAM tabulated text file containing (x, y, z) coordinates
    and scalar/vector fields. Ignores header comments starting with '#'.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    # OpenFOAM lines formatted as: coordinate field_val_1 field_val_2 ...
    data = []
    headers = []
    with open(file_path, "r") as f:
        for line in f:
            line_str = line.strip()
            if line_str.startswith("#"):
                # Check for header names
                if "x" in line_str or "T" in line_str or "coordinate" in line_str:
                    headers = line_str.replace("#", "").split()
                continue
            if line_str:
                parts = line_str.split()
                data.append([float(p) for p in parts])

    arr = np.array(data)
    if headers and len(headers) == arr.shape[1]:
        return pd.DataFrame(arr, columns=headers)
    return pd.DataFrame(arr)


def compute_flame_liftoff_height(axial_coords, temperature_array, threshold_temp=750.0):
    """
    Calculates flame lift-off height based on temperature threshold ignition criterion.
    Default threshold: 750 K (typical threshold for hydrocarbon/hydrogen autoignition marker).
    """
    ignited_idx = np.where(temperature_array >= threshold_temp)[0]
    if len(ignited_idx) > 0:
        first_ignited = ignited_idx[0]
        return axial_coords[first_ignited]
    return np.nan


if __name__ == "__main__":
    print("[INFO] OpenFOAM Flame Profile Extraction Module initialized.")
    print("[INFO] Ready to ingest postProcessing/ sampled line sets.")
