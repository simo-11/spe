# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 20:19:27 2025

@author: simon and copilot
"""
import os
import json
import logging
from hashlib import sha256
from datetime import datetime
import asyncio
import numpy as np
from sectionproperties.pre import geometry
from sectionproperties.analysis import Section

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CACHE_DIR = "cache"
os.makedirs(CACHE_DIR, exist_ok=True)

def make_key(params: dict) -> str:
    """Generate a unique hash key based on sorted input parameters."""
    sorted_items = sorted(params.items())
    key_str = json.dumps(sorted_items)
    return sha256(key_str.encode()).hexdigest()

def get_cache_path(key: str) -> str:
    return os.path.join(CACHE_DIR, f"{key}.json")

def load_cached_result(key: str) -> dict | None:
    path = get_cache_path(key)
    if os.path.exists(path):
        with open(path, "r") as f:
            return json.load(f)
    return None

def save_result_to_cache(key: str, result: dict):
    path = get_cache_path(key)
    with open(path, "w") as f:
        json.dump({
            "result": result,
            "timestamp": datetime.now().isoformat()
        }, f, indent=2)

async def perform_analysis(params: dict) -> dict:
    """Run sectionproperties-based analysis and extract geometric and warping data."""
    logger.info(f"Starting sectionproperties analysis for: {params}")
    width = params["param1"]
    height = params["param2"]
    mesh_size = params.get("mesh_size", 0.01)

    # Create rectangular geometry
    geom = geometry.RectangularGeometry(width=width, depth=height)
    geom.create_mesh(mesh_size=mesh_size)
    section = Section(geometry=geom)

    # Run analysis
    section.calculate_geometric_properties()
    section.calculate_warping()

    # Extract warping displacement and coordinates
    disp = section.get_warping_displacement()
    coords = section.mesh.points
    warping_data = [
        {"x": float(x), "y": float(y), "w": float(w)}
        for (x, y), w in zip(coords, disp)
    ]

    result = {
        "A": float(section.get_area()),
        "Ixx": float(section.get_ixx()),
        "Iyy": float(section.get_iyy()),
        "J": float(section.get_j()),
        "shear_center": float(section.get_shear_center()),
        "warping": warping_data
    }

    logger.info(f"Analysis completed with {len(warping_data)} warping points")
    return result

async def analyze(params: dict) -> dict:
    """Main entry point: check cache, perform analysis if needed."""
    key = make_key(params)
    cached = load_cached_result(key)

    if cached:
        logger.info(f"Returning cached result for key: {key}")
        return {"status": "cached", "result": cached["result"]}

    result = await perform_analysis(params)
    save_result_to_cache(key, result)
    logger.info(f"Result cached under key: {key}")
    return {"status": "computed", "result": result}
