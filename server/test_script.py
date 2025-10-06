# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 20:23:58 2025

@author: simon and copilot
"""
#%% test rectangle
import asyncio
from analyze_core import analyze

params = {
    "param1": 0.2,  # width [m]
    "param2": 0.4,  # height [m]
    "mesh_size": 0.005
}

result = asyncio.run(analyze(params))
print(result["result"]["A"])
print(result["result"]["warping"][:5])  # Show first 5 warping points
