# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 20:36:28 2024

@author: simon
using Copilot base on task
implement TPS in Python
See also
https://grass.osgeo.org/
https://github.com/OSGeo/grass-addons/tree/grass8/src/vector/v.surf.tps
https://se.mathworks.com/matlabcentral/fileexchange/37576-3d-thin-plate-spline-warping-function
https://github.com/topinfrassi01/3d-deformation-algorithms/blob/main/numpy/tps3d.py
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.linalg import solve
import matplotlib.pyplot as plt

def thin_plate_spline(x, y, z, lambda_smooth=0.0):
    """
    Compute thin plate spline deformation.

    Args:
        x (np.ndarray): X coordinates of control points.
        y (np.ndarray): Y coordinates of control points.
        z (np.ndarray): Corresponding values at control points.
        lambda_smooth (float): Smoothing parameter (optional).

    Returns:
        np.ndarray: Deformed values at new points.
    """
    n = len(x)
    K = cdist(np.column_stack((x, y)), np.column_stack((x, y)), metric='euclidean')
    P = np.column_stack((np.ones(n), x, y))
    L = np.vstack((np.hstack((K, P)), np.hstack((P.T, np.zeros((3, 3))))))

    # Solve linear system
    w = solve(L + lambda_smooth * np.eye(n + 3), np.concatenate((z, np.zeros(3))))

    def tps_transform(new_x, new_y):
        r = np.sqrt((x - new_x) ** 2 + (y - new_y) ** 2)
        phi = r ** 2 * np.log(r + 1e-10)  # Avoid log(0)
        affine = np.dot(np.array([1, new_x, new_y]), w[:3])
        deformation = np.dot(phi, w[3:])
        return affine + deformation

    return tps_transform

# Example control points
control_x = np.array([0, 1, 2, 3, 4])
control_y = np.array([0, 1, 2, 3, 4])
control_z = np.array([0, 1, 4, 9, 16])

# Compute TPS transformation
tps = thin_plate_spline(control_x, control_y, control_z, lambda_smooth=0.1)

# Evaluate TPS at new points
x = np.linspace(0, 4, 5)
y = np.linspace(0, 4, 5)
z = tps(x, y)

# Plot control points and deformed surface
plt.figure(figsize=(8, 6))
plt.scatter(control_x, control_y, c='red', label='Control Points')
plt.contourf(x, y, z, levels=20, cmap='viridis', alpha=0.7)
plt.colorbar(label='Deformed Value')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Thin Plate Spline Deformation')
plt.legend()
plt.show()
