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
import math
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
    K = cdist(np.column_stack((x, y)), 
              np.column_stack((x, y)), 
              metric='euclidean')
    P = np.column_stack((np.ones(n), x, y))
    L = np.vstack((np.hstack((K, P)), 
                   np.hstack((P.T, np.zeros((3, 3))))))
    # Solve linear system
    a=L + lambda_smooth * np.eye(n + 3)
    b=np.concatenate((z, np.zeros(3)))
    w = solve(a,b)

    def tps_transform(new_x, new_y):
        r = np.sqrt((x - new_x) ** 2 + (y - new_y) ** 2)
        phi = r ** 2 * np.log(r + 1e-10)  # Avoid log(0)
        affine = np.dot(np.array([1, new_x, new_y]), w[:3])
        deformation = np.dot(phi, w[3:])
        return affine + deformation

    return tps_transform

def countour_plot():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.scatter(c_x, c_y, c_z,c='red', label='Control Points')
    if isinstance(z,np.ndarray) and z.ndim==2:
        plt.contourf(xls, yls, z, levels=50, 
                     cmap=plt.cm.seismic, alpha=alpha)
        plt.colorbar(label='Z Value')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('Thin Plate Spline Deformation')
    plt.legend()
    plt.show()

def rect_psi(x,y):
    """ 
    en.wikiversity.org/wiki/Warping_functions    
    """
    if rect_a<rect_b:
        raise(f"rect_a({rect_a} must be >= rect_b\
 {rect_b}")
    cx=x-rect_x0
    cy=y-rect_y0
    z=cx*cy
    for n in range(rect_nMax):
        kn=(2*n+1)*math.pi/(2*rect_a);
        sv=((math.pow(-1,n)*math.sin(kn*cx)*math.sinh(kn*cy))/
            (math.pow(2*n+1,3)*math.cosh(kn*rect_b)))
        z=z-rect_sm*sv
    return z

def linear_x_y(x,y):
    return x_mul*x+y_mul*y

def get_z_values(f,xa,ya):
    """
    Parameters
    ----------
    f : function to be called
    xa : np.array containing x-values
    ya : np.array containing y-values

    Returns
    -------
    za : np.array containing z-values
        one value for each x,y pair
    """
    za = np.zeros((xa.size,))
    for i in range(xa.size):
        za[i]=f(xa[i],ya[i])
    return za        

def get_z_values_matrix(f,xa,ya):
    """
    Parameters
    ----------
    f : function to be called
    xa : np.array containing x-values
    ya : np.array containing y-values

    Returns
    -------
    za : np.array containing z-values
    """
    za = np.zeros((ya.size,xa.size))
    for i in range(ya.size):
        for j in range(xa.size):
            try:
                za[i,j]=f(xa[j],ya[i])
            except IndexError as e:
                print(f"i={i}, j={j} got {e}")
    return za        


alpha=0.3
# multipliers for _x_y functions
x_mul,y_mul=(1,1)
#f=linear_x_y
f=rect_psi
x_max,y_max=(10,10)
rect_width=10
rect_height=10
rect_a=rect_width/2
rect_b=rect_height/2
rect_nMax=3
rect_x0=rect_a
rect_y0=rect_b
rect_sm=32*math.pow(rect_a,2)/(math.pow(math.pi,3))
x_count=51
y_count=51
plot_analytic=False

if plot_analytic:     
    xls=np.linspace(0, x_max, x_count,dtype=float)
    yls=np.linspace(0, y_max, y_count,dtype=float)
    z=get_z_values_matrix(f,xls,yls)
    countour_plot()
# Compute TPS transformation
#x_mesh,y_mesh=np.meshgrid(xls,yls,indexing='ij')
#x=np.concatenate(x_mesh[:])
#y=np.concatenate(y_mesh[:])
# Example control points
c_x = np.array([0, 1, 1, 2, 4],dtype=float)
c_y = np.array([0, 1, 2, 2, 4],dtype=float)
c_z = get_z_values(f,c_x,c_y)
tps = thin_plate_spline(c_x, c_y, c_z)
# Evaluate TPS at new points
x = np.linspace(0, 4, 5)
y = np.linspace(0, 4, 5)
z = tps(x, y)
countour_plot()
plt.show()