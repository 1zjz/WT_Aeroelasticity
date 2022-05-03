# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 22:48:05 2022

@author: carlo
"""

import numpy as np
import matplotlib.pyplot as plt
# import functions as func # File containing user-defined functions

# Import data
path = "./tangential_load_yaw_30.txt"
values = np.loadtxt(path, skiprows=1, delimiter=',')/25

azimuthal = np.radians(np.loadtxt('azimuthal_discretization.txt', skiprows=1, delimiter=','))
radial = np.loadtxt('radial_discretization.txt', skiprows=1, delimiter=',')

# Using linspace so that the endpoint of 360 is included
# actual = np.radians(np.linspace(0, 360, 20))
# expected = np.arange(0, 70, 10)
 
r, theta = np.meshgrid(radial, azimuthal)
# values = np.random.random((radial.size, azimuthal.size))
 
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
levels = [0, 0.04, 0.08, 0.12, 0.16, 0.2] # 0.12, 0.16, 0.20] #, 2.5, 3.0, 3.5]
CS  = ax.contourf(theta, r, np.transpose(values), levels)
# CS2 = ax.contour(theta, r, np.transpose(values), colors='k')
# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('Tangential Induction')
# cbar.add_lines(CS2)
plt.show()