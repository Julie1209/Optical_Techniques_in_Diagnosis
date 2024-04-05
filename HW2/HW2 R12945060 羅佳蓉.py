# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 17:10:43 2023

@author: cdpss
"""

import numpy as np
import matplotlib.pyplot as plt

ua = 10  
#us = 0    
del_z = 0.025 
N = 10000 
x = np.zeros(N) 

for i in range(N):
    z = 0
    absorb = False
    while not absorb:
        if np.random.rand() <= ua * del_z:
            absorb = True
            x[i] = z
        else:
            z = z + del_z

# Calculate the Beer-Lambert law prediction
depths = np.arange(0, 1 + del_z, del_z)
beer_lambert_prediction = del_z*N*ua*np.exp(-ua * depths)

# Create histogram
bin_edges = np.arange(0, 1 + del_z, del_z)
plt.hist(x, bins=bin_edges, edgecolor='k', alpha=0.7)
plt.plot(depths, beer_lambert_prediction, label="Beer-Lambert")
plt.title('Bar graph is the simulated result of r.n., while the solid line is that from Beer Law')
plt.xlabel('Depth')
plt.ylabel('photons absorbed')
plt.legend()
plt.show()

