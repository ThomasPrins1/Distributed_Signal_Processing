# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 12:59:05 2025

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
np.random.seed(19680806)


#variables:
node_max = 5
max_iterative_nodes = 100
radius = 10
min_coverage = 0.9
covered_area = 0

# Functions:
def plotNodes(locations,node_max):
    for i in range(node_max):
        node_location = np.random.uniform(0,100,[1,2])
        locations[i] = node_location
    return locations

# Code:
for j in range(max_iterative_nodes): # This should allow for a minimal coverage
    if covered_area < min_coverage:
        locations_temp = np.zeros([j,2])
        locations = plotNodes(locations_temp,j)

        covered_area = ((j*np.pi*(radius**2)))/(100*100) # does not account for overlap yet!
        print(covered_area)
x = locations[:,0]
y = locations[:,1]




plt.scatter(x, y, s=radius, alpha=0.5)
plt.show()