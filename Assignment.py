# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 12:59:05 2025

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.geometry import box
np.random.seed(19680806)


""" Variables """
# Tuning variables:
node_max = 100
area_size = 100
radius = 10
min_coverage = 0.9
# Initial variables:
covered_area = 0
# Dependent variables
box_bounds = box(0,0,area_size, area_size)


""" Functions """
def plotNodes(locations,node_max,area_size): # uniformly distributes the location of nodes and returns it.
    for i in range(node_max):
        node_location = np.random.uniform(0,area_size,[1,2])
        locations[i] = node_location
    return locations

def createBoundedPoints(locations,bounds,j): # this creates shapely circles bounded by area size and non-overlap.
    temp_shapes = []
    for i in range(j):
        temp_points = Point(locations[i,:]).buffer(radius)
        temp_shapes.append(temp_points)
    # Remove unions from area calculations (actually more of a merge of shapes):
    covered_area = unary_union(temp_shapes)
    # Remove everything exceeding the allowed area:
    covered_area = covered_area.intersection(bounds)
    return covered_area


""" Main Code """
for j in range(node_max): # Try for number of nodes
    if (covered_area < min_coverage) and (True): # Skip loop if the covered area is good enough
        # we can also add a condition for being well connected in this statement (replace True part with this)
        # to make sure enough area is covered and everything is connected
        temp = np.zeros([j,2])
        locations = plotNodes(temp,j,area_size)
        print("iteration",j,":")
        # Create points with radius for area calculations:
        covered_area = createBoundedPoints(locations,box_bounds,j)
        temp2 = covered_area
        covered_area = covered_area.area/(area_size**2)
        print("Area filed",(covered_area*100),"%")
# x & y might not be needed tbh:
x = locations[:,0]
y = locations[:,1]


""" Chatgpt to plot lol """
from shapely.geometry import MultiPolygon, Polygon

# Set up the plot
fig, ax = plt.subplots(figsize=(8, 8))

# Plot the geometry safely
if isinstance(temp2, Polygon):
    x1, y1 = temp2.exterior.xy
    ax.fill(x1, y1, color='lightblue', alpha=0.5, edgecolor='blue')
elif isinstance(temp2, MultiPolygon):
    for shape in temp2.geoms:
        x1, y1 = shape.exterior.xy
        ax.fill(x1, y1, color='lightblue', alpha=0.5, edgecolor='blue')

# Plot the original node locations
ax.plot(x, y, 'ro')  # red dots for node centers

ax.set_aspect('equal')
plt.title('Clipped and Unioned Coverage Area')
plt.grid(True)
plt.show()

#plt.scatter(x, y, s=radius, alpha=0.5)
#plt.show()