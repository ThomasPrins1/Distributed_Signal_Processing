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
np.random.seed(5885221)


""" Variables """
# Tuning variables:
node_max = 200
area_size = 100
radius = 20
min_coverage = 0.85
# Initial variables:
covered_area = 0
# Dependent variables
box_bounds = box(0,0,area_size, area_size)
connection_type = "Proper_Connected"
connected_graph = False

""" Functions """
def plotNodes(locations,node_max,area_size): # uniformly distributes the location of nodes and returns it.
    for i in range(node_max):
        node_location = np.random.uniform(0,area_size,[1,2])
        locations[i] = node_location
    return locations

def createBoundedPoints(locations,bounds,k): # this creates shapely circles bounded by area size and non-overlap.
    temp_shapes = []
    adjacency_matrix = np.zeros((k,k))
    for i in range(k):
        temp_points = Point(locations[i,:]).buffer(radius)
        temp_shapes.append(temp_points)
        for j in range(k):
            if temp_points.contains(Point(locations[j,:])):
                adjacency_matrix[i,j] = 1
            else:
                adjacency_matrix[i,j] = 0
    
    # Remove unions from area calculations (actually more of a merge of shapes):
    covered_area = unary_union(temp_shapes)
    # Remove everything exceeding the allowed area:
    covered_area = covered_area.intersection(bounds)
    return covered_area,adjacency_matrix

def randomized_gossip(A,input,j):
    #for k in range(j)
        #search index of A matrix to see whos connected
    temp_weights = np.random.uniform(low=0, high=1.0, size=(j,j))
    randomized_weights = np.tril(temp_weights) + np.tril(temp_weights, -1).T
    weight_matrix = np.multiply(randomized_weights,A)
    # fix the weight matrix to sum for rows and collumns to 1
    weight_matrix[0,:] /= np.sum(weight_matrix[0,:])
    weight_matrix[:,0] = weight_matrix[0,:]
    for i in range(1,j):
        toBeNormalised_matrix = ((weight_matrix[i:j,i]/np.sum(weight_matrix[i:j,i]))*(1-np.sum(weight_matrix[0:i-1,i])))
        weight_matrix[i-1:-1,i] = toBeNormalised_matrix
        weight_matrix[i,i-1:-1] = toBeNormalised_matrix.T
    out = weight_matrix@input
    return out,weight_matrix

def find_partner(node, adjacent):
    # I guess we need to find a partner here to combine with node k for an iteration?
    # get the indexes of adjacent for row k to find which it can connect to.
    # then choose one of those indexes by random, remember that the amount of connections is kinda random
    # so this might actually be harder then it sounds lol, maybe keep in mind that node k doesnt send to node k also!
    # return index to connect to!
    print("nothing yet")

""" Main Code """

for j in range(1,node_max): # Try for number of nodes
    print("iteration",j,":")
    if (covered_area < min_coverage) or connected_graph == False: # Skip loop if the covered area is good enough
        # we can also add a condition for being well connected in this statement (replace True part with this)
        # to make sure enough area is covered and everything is connected
        temp = np.zeros([j,2])
        locations = plotNodes(temp,j,area_size)
        np_locations = np.array(locations)
        # Create points with radius for area calculations:
        covered_area,adjacency_matrix = createBoundedPoints(locations,box_bounds,j)
        temp2 = covered_area
        covered_area = covered_area.area/(area_size**2)
        print("Area filed",(covered_area*100),"%")
        
        adjacency_power = []
        for i in range(j):
            adjacency_power.append(np.linalg.matrix_power(adjacency_matrix,i+1))
            #print("non-zeros:",np.sum(np.count_nonzero(np.sum(adjacency_power,axis=0),axis=0)), (i+1)**2)
            if np.sum(adjacency_matrix[i,:]) > 1 and connection_type == "Connected": # always a self loop! so needs to be bigger then 1
                connected_graph = True
            #elif np.count_nonzero(np.sum(adjacency_power,axis=0)) >= ((i+1)**2)  and connection_type == "Proper_Connected":
            elif np.count_nonzero(adjacency_power[i]) >= ((i+1)**2)  and connection_type == "Proper_Connected":
                connected_graph = True
            else:
                connected_graph = False
                print("i needed to connect", i)
                break
            print("i needed to connect", i)
        print("Graph is now connected:",connected_graph)
    else:
        break

num_nodes = j -1
x0_temp = np.random.uniform(low=0, high=20.0, size=(num_nodes,num_nodes))
x0 = x0_temp*np.eye(num_nodes)
x_iterating = x0
print(x0)
x_list = []
for k in range(num_nodes):
    print(adjacency_matrix.shape)
    connecting_node = find_partner(k,adjacency_matrix) # here we need to use find_partner to find the specific index of node k to place the average in!
    x_iterating[k],_ = randomized_gossip(adjacency_matrix,x_iterating[k],num_nodes) # weighted values for x[k]?
    #x_iterating[connecting_node],_ = randomized_gossip(adjacency_matrix,x_iterating[connecting_node],num_nodes) # weighted values for x[connecting node]?
    
    # these 2 nodes now both need to set their values to the average: 1/2(x_k+x_connected) (I think?)
    ### but do they use the same weights or other ones? then is it randomized? or is only the chosen node randomized?
    # something like this?
    #x_iterating[k],x_iterating[connecting_node] = 1/2*(x_iterating[k] + x_iterating[connecting_node])
    x_list.append(x_iterating) # x_list can keep track of all values
# x & y might not be needed tbh:
x = locations[:,0]
y = locations[:,1]
#adjacency = 


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

for i in range(len(adjacency_matrix)):
    for j in range(i + 1, len(adjacency_matrix)):
        if adjacency_matrix[i][j]:
            x_values = [locations[i][0], locations[j][0]]
            y_values = [locations[i][1], locations[j][1]]
            plt.plot(x_values, y_values, 'k-')  # 'k-' denotes a black line
ax.set_aspect('equal')
plt.title('Clipped and Unioned Coverage Area')
plt.grid(True)
plt.show()

#plt.scatter(x, y, s=radius, alpha=0.5)
#plt.show()