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
import cvxpy as cp
np.random.seed(5885221)


""" Variables """
# Tuning variables:
node_max = 200
max_iter = 500
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

def convex_RG(W_ij_list):
    n = len(W_ij_list)
    P = cp.Variable((n,n), nonneg = True)
    sum_PW = np.zeros((n,n))
    constraints = []
    
    for i in range(n):
        constraints.append(cp.sum(P[i]) == 1)
        for j in range(n):
            sum_PW =  sum_PW + (P[i,j]*W_ij_list[i][j])
    W_bar = (1/n)*sum_PW
    constraints.append(W_bar == (1/n)*sum_PW)
    lambda_max = cp.lambda_max(W_bar-(1/n)*np.ones((n,n)))

    prob = cp.Problem(cp.Minimize(lambda_max), constraints)
    prob.solve(solver=cp.CVXOPT)
    print("Solver status:", prob.status)
    return P.value

def convex_PDMM(x,z,n):
    A = cp.Variable((n,n,num_connections), nonneg = True) # each element needs to be a matrix with dimension that fits: A(i,j)*x(i)
    b = cp.Variable((n,n), nonneg = True)
    constraints = []
    for i in range(n):
        for j in range(n): # needs to limit to only connected nodes
            constraints.append(A(i,j)*x(i) + A(j,i)*x(j) == b(i,j))
            cost = f(i) + sum(z(i,j).T * A(i,j)*x(i) + (c/2)*cp.power(2,(cp.norm(A(i,j)*x(i)-(1/2)*b(i,j)))))
            # what is f(i) and c?
    prob = cp.Problem(cp.Minimize(cost), constraints)
    prob.solve(solver=cp.CVXOPT)
    print("Solver status:", prob.status)

def randomized_gossip(A):
    n = A.shape[1]
    W_ij_list = []
    I_mtrx = np.eye(n)
    for i in range(n):
        e_i = np.zeros((n,1))
        e_i[i] = 1
        temp = []
        for j in range(n):
            if A[i,j] == 1 and i!=j:
                e_j = np.zeros((n,1))
                e_j[j] = 1
                W_ij = I_mtrx-(1/2)*(e_i-e_j)*(e_i-e_j).T
                temp.append(W_ij)
            else:
                temp.append(np.zeros((n,n)))
        W_ij_list.append(temp)
    P = convex_RG(W_ij_list)
    return P
"""""
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

def PDMM(k,x): # k is the iteration
    n = A.shape[1]
    I_mtrx = np.eye(n)
    z = np.zeros(2*m,2*m) # m?
    for i in range(n): # total number of nodes
        x_i = convex_PDMM(j) # to find the argmin we need a convex solver! j is the connected nodes
        for j in range(n):
            if A[i,j] == 1 and i!=j:
                y(i,j) = z(i,j) + 2*c*(A_ij*x(i)-(1/2)*b_ij) # A_ij & b_ij?
    for i in range(n):
        for j in range(n):
            if A[i,j] == 1 and i!=j:
                print("sending!")
                # transmit variables somehow
    for i in range(n):
        for j in range(n):
            if A[i,j] == 1 and i!=j:
                # think we update the z value of node j with the y value of node i
                z(j,i) = y(i,j) #slides say k+1 but since its at the end we can just overwrite it i think
    return x,y,z
"""""
def find_partner(P):
    n = P.shape[1]
    rng = np.random.default_rng() # this needs to be a random state
    i = rng.choice(n,1)
    P[i] = np.maximum(P[i],0)
    P[i] /= np.sum(P[i])
    j = rng.choice(n,1,p=P[i][0])
    print(i,j)
    return i,j

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

## Randomized gossip:

num_nodes = j -1
num_edges = int(sum(sum(np.tril(adjacency_matrix)-np.eye((num_nodes)))))

x0 = np.random.uniform(low=0, high=50.0, size=num_nodes)
x_iterating = x0
print(x0)
x_list = []
P = randomized_gossip(adjacency_matrix) # weighted values for x[k]?

print(P)
for k in range(max_iter):
    #print(adjacency_matrix.shape)
    i,j = find_partner(P) # here we need to use find_partner to find the specific index of node k to place the average in!
    #x_iterating[connecting_node],_ = randomized_gossip(adjacency_matrix,x_iterating[connecting_node],num_nodes) # weighted values for x[connecting node]?
    x_mean = (1/2)*(x_iterating[i] + x_iterating[j])
    print(x_iterating)
    x_iterating[i] = x_mean
    x_iterating[j] = x_mean
    # these 2 nodes now both need to set their values to the average: 1/2(x_k+x_connected) (I think?)
    ### but do they use the same weights or other ones? then is it randomized? or is only the chosen node randomized?
    # something like this?
    #x_iterating[k],x_iterating[connecting_node] = 1/2*(x_iterating[k] + x_iterating[connecting_node])
    x_list.append(x_iterating.copy()) # x_list can keep track of all values
# x & y might not be needed tbh:
x = locations[:,0]
y = locations[:,1]


## PDMM:
n = num_nodes
m = num_edges
neighbors_list = []
for i in range(n):
    neighbors_i = []
    for j in range(n):
        if (adjacency_matrix[i,j] == 1):
            neighbors_i.append(j)
    neighbors_list.append(neighbors_i)

A_list = []
B_list = []
d_list = []

c = 5
for i in range(n):
    A_list_j = []
    B_list_j = []
    d_list_j = []
    for j in range(n):
        if (adjacency_matrix[i,j] == 1 and i>j):
            # we need to add the C matrix here still!
            A_list_j.append(-1)
            B_value = np.zeros((m,1))
            B_list_j.append(B_value)
            d_value = (1/2)*np.vstack((B_value.T,B_value.T))
            d_list_j.append(d_value)
        elif (adjacency_matrix[i,j] == 1 and i<j):
            A_list_j.append(1)
            B_value = np.zeros((m,1))
            B_list_j.append(B_value)
            d_value = (1/2)*np.vstack((B_value,B_value))
            d_list_j.append(d_value)
        else:
            A_list_j.append(0)
            B_list_j.append(0)
            d_list_j.append(np.zeros((1,2*m)))
    
    A_list.append(A_list_j)
    B_list.append(B_list_j)
    d_list.append(d_list_j)



print("test")
print(A_list[i][0])
a = np.zeros((n,1))
x_PDMM = np.zeros((n,max_iter+1)) # put x0 as first iteration
x_PDMM[:,0] = x0
m = 1
for i in range(n):
    a[i] = np.average(adjacency_matrix[i,:]*x0)
    num_neighbors = int(sum(adjacency_matrix[i,:])-1)
    z = np.zeros((m,num_neighbors))
    summation = np.zeros((num_neighbors,m))
    for time,k in enumerate(range(max_iter)):
        print(time)
        for index,j in enumerate(neighbors_list[i]):
            summation[index] = (A_list[i][j])*z[index]
        x_PDMM[i,k+1] = (a[i] - np.sum(summation)/(1+c*num_neighbors))
        for index,j in enumerate(neighbors_list[i]):
                y = z[:,index] + 2*c*(A_list[i][j]*x[i,k+1])
                z[:,index] = y

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

x_array = np.array(x_list)  # shape: (iterations+1, num_nodes)

plt.figure(figsize=(10, 6))
num_nodes = x_array.shape[1]

for i in range(num_nodes):
    plt.plot(x_array[:, i])

plt.xlabel('Iteration $k$')
plt.ylabel('Value $x_i[k]$')
plt.title('Randomized Gossip Convergence Over Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Compute the average value of the initial state (consensus target)
x_avg = np.mean(x_array[0])  # constant target value
errors = []

# Compute squared consensus error at each iteration
for x_k in x_array:
    error = np.linalg.norm(x_k - x_avg)**2
    errors.append(error)

# Plot on semilogy scale (like your reference plot)
plt.figure(figsize=(10, 6))
plt.semilogy(errors, label='Randomized Gossip', color='darkorange', linewidth=2)

plt.xlabel('Iteration $k$')
plt.ylabel(r'$||x(k) - x_{\mathrm{avg}}*1||^2$')
plt.title('Consensus Error Over Time')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.show()
#plt.scatter(x, y, s=radius, alpha=0.5)
#plt.show()