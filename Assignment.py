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
max_iter_RG = 4000
max_iter_PDMM = 4000
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

def find_partner(P):
    n = P.shape[1]
    rng = np.random.default_rng() # this needs to be a random state
    i = rng.choice(n,1)
    P[i] = np.maximum(P[i],0)
    P[i] /= np.sum(P[i])
    j = rng.choice(n,1,p=P[i][0])
    print(i,j)
    return i,j

def findIndex(dim):
    if dim != 0:
        index = sum(len(sublist) for sublist in neighbors_list[0:dim])
    else:
        index = 0
    return index

""" Main Code """

for j in range(1,node_max): # Try for number of nodes
    print("iteration",j,":")
    if (covered_area < min_coverage) or connected_graph == False: # Skip loop if the covered area is good enough & proper connected
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
            if np.sum(adjacency_matrix[i,:]) > 1 and connection_type == "Connected": # always a self loop! so needs to be bigger then 1
                connected_graph = True
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

#x0 = np.random.uniform(low=0, high=20.0, size=num_nodes)
x0 = np.random.normal(10, 2, size=num_nodes)
x_iterating = x0.copy()

x_list_RG = []
P = randomized_gossip(adjacency_matrix) # weighted values for x[k]?

for k in range(max_iter_RG):
    i,j = find_partner(P) # here we need to use find_partner to find the specific index of node k to place the average in!
    x_mean = (1/2)*(x_iterating[i] + x_iterating[j])
    x_iterating[i] = x_mean
    x_iterating[j] = x_mean
    x_list_RG.append(x_iterating.copy()) # x_list can keep track of all values

## PDMM:
n = num_nodes
m = num_edges
test = np.zeros((2*m,1))
neighbors_list = []
for i in range(n):
    neighbors_i = []
    for j in range(n):
        if (adjacency_matrix[i,j] == 1 and i!=j):
            neighbors_i.append(j)
    neighbors_list.append(neighbors_i)

print("PDMM")
x_avg = np.mean(x0)
c_max = np.linspace(start=0.1,stop=0.9,num=9)
a = np.zeros((n,1))
x_PDMM = np.zeros((n,1)) # put x0 as first iteration
x_PDMM = x0.copy()
a = x0.copy()
z = np.zeros((2*m,1))
y_PDMM = np.zeros((2*m,1))
varC_x_list_PDMM = []
total_error_PDMM = []
#print(np.sign(5*np.ones(5)-neighbors_list[5])*z[5])
for c in c_max:
    print("Value for c:",c)
    x_list_PDMM = []
    summed_error = 0
    for time,k in enumerate(range(max_iter_PDMM)):
        #print(time)
        i = np.random.choice(n)
        num_neighbors_i = int(sum(adjacency_matrix[i,:])) -1
        # get the edge index for ij pairs
        if i != 0:
            start_idx_ij = sum(len(sublist) for sublist in neighbors_list[0:i])
        else:
            start_idx_ij = 0
        end_idx_ij = start_idx_ij + num_neighbors_i   
        
        
        neighbor_npList = np.array(neighbors_list[i])
        A_ij = np.sign((neighbor_npList - i*np.ones(num_neighbors_i))).reshape(-1, 1).T
        #A_ij = np.ones(num_neighbors_i).T
        x_PDMM[i] = (a[i] - (A_ij@z[start_idx_ij : end_idx_ij]))/(1+c*num_neighbors_i)
        summed_error += np.linalg.norm(x_PDMM - x_avg * np.ones_like(x_PDMM))**2
        
        x_list_PDMM.append(x_PDMM.copy())
        #print(x_PDMM)
        for index,j in enumerate(neighbors_list[i]):
            #y_PDMM[start_idx + index] = z[start_idx + index] + 2*c*(x_PDMM[i])
            if i<j:
                y_PDMM[start_idx_ij + index] = z[start_idx_ij + index] + 2*c*(1*x_PDMM[i])
            else:
                y_PDMM[start_idx_ij + index] = z[start_idx_ij + index] + 2*c*(-1*x_PDMM[i])
        # sending!
        for index,j in enumerate(neighbors_list[i]):        
            # get the edge index for ji pairs
            if j != 0:
                start_idx_ji = sum(len(sublist) for sublist in neighbors_list[0:j])
            else:
                start_idx_ji = 0
            idx_ji = start_idx_ji + neighbors_list[j].index(i) # connection of ji
            z[idx_ji] = y_PDMM[start_idx_ij + index]
    total_error_PDMM.append(summed_error)
    varC_x_list_PDMM.append(x_list_PDMM)
best_c_index = np.argmin(np.array(total_error_PDMM)) # this value represents the index of c basically
print("Best value for c",best_c_index)

print("PDMM using median")
c_max = np.linspace(start=0.1,stop=0.9,num=9)
x_PDMM_median = np.zeros((n,1)) # put x0 as first iteration
x_PDMM_median = x0.copy()
s = x0.copy()
#z = np.zeros((2*m,1))
mu,sigma = [1,1]
z = np.random.normal(1*mu,sigma,2*m)
y_PDMM_median = np.zeros((2*m,1))
varC_x_list_PDMM_median = []
total_error_PDMM_median = []
for c in c_max:
    print("Value for c:",c)
    x_list_PDMM_median = []
    for time,k in enumerate(range(max_iter_PDMM)):
        #print(time)
        summed_error = 0
        for i in range(n):
            num_neighbors_i = int(sum(adjacency_matrix[i,:])) -1
            # get the edge index for ij pairs
            start_idx_ij = findIndex(i)
            end_idx_ij = start_idx_ij + num_neighbors_i
        
            neighbor_npList = np.array(neighbors_list[i])
            A_ij = np.sign((neighbor_npList - i*np.ones(num_neighbors_i))).reshape(-1, 1).T
            summed_factor_Az = A_ij@z[start_idx_ij : end_idx_ij]
            if ((-1-summed_factor_Az)/(c*num_neighbors_i) > s[i]):
                x_PDMM_median[i] = (-1 - (summed_factor_Az))/(c*num_neighbors_i)
            elif ((1-summed_factor_Az)/(c*num_neighbors_i) < s[i]):
                x_PDMM_median[i] = (1 - (summed_factor_Az))/(c*num_neighbors_i)
            else:
                x_PDMM_median[i] = s[i] # just the initial value

            summed_error += np.linalg.norm(x_PDMM_median - x_avg * np.ones_like(x_PDMM_median))**2
            
            for index,j in enumerate(neighbors_list[i]):
                start_idx_ji = findIndex(j)
                idx_ji = start_idx_ji + neighbors_list[j].index(i) # connection of ji
                if i<j:
                    z[idx_ji] = (1/2)*z[idx_ji] + (1/2)*(z[start_idx_ij + index] + 2*c*(1*x_PDMM_median[i]))
                else:
                    z[idx_ji] = (1/2)*z[idx_ji] + (1/2)*(z[start_idx_ij + index] + 2*c*(-1*x_PDMM_median[i]))

        # sending xi to j?
            for index,j in enumerate(neighbors_list[i]):
                start_idx_ji = findIndex(j)
                idx_ji = start_idx_ji + neighbors_list[j].index(i) # connection of ji    
                if i<j:
                    z[idx_ji] = (1/2)*z[idx_ji] + (1/2)*(z[start_idx_ij + index] + 2*c*(1*x_PDMM_median[i]))
                else:
                    z[idx_ji] = (1/2)*z[idx_ji] + (1/2)*(z[start_idx_ij + index] + 2*c*(-1*x_PDMM_median[i]))
                #z[start_idx_ij+index] = z[idx_ji]
        x_list_PDMM_median.append(x_PDMM_median.copy())
    total_error_PDMM_median.append(summed_error)
    varC_x_list_PDMM_median.append(x_list_PDMM_median)
best_c_index_median = np.argmin(np.array(total_error_PDMM_median)) # this value represents the index of c basically






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

x_array = np.array(x_list_RG)  # shape: (iterations+1, num_nodes)

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

x_array = np.array(varC_x_list_PDMM[best_c_index])  # shape: (iterations+1, num_nodes)

plt.figure(figsize=(10, 6))
num_nodes = x_array.shape[1]

for i in range(num_nodes):
    plt.plot(x_array[:, i])

plt.xlabel('Iteration $k$')
plt.ylabel('Value $x_i[k]$')
plt.title('PDMM Convergence Over Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

x_array = np.array(varC_x_list_PDMM_median[best_c_index_median])  # shape: (iterations+1, num_nodes)

plt.figure(figsize=(10, 6))
num_nodes = x_array.shape[1]

for i in range(num_nodes):
    plt.plot(x_array[:, i])

plt.xlabel('Iteration $k$')
plt.ylabel('Value $x_i[k]$')
plt.title('Median PDMM Convergence Over Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

print("average of x0:",x_avg)
# Compute errors, where PDMM error is plotted using best value for c:

#x_list_PDMM_best = [min(np.linalg.norm(x_k - x_avg * np.ones_like(x_k))**2) for x_k in varC_x_list_PDMM[i] for i in varC_x_list_PDMM]
errors_PDMM = [np.linalg.norm(x_k - x_avg * np.ones_like(x_k))**2 for x_k in varC_x_list_PDMM[best_c_index]]
errors_PDMM_median = [np.linalg.norm(x_k - x_avg * np.ones_like(x_k))**2 for x_k in varC_x_list_PDMM_median[best_c_index_median]]
errors_gossip = [np.linalg.norm(x_k - x_avg * np.ones_like(x_k))**2 for x_k in x_list_RG]

plt.figure(figsize=(10, 6))
plt.semilogy(errors_PDMM, label=f'PDMM (best $c$ = {c_max[best_c_index]:.2f})', linewidth=2)
plt.semilogy(errors_PDMM_median, label=f'Median PDMM (best $c$ = {c_max[best_c_index_median]:.2f})', linewidth=2)
plt.semilogy(errors_gossip, label='Randomized Gossip', linewidth=2, linestyle='--')

plt.xlabel('Iteration $k$')
plt.ylabel(r'$||x(k) - x_{\mathrm{avg}}\cdot\mathbf{1}||^2$')
plt.yscale('log')
plt.title('Consensus Error Over Time: PDMM vs Gossip')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 6))

for idx, x_list in enumerate(varC_x_list_PDMM):
    c_val = c_max[idx]
    errors = [np.linalg.norm(x_k - x_avg * np.ones_like(x_k))**2 for x_k in x_list]
    plt.semilogy(errors, label=f'c = {c_val:.2f}')

plt.xlabel('Iteration $k$')
plt.ylabel(r'$||x(k) - x_{\mathrm{avg}}\cdot\mathbf{1}||^2$')
plt.title('Consensus Error Over Time for Different $c$ Values (PDMM)')
plt.grid(True, which="both", linestyle="--")
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 6))

for idx, x_list in enumerate(varC_x_list_PDMM_median):
    c_val = c_max[idx]
    errors = [np.linalg.norm(x_k - x_avg * np.ones_like(x_k))**2 for x_k in x_list]
    plt.semilogy(errors, label=f'c = {c_val:.2f}')

plt.xlabel('Iteration $k$')
plt.ylabel(r'$||x(k) - x_{\mathrm{avg}}\cdot\mathbf{1}||^2$')
plt.title('Consensus Error Over Time for Different $c$ Values (Median PDMM)')
plt.grid(True, which="both", linestyle="--")
plt.legend()
plt.tight_layout()
plt.show()