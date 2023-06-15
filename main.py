import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import random
import copy
import scipy
from pyvis.network import Network


G = nx.Graph()

#CONSTANTS

B1 = 1
B2 = 2.5
B = 5
Di = 3
K = 0.2
NODES = 50
AR = 0.2
BR = AR
CR = 9
Q = 0.1
DAMP = 0.2

#Initialise time
t = 5

G = nx.barabasi_albert_graph(NODES, 2)
pos = nx.spring_layout(G)

Nodes_list_B = list(G.nodes)


#Create Adaptive graph
def is_adaptive(number):
    if random.random() < number:
        return(1)
    else:
        return(0)

for i in (Nodes_list_B):
    G.nodes[i]['Adaptive'] = is_adaptive(0.5)
    G.nodes[i]['x_state'] =random.random() 
    G.nodes[i]['y_state'] =random.random() 
    G.nodes[i]['z_state'] =random.random() 

G_adaptive = copy.deepcopy(G)
G_adaptive.remove_edges_from(G_adaptive.edges())

AdaptiveList = []
for i in (Nodes_list_B):
    AdaptiveList.append(G.nodes[i]['Adaptive'])

print(AdaptiveList)
avg_b_degree = (2*(G.number_of_edges()))/G.number_of_nodes()

for i in range(len(AdaptiveList)):
    if AdaptiveList[i] == 0:
        G_adaptive.remove_node(i)

#ADJACENCY MATRIX FOR BACKBONE
Ab = nx.adjacency_matrix(G)

#ADJACENCY MATRIX FOR CONTROL
Ac = nx.adjacency_matrix(G_adaptive)

#CALCULATE Uij

def calculate_deltai(node, graph):
    num_of_links = 0
    for i in range(NODES):
        if graph.has_edge(node,i):
            num_of_links += 1
    deltai = max(0, num_of_links - Di)
    return(deltai)

def Heaviside(num):
    if num < 0:
        return 0
    elif num > 0:
        return 1

def calculate_deltaij(node1, node2, graph):
    deltaij = calculate_deltai(node1, graph) + calculate_deltai(node2, graph)
    return deltaij

def calculate_eij(node1, node2, graph):
    eij = graph.nodes[node1]['x_state'] - graph.nodes[node2]['x_state']
    return(eij)

def calculate_Uij(node1, node2, graph):
    Uij = B1*abs(calculate_eij(node1, node2, graph)) - B2*calculate_deltaij(node1, node2, graph)
    return Uij

def calculate_Sigma(node1, node2, graph, t):
    uij = calculate_Uij(node1, node2, graph)
    h = t/4
    k1 = 0
    l1 = h*(uij)
    k2 = h*(l1/2)
    l2 = h*(uij - DAMP*(l1/2) + 2*B*(k1/2) + 6*B*((k1/2)**2) + 4*B*(k1/2))
    k3 = h*(l2/2)
    l3 = h*(uij - (DAMP*(l2/2) + 2*B*(k2/2) + 6*B*((k2/2)**2) + 4*B*(k2/2)))
    k4 = l3/2
    print("h: " + str(h))
    print("l1: " + str(l1))
    print("k2: " + str(k2))
    print("k3: " + str(k3))
    sigma = 0
    iter_thru = np.linspace(0, t, 4)
    for time in iter_thru:
        sigma = sigma + 1/6*(k1 + 2*k2 + 2*k3 +k4)

    return(uij, sigma)



def Sync_E(graph):
    eij_sum = 0
    #NESTED FOR LOOPS
    #FOR LOOP TO GO THROUGH I
    for i in enumerate(list(graph.nodes)):
        #FOR LOOP TO GO THROUGH J, EXCEPT WHEN J = I
        for j in enumerate(list(graph.nodes)):
            if i != j:
            #ADD EIJ
                eij_sum += abs(calculate_eij(i, j, graph))
            else:
                eij_sum = eij_sum
    et = np.sqrt((1/(NODES(NODES)))*eij_sum)

def Synchronisation(graph):
    for i in range(NODES):
        for j in range(NODES):


            uij, sigma = calculate_Sigma(i,j, graph, t)
            print("uij: " + str(uij))
            print("sigma: " + str(sigma))
            if abs(sigma) > 0.5:
                if (i in G_adaptive and j in G_adaptive):
                    G_adaptive.add_edge(i, j)
            

def y_der(x, y, z, i):
    sum_b = 0
    sum_c = 0

    for j in range(NODES):
        diff = (G.nodes[j]['y_state'] - y)
        sum_b = sum_b + (G.has_edge(i, j)*diff)
        sum_c = sum_c + (G_adaptive.has_edge(i,j)*diff)

    y_dot = x + AR*y + K*(sum_b) + Q(sum_c)
    return(y_dot)    

def x_der(x,y,z):
    x_dot = -y-z
    return(x_dot)

def z_der(x,y,z):
    z_dot = BR + z(z - CR)
    return(z_dot)

def draw_graph(nx_graph):
    fig, axes = plt.subplots(1,1,dpi=72)
    nx.draw(nx_graph, pos=pos, ax=axes, with_labels=True)



def run_sync():
    Synchronisation(G)




run_sync()
draw_graph(G)
avg_c_degree = (2*G_adaptive.number_of_edges())/G_adaptive.number_of_nodes()
print("Average degree (Backbone): " + str(avg_b_degree))
print("Average degree (Control): " + str(avg_c_degree))
draw_graph(G_adaptive)
plt.show()




