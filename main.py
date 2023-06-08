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
Di = 4
NODES = 8

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
    G.nodes[i]['x_state'] =random.random() 
    G.nodes[i]['Adaptive'] = is_adaptive(0.5)

G_adaptive = copy.deepcopy(G)

AdaptiveList = []
for i in (Nodes_list_B):
    AdaptiveList.append(G.nodes[i]['Adaptive'])

print(AdaptiveList)
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

def calculate_deltaij(node1, node2, graph):
    deltaij = calculate_deltai(node1, graph) + calculate_deltai(node2, graph)
    return deltaij

def calculate_eij(node1, node2, graph):
    eij = graph.nodes[node1]['x_state'] - graph.nodes[node2]['x_state']
    return(eij)

def calculate_Uij(node1, node2, graph):
    Uij = B1*abs(calculate_eij(node1, node2, graph)) - B2*calculate_deltaij(node1, node2, graph)
    return Uij

def Sync_E(graph):
    #NESTED FOR LOOPS
    #FOR LOOP TO GO THROUGH I
        #FOR LOOP TO GO THROUGH J, EXCEPT WHEN J = I
            #ADD EIJ
    pass

def draw_graph(nx_graph):
    fig, axes = plt.subplots(1,1,dpi=72)
    nx.draw(nx_graph, pos=pos, ax=axes, with_labels=True)
  

draw_graph(G)
draw_graph(G_adaptive)
plt.show()




