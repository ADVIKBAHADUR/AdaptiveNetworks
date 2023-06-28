import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import random
import copy
import csv
from scipy.integrate import odeint
from pyvis.network import Network

g = nx.Graph()
G_adaptive = nx.Graph()
global Central_Nodes, Outer_Nodes, adaptiveOuter, adaptiveCentral, pos

Central_Nodes, Outer_Nodes, adaptiveCentral, adaptiveOuter = [], [], [], []
pos = nx.spring_layout(g)
class AdaptiveNetwork:
    def __init__(self, NumberNodes, q, AdaptibilityRatio, CentralityRatio, times, d):
        self.NumberNodes = NumberNodes
        self.d = d
        self.q = q
        self.AdaptibilityRatio = AdaptibilityRatio
        self.CentralityRatio = CentralityRatio
        self.times = times

    def create_B_graph(self):
        g = nx.barabasi_albert_graph(self.NumberNodes, self.d)
        self.draw_graph(g)

    def create_A_graph(self):
        degree_list = list(g.degree())
        print(degree_list)
        for i in range(self.NumberNodes):
            degree = degree_list[i][1]
            if int(degree) > self.d:
                Central_Nodes.append(i)
            else:
                Outer_Nodes.append(i)
            for i in range(self.NumberNodes):
                if self.is_adaptive(i):
                    G_adaptive.add_node(i)
                
    def is_adaptive(self, node):
        num_adaptive_nodes = int(self.NumberNodes/self.AdaptibilityRatio)
        num_central_adaptive = int(num_adaptive_nodes * self.CentralityRatio)
        adaptiveCentral = random.sample(Central_Nodes, num_central_adaptive)
        adaptiveOuter = random.sample(Outer_Nodes, self.NumberNodes - num_central_adaptive)

        if node in adaptiveCentral or node in adaptiveOuter:
            return True
        else:
            return False
    def draw_graph(self, nx_graph):
        fig, axes = plt.subplots(1,1,dpi=72)
        pos = nx.spring_layout(g)
        nx.draw(nx_graph, pos=pos, ax=axes, with_labels=True)


A = AdaptiveNetwork(100, 0.2, 0.9, 0.9, 10, 4)
A.create_B_graph()
plt.show()
    


