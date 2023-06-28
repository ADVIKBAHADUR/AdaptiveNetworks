import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import random
import copy
import csv
from scipy.integrate import odeint
from pyvis.network import Network

ratio = 0.1
G = nx.barabasi_albert_graph(10, 3)
pos = nx.spring_layout(G)

name_b = str(ratio)
name_a = "plots_ratio/" + name_b + "A.png"
name_b = "plots_ratio/" + name_b + ".png"

def draw_graph(nx_graph):
    fig, axes = plt.subplots(1,1,dpi=72)
    nx.draw(nx_graph, pos=pos, ax=axes, with_labels=True)
draw_graph(G)
plt.savefig(name_b)