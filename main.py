import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import random
import copy
import csv
from scipy.integrate import odeint
from pyvis.network import Network

def main(nodes, q, iters, d, ratio):
    ratio = float(ratio)
    G = nx.Graph()
    nodes = int(nodes)
    q = float(q)
    iters = int(iters)
    d = int(d)

    #SETUP CSV DATA LOGGER
    f = open('data.csv','a', encoding='UTF8', newline='')
    writer = csv.writer(f)
    # header = ['Number on Nodes Backbone', 'Ratio nodes adaptive', 'avg degree, bb', 'avg degree. c', 'avg_err']
    # writer.writerow(header)

    #CONSTANTS
    global avg_c_degree
    avg_c_degree = 0
    B1 = 1
    B2 = 1.5
    B = 5
    Di = d
    K = 0.2
    NODES = nodes
    AR = 0.2
    BR = AR
    CR = 9
    Q = q
    DAMP = 0.2
    #Initialise time
    global t, synce
    t, synce = 0, 0.0
    ratio = ratio
    global num_adaptive
    num_adaptive = 0

    G = nx.barabasi_albert_graph(NODES, 3)
    pos = nx.spring_layout(G)

    Nodes_list_B = list(G.nodes)
    global Central_Nodes, Outer_Nodes
    Central_Nodes = []
    Outer_Nodes = []

    for i in range(NODES):
        degree = G.degree[i]
        if int(degree) > Di:
            Central_Nodes.append(i)
        else:
            Outer_Nodes.append(i)


    #Create Adaptive graph
    def is_adaptive(node, ratio):
        num_of_total_adaptive = int(NODES * 0.7)
        num_of_central_adaptive = int(num_of_total_adaptive*ratio)
        percent_of_central = 0
        if num_of_central_adaptive != 0:
            percent_of_central = len(Central_Nodes)/num_of_central_adaptive
        else:
            percent_of_central = 0
        if node in Central_Nodes:
            if random.random() < percent_of_central:
                return(1)
            else:
                return(0)
            
        else:
            if random.random() > percent_of_central:
                return(1)
            else:
                return(0)

    for i in (Nodes_list_B):
        G.nodes[i]['Adaptive'] = is_adaptive(i, ratio)
        G.nodes[i]['x_state'] = 0
        G.nodes[i]['y_state'] = 0 
        G.nodes[i]['z_state'] = 0

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
        for p in range(NODES):
            if graph.has_edge(node,p):
                num_of_links += 1
        deltai = max(0, num_of_links - Di)
        #print("numi: " + str(num_of_links))
        return(deltai)

    def calculate_deltaj(node, graph):
        num_of_links = 0
        for p in range(NODES):
            if graph.has_edge(p,node):
                num_of_links += 1
        deltai = max(0, num_of_links - Di)
        #print("numj: " + str(num_of_links))
        return(deltai)


    def Heaviside(num):
        if num < 0:
            return 0
        elif num > 0:
            return 1

    def calculate_deltaij(node1, node2, graph):
        deltaij = calculate_deltai(node1, graph) + calculate_deltaj(node2, graph)
        return deltaij

    def calculate_eij(node1, node2, graph):
        eij = graph.nodes[node1]['x_state'] - graph.nodes[node2]['x_state']
        # print("eij: " + str(eij))
        return(eij)

    def calculate_Uij(node1, node2, graph):
        eij = abs(calculate_eij(node1, node2, graph))
        deltij = calculate_deltaij(node1, node2, G_adaptive)
        #print("deltij: " + str(deltij))
        #print("eij: "  + str(eij))
        Uij = B1*eij - B2*deltij
        return Uij

    def calculate_Sigma(node1, node2, graph, t):
        uij = calculate_Uij(node1, node2, graph)
        # h = t/30
        # k1 = 0
        # l1 = h*(uij)
        # k2 = h*(l1/2)
        # l2 = h*(uij - DAMP*(l1/2) + 2*B*(k1/2) + 6*B*((k1/2)**2) + 4*B*(k1/2))
        # k3 = h*(l2/2)
        # l3 = h*(uij - (DAMP*(l2/2) + 2*B*(k2/2) + 6*B*((k2/2)**2) + 4*B*(k2/2)))
        # k4 = l3/2
        # # print("h: " + str(h))
        # # print("l1: " + str(l1))
        # # print("k2: " + str(k2))
        # # print("k3: " + str(k3))
        # sigma = 0
        # iter_thru = np.linspace(0, t, 30)
        # for time in iter_thru:
        #     sigma = sigma + 1/6*(k1 + 2*k2 + 2*k3 +k4)

        def f(sigma, t):
            return (sigma[1], uij - (DAMP*(sigma[1])) - ((2*B*(sigma[0]-1))*sigma[0]*(2*(sigma[0])-1)))

        y0 = [0.1,0.1]
        step = 0

        xs = np.arange(0, t+1, 1)
        sigma, sigma_dev = odeint(f, y0, xs).T
        sigma_ij = sigma[-1]
        return(uij, sigma_ij)

    def Sync_E(graph):
        eij_sum = 0
        #NESTED FOR LOOPS
        #FOR LOOP TO GO THROUGH I
        for i in range(NODES):
            #FOR LOOP TO GO THROUGH J, EXCEPT WHEN J = I
            for j in range(NODES):
                if i != j:
                #ADD EIJ
                    eij_sum += abs(calculate_eij(i, j, graph))
                else:
                    eij_sum = eij_sum
        et = np.sqrt((1/(NODES*(NODES - 1)))*eij_sum)
        return(et)

    def Synchronisation(graph):
        global synce, t
        avg_c_degree = 0
        xplus, yplus, zplus = 0, 0, 0
        x_dot, y_dot, z_dot = 0,0,0
        for m in range(iters):
            #print("iteration = " + str(m))

            tplus = t + 1

            for i in range(NODES):
                iter_thru = np.linspace(t, tplus, 4)
                for l in iter_thru:
                    x_dot = x_der(graph.nodes[i]['x_state'], graph.nodes[i]['y_state'], graph.nodes[i]['z_state'])
                    #print(x_dot)
                    y_dot = y_der(graph.nodes[i]['x_state'], graph.nodes[i]['y_state'], graph.nodes[i]['z_state'], i)
                    #print(y_dot)
                    z_dot = z_der(graph.nodes[i]['x_state'], graph.nodes[i]['y_state'], graph.nodes[i]['z_state'])
                    #print(z_dot)
                    tplus = l
                    xplus = graph.nodes[i]['x_state'] + 0.25*(x_dot)
                    yplus = graph.nodes[i]['y_state'] + 0.25*(y_dot)
                    zplus = graph.nodes[i]['z_state'] + 0.25*(z_dot)   
                    graph.nodes[i]['x_state'] = xplus
                    graph.nodes[i]['y_state'] = yplus
                    graph.nodes[i]['z_state'] = zplus

                    t = tplus
                for j in range(NODES):
                    if(G_adaptive.has_node(i) and G_adaptive.has_node(j)):
                        uij, sigma = calculate_Sigma(i,j, graph, t)
                        #print("uij: " + str(uij))
                        print("sigma: " + str(sigma))
                        if (sigma) > 0.5:
                            if (i != j and not G_adaptive.has_edge(i,j)) :
                                G_adaptive.add_edge(i, j)
                                print("edge created!")
                        elif (sigma) < 0.5 and i != j and G_adaptive.has_edge(i,j):
                            G_adaptive.remove_edge(i,j)
                            print("edge broken!")

            synce = synce + Sync_E(G)
            degree_of_graph = 2*G_adaptive.number_of_edges()/G_adaptive.number_of_nodes()
            avg_c_degree = avg_c_degree + degree_of_graph
        avg_c_degree = avg_c_degree/iters
        synce = synce/iters
                

    def y_der(x, y, z, i):
        sum_b = 0
        sum_c = 0

        for j in range(NODES):
            diff = (G.nodes[j]['y_state'] - y)
            sum_b = sum_b + (G.has_edge(i, j)*diff)
            sum_c = sum_c + (G_adaptive.has_edge(i,j)*diff)

        y_dot = x + AR*y + K*(sum_b) + Q*(sum_c)
        return(y_dot)    

    def x_der(x,y,z):
        x_dot = -y-z
        return(x_dot)

    def z_der(x,y,z):
        z_dot = BR + z*(z - CR)
        return(z_dot)

    def draw_graph(nx_graph):
        fig, axes = plt.subplots(1,1,dpi=72)
        nx.draw(nx_graph, pos=pos, ax=axes, with_labels=True)



    def run_sync():
        Synchronisation(G)


    run_sync()
    draw_graph(G)
    #avg_c_degree = (2*G_adaptive.number_of_edges())/G_adaptive.number_of_nodes()
    print("Average degree (Backbone): " + str(avg_b_degree))
    print("Average degree (Control): " + str(avg_c_degree))
    print("Average Error: " + str(synce))
    data = [NODES, ratio, Di, avg_b_degree, avg_c_degree, synce]
    
    writer.writerow(data)
    f.close()
    name_b = str(q)
    name_a = "plots_ratio_5/" + name_b + "A.png"
    name_b = "plots_ratio_5/" + name_b + ".png"
    plt.savefig(name_b)
    draw_graph(G_adaptive)
    plt.savefig(name_a)
    return(synce, avg_c_degree)



