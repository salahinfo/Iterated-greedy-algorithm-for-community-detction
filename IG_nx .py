import random
from collections import deque
import copy
import time 
import sys
import math
import networkx.algorithms.community as nx_comm
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls


class IG(GraphTolls) :
    def __init__(self,G,Nb,Beta):    
        self.G = G
        self.Nb = Nb
        self.Beta = Beta
        self.m = G.number_of_edges()
        self.n = G.number_of_nodes()
        self.Mod_val = 0
     
    
    def GCH(self):
        community = []
        vertex_list = [i for i in self.G.nodes()]
        node = random.choice(vertex_list)
        vertex_list.remove(node)
        community.append({node})
        while vertex_list != []:    
            node = random.choice(vertex_list)
            vertex_list.remove(node)
            MAX_Q = 0
            pos = -1
            for index,clusters in enumerate(community):
                if super().is_edge_betw(self.G,node,clusters):        
                    Kbv = super().select_edge_betw(self.G,node,clusters)
                    db = sum([j for k,j in self.G.degree(clusters)])
                    delta_Q = 1/self.m * Kbv -self.G.degree(node)/(2*self.m**2)*db
                    if delta_Q > MAX_Q:
                        MAX_Q = delta_Q
                        pos = index
                else :
                    delta_Q = 0

            if MAX_Q > 0:
                community[pos].add(node)
            else:
                community.append({node})
            
        return  community 
    
    def Destruction(self,community):
        vertex_list = [i for i in self.G.nodes()]
        drop_node = []
        cut_len = int(float(self.n)* float(self.Beta))
        random.shuffle(vertex_list)
        drop = vertex_list[self.n-cut_len:]
        pres_node = vertex_list[ :self.n-cut_len]
        for i in range(cut_len) :
            v = drop.pop()
            drop_node.append(v)
            for cluster in community:
                if v in cluster:
                    cluster.remove(v)
                    if len(cluster) == 0: 
                        del cluster 

        return  community , drop_node, pres_node
    
    def Reconstruction(self,community,drop_node):
    
        random.shuffle(drop_node)
        for node in drop_node:
            MAX_Q = 0
            pos = -1
            for index,clusters in enumerate(community):  
                if super().is_edge_betw(self.G,node,clusters): 
                    Kbv = super().select_edge_betw(self.G,node,clusters)
                    db = sum([j for k,j in self.G.degree(clusters)])
                    delta_Q = 1/self.m * Kbv - self.G.degree(node)/(2*self.m**2)*db
                    if delta_Q > MAX_Q:
                        MAX_Q = delta_Q
                        pos = index
    
            if MAX_Q > 0:
                community[pos].add(node)
            else:
                community.append({node})
                
                
        return community
    
    
    def lebel_node (self,community):
        label = sorted([i for i in self.G.nodes()])
        for index,no in enumerate (label):
            for i in range(len(community)):
                if no in community[i]:
                    label[index] = i
        
        return label    
                               
    def Run_IG (self):
        start = time.time()
        soltion = self.GCH()
        best_solution = copy.deepcopy(soltion)
        best_Q = nx_comm.modularity(self.G, soltion)
        nb_iter = 0
        while nb_iter < self.Nb:

            soltion,drop_nodes,preserve_node = self.Destruction(soltion)
            soltion = self.Reconstruction(soltion,drop_nodes) 
            Q2 = nx_comm.modularity(self.G, soltion)
            if Q2 > best_Q:
                best_solution = copy.deepcopy(soltion)
                best_Q = Q2
           
        
            nb_iter = nb_iter + 1
        
       
        end = time.time()
        t = end-start
        
        return best_Q, best_solution,t                               

def de_main():
    path =  '/home/yacine/Desktop/real_network/football.gml'
    Number_iter = 200
    Beta = 0.3
    data = GraphTolls()
    graph = data.Read_Graph(path)
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < 10 :
        communities = IG(graph, Number_iter, Beta)
        mod,community,tim = communities.Run_IG() 
        Q_list.append(mod)
        Time_list.append(tim)
        label = communities.lebel_node(community)  
        True_partition = data.Read_GroundTruth('/home/yacine/Desktop/real_network/football.txt')
        NMI = normalized_mutual_info_score(True_partition,label)
        NMI_list.append(NMI)      
        nb_run = nb_run +1
    
    
    Q_avg = communities.avg(Q_list)
    Q_max = communities.max(Q_list)
    Q_std = communities.stdev(Q_list)
    NMI_max = communities.max(NMI_list)
    time_run = communities.avg(Time_list)
    print("Q_avg",Q_avg,"Q_max",Q_max,"NMI",NMI_max,time_run)
    
if __name__ == '__main__':
    
    de_main()









    
                
