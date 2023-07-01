import random
from collections import deque
import copy
import time 
import sys
import math
import networkx.algorithms.community as nx_comm
from scipy.stats import expon
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls


class ICG(GraphTolls) :
    def __init__(self,G,Nb,Beta):    
        self.G = G
        self.Nb = Nb
        self.Beta = Beta
        self.m = G.number_of_edges()
        self.n = G.number_of_nodes()
        self.Mod_val = 0

    def expon (self , value , teta):
        x = (1/teta)*(value)
        p =float(math.exp(x))
        return p      
    
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
    
    def crousel(self,clusters,preserve_node,drop_node, alpha):
        iterations = int(alpha*self.n)
        cluster = clusters.copy()
        for i in range(iterations):
            node = preserve_node.pop(0)
            drop_node.append(node)
            
            for community in clusters:
                communityy = community.copy()
                for v in community.copy():
                    if node == v:
                        community.discard(v)
                        if community == {}:
                            clusters.remove({})
        
            selected_node = random.choice(drop_node)
            drop_node.remove(selected_node)
            preserve_node.append(selected_node)
        
            MAX_Q = 0
            pos = -1
            for index,community in enumerate(clusters):    
                Kbv = super().select_edge_betw(self.G,node,community)
                db = sum([j for k,j in self.G.degree(community)])
                delta_Q = 1/self.m * Kbv -self.G.degree(node)/(2*self.m**2)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = index
                    
            if MAX_Q > 0:
                if selected_node not in clusters[pos]:    
                    clusters[pos].add(selected_node)
            else:
                clusters.append({selected_node})
            
            
        return clusters,drop_node 
    
    def localsearch(self, clusters):
        
        node_list = [i for i in range(self.n)]
        count = 0
        while count < self.n:
            
            random.shuffle(node_list)
            for node in node_list:
                Q = []
                degree = self.G.degree(node)
                for community in clusters:
                    if node in community:
                        self_community = community
                        Kav = super().select_edge_betw(self.G,node,self_community)
                        da = sum([j for k,j in self.G.degree(self_community)])
                        break            
                for index,community in enumerate(clusters):
                    if node in community:
                        delta_Q = 0
                        before = index
                        Q.append(delta_Q)
                    else:
                        Kbv = super().select_edge_betw(self.G,node,community)
                        db = sum([j for k,j in self.G.degree(community)])
                        delta_Q = (1/self.m) * (Kbv-Kav) + (degree/(2*self.m**2))*(da - db - degree)
                        Q.append(delta_Q)

                if max(Q) > 0:
                    pos = Q.index(max(Q))
                    clusters[pos].add(node)
                    clusters[before].remove(node)
                    if clusters[before] == []:
                        del clusters[before]

                else:
                    count += 1

        membership = [0 for i in range(self.n)]
        for i in range(len(clusters)):
            for index in clusters[i]:
                membership[index] = i
        
        return membership, clusters
        
    def lebel_node (self,community):
        label = sorted([i for i in self.G.nodes()])
        for index,no in enumerate (label):
            for i in range(len(community)):
                if no in community[i]:
                    label[index] = i
        
        return label    
                               
    def Run_ICG (self):
        start = time.time()
        soltion = self.GCH()
        solutionmm, soltion = self.localsearch(soltion)
        best_solution = copy.deepcopy(soltion)
        T_init = 0.025*nx_comm.modularity(self.G, soltion)
        T = T_init
        nb_iter = 0
        while nb_iter < self.Nb:
            Q1 = nx_comm.modularity(self.G, soltion)
            incumbent_solution = copy.deepcopy(soltion)
            soltion,drop_nodes,preserve_node = self.Destruction(soltion)
            soltion, drop_nodes = self.crousel(soltion,preserve_node, drop_nodes, 0.5) 
            soltion = self.Reconstruction(soltion,drop_nodes)
            solutionmm, soltion = self.localsearch(soltion) 
            Q2 = nx_comm.modularity(self.G, soltion)
            if Q2 > nx_comm.modularity(self.G, best_solution):
                best_solution = copy.deepcopy(soltion)
                
            P = random.random()
            if Q2 < Q1 and P > math.exp((Q2 - Q1)//T):
                soltion = copy.deepcopy(incumbent_solution)
                T = T_init

            elif Q2  >=  Q1:
                T = T_init
              
            else:
                T = T*0.9
        
            nb_iter = nb_iter + 1
        
        self.Mod_val = nx_comm.modularity(self.G, best_solution)
        end = time.time()
        t = end-start
        
        return self.Mod_val, best_solution,t                               

def de_main():
    path =  '/home/yacine/Desktop/real_network/football.gml'
    Number_iter = 200
    Beta = 0.4
    data = GraphTolls()
    graph = data.Read_Graph(path)
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < 2 :
        communities = ICG(graph, Number_iter, Beta)
        mod,community,tim = communities.Run_ICG() 
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









    
                
