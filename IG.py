import random
from collections import deque
import copy
import time 
import sys
import math
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls
from cdlib import algorithms, evaluation
import networkx as nx
import igraph as ig
from cdlib import NodeClustering, datasets
import leidenalg as la
 


class IG(GraphTolls) :
    def __init__(self,graph,Nb,Beta, path):    
        self.Nb = Nb
        self.Beta = Beta
        self.Mod_val = 0
        super( IG, self).__init__(path)
     
    
    def GCH( self):
        
        membership = { i : None for i in self.graph.nodes() }
        print(membership)
        vertex_list = list(self.graph.nodes())
        node = random.choice(vertex_list)
        #print("nnnn",node)
        com_id = 0
        membership[node] = com_id
        self.DegCom[com_id] = self.Degree[node]
        #print(self.Degree)
        vertex_list.remove(node)
        for node in  vertex_list:    
            comm_ngh = super().neigh_comm( node, membership)
            MAX_Q = 0
            pos = -1
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q =  Kbv - self.Degree[node]/(2.*self.m)*db
                #print(delta_Q)
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                membership = super().insert_node( node, membership, pos, comm_ngh.get( pos , 0))
            else:
                com_id = com_id + 1
                membership = super().insert_node(node , membership, com_id, comm_ngh.get(com_id, 0))
                                     
        
        return  membership
    
   # def Destruction( self, membership, graph):
        drop_node= []
        #s = self.renumber(membership)
        membership = super().init(  membership , weight='weight') 
        cut_len = int(len(membership)* float(self.Beta)) 
        drop_node = random.sample( list(membership.keys()), cut_len )
        for al in drop_node:
            com_id = membership[al]
            wgh = super().neigh_comm(al, membership)    
            membership = super().delet_node( al, membership,  com_id, wgh.get( com_id, 0.))
            if al not in set(membership.values()):
                membership = super().insert_node( al, membership, al, wgh.get( al, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(membership.values()))
                membership = super().insert_node( al, membership, com_id, wgh.get( com_id, 0.))     
        
        #membership = super().renumber( membership)
        return  membership, drop_node 

     
    
    def Destruction( self, membership):       
        drop_node = []
        merg_node = []
        cut_len = int(len(membership)* float(self.Beta)) 
        index_community = random.sample( list(membership.keys()), cut_len )
        #print("list", cut_len, index_community)
        #print("degcomunity", self.DegCom)
        #print("mm",self.membership)
        for al in index_community:
            com_id = membership[al]
            wgh = super().neigh_comm(al, membership)    
            membership = super().delet_node(al, membership, com_id, wgh.get(com_id, 0.))
            if self.internal[com_id] == 0.:
                del self.internal[com_id]            
                      
        #merg_node = [ nod for nod in self.Node_list if nod not in index_community] 
        return  membership, index_community
    

    def Leiden_Recon(self , memb):
        #for key in memb.keys():
         #   print(f"Key: {key}, Type: {type(key)}")
        
        #print(len(memb))
        #memb = {node-1: node-1 % 2 for node in self.graph.nodes()}
        memb = dict(sorted(memb.items()))
        membershi = list(memb.values())
        print("mmm",membershi)
        cl = ig.Clustering(membershi)
        print("mmm",cl.membership)
        graph = ig.Graph.from_networkx(self.graph)
        communities = la.find_partition( graph, la.ModularityVertexPartition, initial_membership = cl.membership)
        #communities = algorithms.leiden( self.graph , initial_membership = initial_membership)
        node_to_community = {}
        for community_id, community in enumerate(communities):
            for node in community:
                node_to_community[node] = community_id

        return node_to_community

    def reconcstruction(self, drop_node, solution):
        random.shuffle(drop_node)
        for node in  drop_node:    
            comm_ngh = super().neigh_comm( node, solution)
            MAX_Q = 0
            pos = -1
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q =  Kbv - self.Degree[node]/(2.*self.m)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                solution = super().insert_node( node, solution, pos, comm_ngh.get( pos, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(self.membership.values()))
                solution = super().insert_node( node, solution, com_id, comm_ngh.get( com_id, 0.))

        return solution


    def Run_IG (self):
        start = time.time()
        soltion = self.GCH()
        print(time.time()- start)
        best_solution = copy.deepcopy(soltion)
        best_Q = self.modularity(soltion)
        print(best_Q)
        nb_iter = 0
        while nb_iter < self.Nb:

            soltion,drop_nodes = self.Destruction(soltion)
            #print(soltion)
            soltion = self.reconcstruction(drop_nodes, soltion) 
            #soltion = self.Leiden_Recon(soltion)
            print(soltion)
            soltion = super().init( soltion, weight='weight')
            Q2 = self.modularity(soltion)
            print(" the value of modularity and time ",Q2, time.time()- start)
            if Q2 > best_Q:
                best_solution = copy.deepcopy(soltion)
                best_Q = Q2
           
        
            nb_iter = nb_iter + 1
        
       
        end = time.time()
        t = end-start
        
        return best_Q, best_solution,t                               

def de_main():
    #path =  '/home/yacine/Desktop/LFR/network_mu_0.8.dat'
    path = sys.argv[1]
    #'/home/yacine/Desktop/real_network/louvain.txt'
    Number_iter = int(sys.argv[2])
    Beta = sys.argv[3]
    data = GraphTolls(path)
    graph = data.Read_Graph()
    #graph = datasets.fetch_network_data(net_name="karate_club", net_type="igraph")
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < int(sys.argv[5]):
        communities = IG(graph, Number_iter, Beta,path)
        mod,community,tim = communities.Run_IG() 
        Q_list.append(mod)
        Time_list.append(tim)
        #label = communities.lebel_node(community)  
        if sys.argv[4] != 'None':
            community = dict(sorted(community.items()))
            True_partition = data.Read_GroundTruth(sys.argv[4])
            NMI = normalized_mutual_info_score(True_partition, list(community.values()))
            NMI_list.append(NMI)      
        nb_run = nb_run +1
    
    
    #data.writefile(Q_list)
    Q_avg = communities.avg(Q_list)
    Q_max = communities.max(Q_list)
    Q_std = communities.stdev(Q_list)
    NMI_max = communities.max(NMI_list)
    time_run = communities.avg(Time_list)
    print("Q_avg",Q_avg,"Q_max",Q_max,"NMI",NMI_max,time_run,Q_std)
    
if __name__ == '__main__':
    
    de_main()









    
                
