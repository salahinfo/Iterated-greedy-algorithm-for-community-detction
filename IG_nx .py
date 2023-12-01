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
    def __init__(self,graph,Nb,Beta, path):    
        self.Nb = Nb
        self.Beta = Beta
        self.Mod_val = 0
        super( IG, self).__init__(path)
     
    
    def GCH( self):
         
        vertex_list = list(self.graph.nodes())
        node = random.choice(vertex_list)
        #print("nnnn",node)
        com_id = 0
        self.membership[node] = com_id
        self.DegCom[com_id] = self.Degree[node]
        #print(self.Degree)
        vertex_list.remove(node)
        for node in  vertex_list:    
            comm_ngh = super().neigh_comm( node)
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
                super().insert_node( node, pos, comm_ngh.get( pos , 0))
            else:
                com_id = com_id + 1
                super().insert_node(node , com_id, comm_ngh.get(com_id, 0))
                                     
        
        return  self.membership 
    
    def Destruction( self):       
        drop_node = []
        merg_node = []
        cut_len = int(len(self.membership)* float(self.Beta)) 
        index_community = random.sample( list(self.membership.keys()), cut_len )
        #print("list", cut_len, index_community)
        #print("degcomunity", self.DegCom)
        #print("mm",self.membership)
        for al in index_community:
            com_id = self.membership[al]
            wgh = super().neigh_comm(al)    
            super().delet_node(al, com_id, wgh.get(com_id, 0.))
            if self.internal[com_id] == 0.:
                del self.internal[com_id]            
                      
        #merg_node = [ nod for nod in self.Node_list if nod not in index_community] 
        return  self.membership, index_community

    def reconcstruction(self, drop_node):
        random.shuffle(drop_node)
        for node in  drop_node:    
            comm_ngh = super().neigh_comm( node)
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
                super().insert_node( node, pos, comm_ngh.get( pos, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(self.membership.values()))
                super().insert_node( node, com_id, comm_ngh.get( com_id, 0.))

        return self.membership


    def Run_IG (self):
        start = time.time()
        soltion = self.GCH()
        print(time.time()- start)
        best_solution = copy.deepcopy(soltion)
        best_Q = self.modularity()
        print(best_Q)
        nb_iter = 0
        while nb_iter < self.Nb:

            soltion,drop_nodes = self.Destruction()
            soltion = self.reconcstruction(drop_nodes) 
            super().init( soltion, weight='weight')
            Q2 = self.modularity()
            print(" the value of modularity and time ",Q2, time.time()- start)
            if Q2 > best_Q:
                best_solution = copy.deepcopy(soltion)
                best_Q = Q2
           
        
            nb_iter = nb_iter + 1
        
       
        end = time.time()
        t = end-start
        
        return best_Q, best_solution,t                               

def de_main():
    path =  '/home/yacine/Desktop/real_network/netscience.gml'
    Number_iter = 200
    Beta = 0.4
    data = GraphTolls(path)
    graph = data.Read_Graph()
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < 5:
        communities = IG(graph, Number_iter, Beta,path)
        mod,community,tim = communities.Run_IG() 
        Q_list.append(mod)
        Time_list.append(tim)
        #label = communities.lebel_node(community)  
        #True_partition = data.Read_GroundTruth('/home/yacine/Desktop/real_network/polbooks.txt')
        #NMI = normalized_mutual_info_score(True_partition, list(community.values()))
        #NMI_list.append(NMI)      
        nb_run = nb_run +1
    
    
    Q_avg = communities.avg(Q_list)
    Q_max = communities.max(Q_list)
    Q_std = communities.stdev(Q_list)
    NMI_max = communities.max(NMI_list)
    time_run = communities.avg(Time_list)
    print("Q_avg",Q_avg,"Q_max",Q_max,"NMI",NMI_max,time_run,Q_std)
    
if __name__ == '__main__':
    
    de_main()









    
                
