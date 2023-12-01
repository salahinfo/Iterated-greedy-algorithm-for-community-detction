import random
from collections import deque
import copy
import time 
import sys
import math
from scipy.stats import expon
from sklearn.metrics.cluster import normalized_mutual_info_score
from GraphTools import GraphTolls


class ICG(GraphTolls) :
    def __init__( self, G, Nb, Beta,path):    
        self.G = G
        self.Nb = Nb
        self.Beta = Beta
        self.Mod_val = 0
        super( ICG, self).__init__(path)

    def expon (self , value , teta):
        x = (1/teta)*(value)
        p =float(math.exp(x))
        return p      
    
    def GCH( self): 
        vertex_list = list(self.G.nodes())
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
                super().insert_node( node , com_id, comm_ngh.get(com_id, 0))
                                     
        
        return  self.membership 
    
    def Destruction( self):       
        node_list = list(self.G.nodes())
        random.shuffle(node_list)
        cut_len = int(len(self.membership)* float( self.Beta)) 
        preserve_node = node_list[  : cut_len]
        drop_node = node_list[cut_len : ]
        for al in drop_node:
            com_id = self.membership[al]
            wgh = super().neigh_comm(al)    
            super().delet_node(al, com_id, wgh.get(com_id, 0.))
            if self.internal[com_id] == 0.:
                del self.internal[com_id]            
                      
        #merg_node = [ nod for nod in self.Node_list if nod not in index_community] 
        return  self.membership, preserve_node, drop_node
    
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

    
    def crousel(self, preserve_node, drop_node, alpha):
        iterations = int(alpha*self.n)
        for i in range(iterations):
            node = preserve_node.pop(0)
            drop_node.append(node)
            for vertex,com in self.membership.items():
                if node == vertex:
                    wgh = super().neigh_comm(node)  
                    super().delet_node(node, self.membership[node],wgh.get( self.membership[node], 0.))
                    
            selected_node = random.choice(drop_node)
            drop_node.remove(selected_node)
            preserve_node.append(selected_node)
            MAX_Q = 0
            pos = -1    
            comm_ngh = super().neigh_comm( selected_node)
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q =  Kbv - self.Degree[selected_node]/(2.*self.m)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                if self.membership[selected_node] != pos :
                    super().insert_node( selected_node, pos, comm_ngh.get( pos, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(self.membership.values()))
                super().insert_node( selected_node, com_id, comm_ngh.get( com_id, 0.))
        
        return self.membership, drop_node 
    
    def localsearch(self, clusters):
        super().init(clusters)
        node_list = [i for i in self.G.nodes()]
        count = 0
        modified = True
        while modified:
            modified = False
            for vsele in node_list :
                degree = self.Degree[vsele]
                com_befor = self.membership[vsele]
                ngh_com = super().neigh_comm(vsele)
                dvc = super().ngh_node(vsele, com_befor)
                devc = self.DegCom[com_befor]
                maxq = 0
                m_com = com_befor             
                for com, dvcp in ngh_com.items() : 
                    devcp = self.DegCom[com]
                    deq = (1/self.m) * ( dvcp-dvc )- ((degree)/(2.*self.m**2.)) *( devcp-devc+degree )
                
                    if deq > maxq :
                        maxq = deq
                        m_com = com
                        modified = True
                        
                super().delet_node( vsele, com_befor, ngh_com.get( com_befor, 0.))            
                super().insert_node( vsele, m_com, ngh_com.get( m_com, 0.))   
           
        return clusters
    
    def reconcstruction( self, drop_node):
        random.shuffle( drop_node )
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

                         
    def Run_ICG (self):
        start = time.time()
        soltion = self.GCH()
        soltion = self.localsearch(soltion)
        best_solution = copy.deepcopy(soltion)
        q = super().modularity()
        print(q, time.time()-start)
        T_init = 0.025*q
        T = T_init
        nb_iter = 0
        while nb_iter < self.Nb:
            Q1 = super().modularity()
            incumbent_solution = copy.deepcopy( soltion)
            soltion,drop_nodes,preserve_node = self.Destruction()
            soltion, drop_nodes = self.crousel( preserve_node, drop_nodes, 0.7) 
            soltion = self.reconcstruction( drop_nodes)
            soltion = self.localsearch( soltion) 
            Q2 = super().modularity()
            print("q and time ",Q2, time.time() - start)
            if Q2 > super().modularity():
                best_solution = copy.deepcopy( soltion)
                
            P = random.random()
            if Q2 < Q1 and P > math.exp((Q2 - Q1)//T):
                soltion = copy.deepcopy(incumbent_solution)
                T = T_init

            elif Q2  >=  Q1:
                T = T_init
              
            else:
                T = T*0.9
        
            nb_iter = nb_iter + 1
        
        self.Mod_val = super().modularity()
        end = time.time()
        t = end-start
        
        return self.Mod_val, best_solution,t                               

def de_main():
    path =  '/home/yacine/Desktop/real_network/polbooks.gml'
    Number_iter = 300
    Beta = 0.3
    data = GraphTolls(path)
    graph = data.Read_Graph()
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < 1 :
        communities = ICG( graph, Number_iter, Beta,path)
        mod, community, tim = communities.Run_ICG() 
        Q_list.append(mod)
        Time_list.append(tim)
        #label = communities.lebel_node(community)  
        #True_partition = data.Read_GroundTruth('/home/yacine/Desktop/real_network/polbooks.txt')
        #NMI = normalized_mutual_info_score(True_partition,list(community.values()))
        #NMI_list.append(NMI)      
        nb_run = nb_run +1
    
    
    Q_avg = communities.avg(Q_list)
    Q_max = communities.max(Q_list)
    Q_std = communities.stdev(Q_list)
    NMI_max = communities.max(NMI_list)
    time_run = communities.avg(Time_list)
    print("Q_avg",Q_avg,"Q_max",Q_max,"NMI",NMI_max,time_run)
    
if __name__ == '__main__':
    
    de_main()









    
                
