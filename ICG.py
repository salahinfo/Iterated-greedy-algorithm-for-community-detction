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
        membership = { i : None for i in self.graph.nodes() }
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
                membership = super().insert_node( node , membership, com_id, comm_ngh.get(com_id, 0))
                                     
        
        return  membership 
    
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
                      
        merg_node = [ nod for nod in membership if nod not in index_community] 
        return  membership, index_community , merg_node
    
    def reconcstruction(self, membership, drop_node):
        membership = super().init(membership)
        random.shuffle(drop_node)
        for node in  drop_node:    
            comm_ngh = super().neigh_comm( node, membership)
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
                membership = super().insert_node( node, membership,  pos, comm_ngh.get( pos, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(membership.values()))
                membership = super().insert_node( node, membership, com_id, comm_ngh.get( com_id, 0.))

        return membership

    def crousel(self, membership,preserve_node, drop_node, alpha):
        membership = super().init(membership)
        iterations = int(alpha*self.n)
        for i in range(iterations):
            node = preserve_node.pop(0)
            drop_node.append(node)
            for vertex,com in membership.items():
                if node == vertex:
                    wgh = super().neigh_comm(node, membership)  
                    membership = super().delet_node( node, membership, membership[node], wgh.get( membership[node], 0.))
                    
            selected_node = random.choice(drop_node)
            drop_node.remove(selected_node)
            preserve_node.append(selected_node)
            MAX_Q = 0
            pos = -1    
            comm_ngh = super().neigh_comm( selected_node, membership)
            for com, Kbv in comm_ngh.items():    
                db = self.DegCom[com]
                delta_Q =  Kbv - self.Degree[selected_node]/(2.*self.m)*db
                if delta_Q > MAX_Q:
                    MAX_Q = delta_Q
                    pos = com
                else :
                    delta_Q = 0
            
            if MAX_Q > 0:
                if  membership[selected_node] != pos :
                    membership = super().insert_node( selected_node, membership, pos, comm_ngh.get( pos, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(membership.values()))
                membership = super().insert_node( selected_node, membership, com_id, comm_ngh.get( com_id, 0.))
        
        return membership, drop_node 
    
    def localsearch(self, membership):
        membership = super().init(membership)
        node_list = [i for i in self.G.nodes()]
        count = 0
        modified = True
        while modified:
            modified = False
            for vsele in node_list :
                degree = self.Degree[vsele]
                com_befor = membership[vsele]
                ngh_com = super().neigh_comm(vsele , membership)
                dvc = super().ngh_node(vsele, membership, com_befor)
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
                        
                membership = super().delet_node( vsele, membership,  com_befor, ngh_com.get( com_befor, 0.))            
                membership = super().insert_node( vsele, membership, m_com, ngh_com.get( m_com, 0.))   
           
        return membership
    
    def reconcstruction( self, membership, drop_node):
        membership = super().init(membership)
        random.shuffle( drop_node )
        for node in  drop_node:    
            comm_ngh = super().neigh_comm( node, membership)
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
                membership = super().insert_node( node, membership, pos, comm_ngh.get( pos, 0.))
            else:
                com_id = super().generate_random_not_in_list(set(membership.values()))
                membership = super().insert_node( node, membership, com_id, comm_ngh.get( com_id, 0.))

        return membership

                         
    def Run_ICG (self):
        start = time.time()
        soltion = self.GCH()
        soltion = self.localsearch(soltion)
        best_solution = copy.deepcopy(soltion)
        best_q = super().modularity(best_solution)
        print(best_q, time.time()-start)
        T_init = 0.025* best_q
        T = T_init
        nb_iter = 0
        while nb_iter < self.Nb:
            Q1 = super().modularity(soltion)
            incumbent_solution = copy.deepcopy( soltion)
            soltion,drop_nodes,preserve_node = self.Destruction(soltion)
            soltion, drop_nodes = self.crousel( soltion,preserve_node, drop_nodes, 0.6) 
            soltion = self.reconcstruction( soltion,drop_nodes)
            soltion = self.localsearch( soltion) 
            Q2 = super().modularity(soltion)
            print("q and time ",Q2, nb_iter)
            if Q2 > best_q:
                best_solution = copy.deepcopy( soltion)
                best_q = Q2
                
            P = random.random()
            if Q2 < Q1 and P > math.exp((Q2 - Q1)//T):
                soltion = copy.deepcopy(incumbent_solution)
                T = T_init

            elif Q2  >=  Q1:
                T = T_init
              
            else:
                T = T*0.9
        
            nb_iter = nb_iter + 1
        
        self.Mod_val = super().modularity(best_solution)
        end = time.time()
        t = end-start
        
        return best_q, best_solution,t                               

def de_main():
    path = sys.argv[1] 
    #'/home/yacine/Desktop/LFR/network_mu_0.8.dat'
    Number_iter = int(sys.argv[2])
    Beta = sys.argv[3]
    data = GraphTolls(path)
    graph = data.Read_Graph()
    NMI_list = []
    Time_list = [] 
    Q_list = []
    nb_run = 0
    while nb_run < int(sys.argv[5]) :
        communities = ICG( graph, Number_iter, Beta,path)
        mod, community, tim = communities.Run_ICG() 
        Q_list.append(mod)
        Time_list.append(tim)
        #label = communities.lebel_node(community)  
        if sys.argv[4]!= 'None':
            community = dict(sorted(community.items()))
            True_partition = data.Read_GroundTruth(sys.argv[4])
            NMI = normalized_mutual_info_score(True_partition,list(community.values()))
            NMI_list.append(NMI)      
        
        nb_run = nb_run +1
    
    #data.writefile(Q_list)
    Q_avg = communities.avg(Q_list)
    Q_max = communities.max(Q_list)
    Q_std = communities.stdev(Q_list)
    NMI_max = communities.max(NMI_list)
    time_run = communities.avg(Time_list)
    print("Q_avg",Q_avg,"Q_max",Q_max,"NMI",NMI_max,Q_std,time_run)
    
if __name__ == '__main__':
    
    de_main()









    
                
