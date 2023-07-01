import math
import re
import sys 
import networkx as nx
import random

class GraphTolls:
     
    def Read_Graph(self, Path):
        
        if Path[len(Path)-3: ] == 'txt' or Path[len(Path)-3: ] == 'dat':
            Graph = nx.read_edgelist(Path, nodetype = int)
        elif Path[len(Path)-3: ] == 'gml':
            Graph = nx.read_gml(Path,label = 'id')
        else :
            raise TypeError (" the type of graph is not suportable or not no such file or directory")

        return Graph
    
    def Reve(self,x):
        sa = x.split()[::-1]
        l = []
        for i in sa:
            l.append(i)

        l=('  '.join(l))
        return l

    def Remove_Revers(self,path):
        with open(path, "r") as file:
            lines = file.readlines()
            result=[]
            for xa in lines:
                xa=re.sub(r'\s','  ',xa)
                if self.reve(xa) not in result:
                   result.append(xa.strip())   

        return(result)

    def Remove_Dublicate (self,path):
        pathw = path[ : -3]+'txt' 
        lines = self.remove_revers(path)
        with open(pathw, 'w') as f:
            for line in lines:
               f.write(line)
               f.write('\n')


    def Read_GroundTruth(self,path):
        with open(path, "r") as file:
            lines = file.readlines()
            result = []
            for x in lines:
                x = x.rstrip()
                result.append(x.split()[1])

        true_partion = [int(x)for x in result]
        return true_partion
    
     
    def Is_Intersiction(self,communities):
        dupes = []
        flat = [item for sublist in communities for item in sublist]
        for f in flat:
            if flat.count(f) > 1:
                if f not in dupes:
                    dupes.append(f)

        if dupes:
            return True
        else:
            return False   
        
    def sum(self,arg):
        if len(arg) < 1:
            return None
        else:
            return sum(arg)

    def count(self,arg):
        return len(arg)
  
    def min(self,arg):
        if len(arg) < 1:
            return None
        else:
            return min(arg)
  
    def max(self,arg):
        if len(arg) < 1:
            return None
        else:
            return max(arg)
  
    def avg(self,arg):
        if len(arg) < 1:
            return None
        else:
            return sum(arg) / len(arg)   
  
    def median(self,arg):
        if len(arg) < 1:
            return None
        else:
            arg.sort()
            return  arg[len(arg) // 2]
  
    def stdev(self,arg):
        if len(arg) < 1 or len(arg) == 1:
            return None
        else:
            avg = self.avg(arg)
            sdsq = sum([(i - avg) ** 2 for i in arg])
            stdev = (sdsq / (len(arg) - 1)) ** .5
            return stdev
  
    def percentile(self, arg):
        if len(arg) < 1:
            value = None
        elif (arg >= 100):
            sys.stderr.write('ERROR: percentile must be < 100.  you supplied: %s\n'% arg)
            value = None
        else:
            element_idx = int(len(arg) * (arg / 100.0))
            self.arg.sort()
            value = self.arg[element_idx]
        return value  
    
    def select_edge_betw(self,g,* arg):
        Edg_betw = 0
        for i in list(g.neighbors(arg[0])):
            if i in arg[1] :
                Edg_betw += 1
        
        
        return Edg_betw
    def is_edge_betw(self,g,vert,commu):
        Edg_betw = []
        for i in list(g.neighbors(vert)):
            if i in commu :
                return True
            
        return False
   
