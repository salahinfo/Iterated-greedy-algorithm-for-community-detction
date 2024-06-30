# The iterated greedy algorithm for community detection in complex networks:
# Introduction 
This repository includes the implementation of an iterated greedy algorithm [1] in Python to define the community structure in complex networks. The iterated greedy algorithm consists of two main components: destruction and reconstruction procedures 
the iterated algorithm consists of different steps:
```
1- initializeation procedures  
2- Iteratively aplly the destruction and reconstruction phases until no further improvment
```
## Hybrid carousell Iterated greedy algorithm (ICG) :
ICG [2] is considered an imporvment of IG algorithm, it integrate the local search procdure and the caroussel procedure to enhance the quality of solution 

# Usage :
 Note that the algorithm has been executed ten times to give you the metric value's maximum, average, and standard deviation. 
 Install all requirements packages mentioned in the file requirments.txt and download the datasets you want to apply the algorithm.
 Then, execute the below command line 
 if the networks have no ground truth. Execute this command.  
 ```
python3  IG.py  path_of_datasets  100  0.5  None  number of run algorihm
python3  ICG.py  path_of_datasets  100  0.5  None  number of run algorihm 

 ```
if the networks have ground-truth excute this commmnad 
```
python3  IG.py  path_of_datasets  100  0.5  path_of_ground-truth  number of run algorihm
python3  ICG.py  path_of_datasets  100  0.5  path_of_ground-truth  number of run algorihm 
```
# Refernce:
[1] : https://www.sciencedirect.com/science/article/abs/pii/S0167739X17323932
[2] : https://www.sciencedirect.com/science/article/abs/pii/S037843711931235X
