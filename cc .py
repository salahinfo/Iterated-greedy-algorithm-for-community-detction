
import networkx as nx


G = nx.read_gml('/home/yacine/Desktop/real_network/karate.gml')

average_degree = sum(dict(G.degree()).values()) / len(G)
print("Average Degree:", average_degree)

# Calculate the average clustering coefficient
average_clustering_coefficient = nx.average_clustering(G)
print("Average Clustering Coefficient:", average_clustering_coefficient)