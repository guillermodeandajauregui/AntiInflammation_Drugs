import networkx as nx
from networkx.algorithms import bipartite
import pandas as pd

edges = [('A1','B1',3),
         ('A1','B2',7),
         ('A2','B1',2),
         ('A2','B2',4),
         ]

B = nx.Graph()
B.add_weighted_edges_from(edges)

def my_weight(G, u, v, weight='weight'):
  w = 0
for nbr in set(G[u]) & set(G[v]):
  w += G.edge[u][nbr].get(weight, 1) + G.edge[v][nbr].get(weight,1)
return w

G = bipartite.generic_weighted_projected_graph(B, ['A1', 'A2'], weight_function=my_weight)


print G.edges(data=True)
