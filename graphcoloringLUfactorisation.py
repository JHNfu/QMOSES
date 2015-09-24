from qmoseslib import *
import networkx as nx
from dwave_sapi import local_connection
import matplotlib.pyplot as plt
import scipy.sparse as sps

no_nodes = 10
no_col = 4 # Number of colours

graph = nx.erdos_renyi_graph(no_nodes, 0.3)
pos1 = nx.spring_layout(graph)

edgelist = graph.edges()
adjlist = graph.adjacency_list()
adj_matrix = nx.adjacency_matrix(graph)
original_matrix = sps.csr_matrix(adj_matrix)

# solver = local_connection.get_solver("c4-sw_optimize")
# opt_sol = colourgraph(no_col,adjlist, solver_type='dwave', dwave_solver = solver)

# print opt_sol, '\n'

solver_address = r"C:\Users\John\PycharmProjects\public-isakov\public-isakov\build\visual_studio_2013\bin\an_ss_ge_fi_vdeg_omp.exe"
opt_sol = colourgraph(no_col, adjlist, solver_type='isakov', isakov_address = solver_address)

print opt_sol, '\n'

subgroups = [[] for i in range(no_col+1)]
for idx, result in enumerate(opt_sol):
    subgroups[result].append(idx)

colours = ['#ffc0cb', '#0000ff', '#ff0000', '#40e0d0', '#800080', '#f6546a', '#ffa500', '#088da5', '#008000']

plt.figure(1)
plt.subplot2grid((2,1),(0,0), colspan = 1)

for idx in range(no_col + 1):
    nx.draw_networkx_nodes(graph,pos1,
                           nodelist=subgroups[idx],
                           node_color=colours[idx],
                           node_size=500,
                           alpha=1)

nx.draw_networkx_edges(graph,pos1,
                       edgelist = edgelist)
labels=nx.draw_networkx_labels(graph,pos=pos1,
                               font_size=15,
                               font_color='w',
                               font_family='sans-serif',)
plt.subplot2grid((2,1),(1,0))
plt.spy(original_matrix, markersize=14)

plt.show()
