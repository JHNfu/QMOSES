from __future__ import division
from __future__ import absolute_import
__author__ = 'John'
import scipy.sparse as sps
import matplotlib.pyplot as plt
from six.moves import range
from cvxopt import spmatrix, amd
from networkx.algorithms import bipartite

from qmoseslib import *

test = nx.erdos_renyi_graph(20,0.5)
adjlist = test.adjacency_list()
#adjlist = [[12, 5], [12, 2, 4], [1, 9], [12, 4, 6], [8, 1, 3, 11], [0, 9], [11, 3], [10, 13, 14], [4, 14], [2, 5], [11, 12, 7], [10, 4, 14, 6], [0, 1, 10, 3], [14, 7], [8, 11, 13, 7]]
#print adjlist

# Pymetis nested dissection
n = 1
type_solver = 'isakov'
solver_address = r"C:\Users\John\PycharmProjects\public-isakov\public-isakov\build\visual_studio_2013\bin\an_ss_ge_fi_vdeg_omp.exe"

# FIRST IDENTIFY EDGE SEPARATORS
optimal_solution = recursive_bisection(n,adjlist, solver_type=type_solver, isakov_address=solver_address)

optimal_solution_partition = partition(adjlist,solver_type=type_solver, isakov_address=solver_address)
optimal_solution = optimal_solution_partition
print 'partition', optimal_solution_partition

print optimal_solution_partition['node_list']

node_list = optimal_solution_partition['node_list']
opt_list = qmosesnodelist_to_pymetisoutput(node_list)

# Finding edge separator Separators
edgelist = adjlist_to_edgelist(adjlist)
edge_separators, possible_node_separators = find_edgeseparators(node_list,edgelist)

# Finding node separator from edge separator using networkx function
G = nx.Graph()
G.add_edges_from(edge_separators)
print 'Is the graph bipartite?'
print(bipartite.is_bipartite(G))
matching = nx.bipartite.maximum_matching(G)
vertex_cover = list(nx.bipartite.to_vertex_cover(G, matching))
# print 'edge separators:', edge_separators
# print 'vertex cover:', vertex_cover
S = vertex_cover

# calculating quantum vertex cover
edge_sep_adjlist = edgelist_to_adjlist(edge_separators)
h, J, offset = vertexcover(edge_sep_adjlist)
answer = isakovlib.solve_Ising(h, J, offset)
opt_sol = answer['solutions'][0]

vertex_cover_qc = []
for idx, node in enumerate(opt_sol):
    if node == 1:
        vertex_cover_qc.append(idx)

print 'Possible Node Separators', possible_node_separators
print 'Quantum Node Separator', vertex_cover_qc
print 'Node separators:', S

# REMOVE NODE SEPARATORS FROM Q and P
Q = [x for x in node_list[0] if x not in S] #tick
P = [x for x in node_list[1] if x not in S] #tick

print 'Q', Q
print 'P', P

graph = nx.Graph()
graph.add_edges_from(edgelist)
graph.add_nodes_from(range(len(adjlist)))
adj_matrix = nx.adjacency_matrix(graph)
original_matrix = sps.csr_matrix(adj_matrix)

if S == node_list[0]:
    print 'ERROR'
    exit()


# Building subgraph based on Q and P nodes
Q_subgraph = graph.subgraph(Q)
P_subgraph = graph.subgraph(P)
Q_adjmatrix = nx.adjacency_matrix(Q_subgraph)
P_adjmatrix = nx.adjacency_matrix(P_subgraph)

print 'PSUBGRAPH'
print P_subgraph

print 'PMATRIX'
print P_adjmatrix
print type(P_adjmatrix)
print P_subgraph.nodes()
print 'QMATRIX'
print Q_adjmatrix
print Q_subgraph.nodes()

#converting sparse matrix dict to cxvopt format
def adjmatrix_to_cvxformat(node_list,adjmatrix):
    nonzero = []
    idx_loc = []
    jdx_loc = []
    for idx in range(len(node_list)):
        for jdx in range(len(node_list)):
            if adjmatrix[(idx,jdx)] == 1:
                nonzero.append(adjmatrix[(idx,jdx)])
                idx_loc.append(idx)
                jdx_loc.append(jdx)
    return(nonzero, idx_loc, jdx_loc)

# NOW TO NEED FIND THE ORDERING OF P and Q using minimum spanning tree to improve results
(Q_nonzero,Q_idx_loc,Q_jdx_loc) = adjmatrix_to_cvxformat(Q,Q_adjmatrix)
(P_nonzero,P_idx_loc,P_jdx_loc) = adjmatrix_to_cvxformat(P,P_adjmatrix)
Q_spmatrix = spmatrix(Q_nonzero,Q_idx_loc,Q_jdx_loc)
P_spmatrix = spmatrix(P_nonzero,P_idx_loc,P_jdx_loc)
Q_order_final = [x for x in amd.order(Q_spmatrix)]
P_order_final = [x for x in amd.order(P_spmatrix)]

print amd.order(Q_spmatrix)
print 'PMATRIXX'
print P_spmatrix
print 'QMATRIXX'
print Q_spmatrix

# Change numbering of orderings in original_adjlist
# S becomes N - len(S) + Element In S
S_new = [(len(adjlist) - len(S) + x) for x in range(len(S))]
Q_new = Q_order_final
print 'P_order_dinal', P_order_final
print 'Q_order_dinal', Q_order_final
P_new = [x+len(Q) for x in P_order_final]

print 'Q_new: %s P_new: %s S_new: %s' % (Q_new, P_new, S_new)
print 'Q_old: %s P_old: %s S_old: %s' % (Q, P, S)

NodeList_old = Q + P + S
NodeList_new = Q_new + P_new + S_new

print NodeList_old
print NodeList_new

print 'A_old',

def update_adjlist(A_old,A_new,A_adjlist):
    '''
    reorders adjacency list based on a new and old nodelist
    :param A_old:
    :param A_new:
    :param A_adjlist:
    :return:
    '''
    new_adjlist = list(A_adjlist)
    for ndx, node in enumerate(A_adjlist):
        for adjndx, adjnode in enumerate(node):
            for idx, sep_node in enumerate(A_old):
                #print node, adjnode, idx, sep_node
                if adjnode == sep_node:
                    #print 'well then better change %s to a %s' % (adjlist[ndx][adjndx],A_new[idx])
                    new_adjlist[ndx][adjndx] = A_new[idx]
    # Now reordering list of lists
    for ndx, node in enumerate(A_adjlist):
        for idx, sep_node in enumerate(A_old):
            if ndx == sep_node:
                new_adjlist[A_new[idx]] = node
    return new_adjlist

print 'Old_Adjlist List', adjlist
print len(adjlist)

final_adjlist = update_adjlist(NodeList_old,NodeList_new,adjlist)

print 'New_Adjlist List', final_adjlist
# P_new =
# [x for element in [y for node in adjlist]]
# Q_new = [x for x in range(len(Q))]

nodesep = nx.Graph()
final_edgelist = adjlist_to_edgelist(final_adjlist)
nodesep.add_edges_from(final_edgelist)
nodesep.add_nodes_from(range(len(adjlist)))
adj_matrix = nx.adjacency_matrix(nodesep)
nodesep_matrix = sps.csr_matrix(adj_matrix)
edgesep_matrix = original_matrix
# nodes
plt.figure(1)
plt.subplot2grid((2,3),(0,0), colspan = 3, rowspan = 1)
plt.axis('off')
#pos = nx.spring_layout(graph)
pos = nx.shell_layout(graph)
# pos = nx.graphviz_layout(G, prog='neato', root=None, args='')
plt.title('Fill-reducing Ordering Based on Graph Partition')
nx.draw_networkx_nodes(graph,pos,
                       nodelist=Q,
                       node_color='r',
                       node_size=500,
                       alpha=1)
nx.draw_networkx_nodes(G,pos,
                       nodelist=P,
                       node_color='b',
                       node_size=500,
                       alpha=1)
nx.draw_networkx_nodes(G,pos,
                       nodelist=S,
                       node_color='k',
                       node_size=500,
                       alpha=1)
nx.draw_networkx_labels(G,pos,font_color='w')
nx.draw_networkx_edges(G,pos,
                       edgelist=edge_separators,
                       width=8,alpha=0.8,edge_color='y')
nx.draw_networkx_edges(G,pos,
                       edgelist=edgelist,
                       width=1,alpha=1,edge_color='k')

### PLOT SO BLACK CIRCLES, it'll look way better
### SUB PLOT BEFORE AND AFTER
plt.subplot2grid((2,2),(1,0))
plt.axis('off')
plt.title('Original')
plt.spy(original_matrix , markersize=14) ### THIS ONE IS CLEARLY BEFORE

# plt.subplot2grid((2,3),(1,1))
# plt.title('B')
# plt.axis('off')
# plt.spy(edgesep_matrix, markersize=14) ### THIS ONE IS CLEARLY BEFORE

plt.subplot2grid((2,2),(1,1))
plt.title('Reduced Fill-in')
plt.axis('off')
plt.spy(nodesep_matrix, markersize=14) ### THIS ONE IS CLEARLY BEFORE
# NOW TO REORDER MATRIX

plt.show()




# Nested dissection algorithm using partition.