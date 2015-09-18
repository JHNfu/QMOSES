from __future__ import division
from __future__ import absolute_import
__author__ = 'John'
import scipy.sparse as sps
import pymetis as pm
from qmoseslib import *
import networkx as nx
import meshpy.triangle as triangle
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from six.moves import range
from cvxopt import spmatrix, amd

from qmoseslib import *

def needs_refinement(vertices, area):
    bary = np.sum(np.array(vertices), axis=0)/3
    max_area = 0.1 + (la.norm(bary, np.inf)-1)*0.1
    return bool(area > max_area)

# DOING IT ALL MANUALLY
# HOPE TO GET A SQUARE MESH WITHIN A BIGGER SQUARE
#
# # Points define the shape
# points_outer = [(-2,2),(2,2),(2,-2),(-2,-2)]
# points_inner = [(-1,1),(1,1),(1,-1),(-1,-1)]
# points = []
# points.extend(points_inner)
# points.extend(points_outer)
#
# # print 'Points:', points
# # print 'Length:', len(points)
#
# #Facets define the line between two points, denoted by the number of node as
# # outlined in the points vector
# facets_outer = [(4,5),(5,6),(6,7),(7,4)]
# facets_inner = [(0,1),(1,2),(2,3),(3,0)]
#
# #end = (2*len(points_outer))-1
# #start = len(points_inner)
# #facets_outer = [(i, i+1) for i in range(start, end)] + [(end, start)]
#
# facets = []
# facets.extend(facets_inner)
# facets.extend(facets_outer)
#
# meshinfo = triangle.MeshInfo()
# # Holes define which shape (series of faces creating a closed shape)
# meshinfo.set_holes([(0,0)])
# meshinfo.set_points(points)
# meshinfo.set_facets(facets)
#
# mesh = triangle.build(meshinfo, refinement_func=needs_refinement)
# #mesh = triangle.build(meshinfo)
#
# def round_trip_connect(start, end):
#     return [(i, i+1) for i in range(start, end)] + [(end, start)]
#
# mesh_points = np.array(mesh.points)
# mesh_tris = np.array(mesh.elements)
#
# # print '\nMesh Points NP Array:\n', mesh_points
# # print 'Type:', mesh.points
# # print '\nMesh Elements NP Array:\n', mesh_tris
# # print 'Type:', mesh.elements
#
# mesh_tris_list = mesh_tris.tolist()
# adjlist = meshpytrielements_to_adjlist(mesh_tris_list)
# print 'Final Adjacency List', adjlist
#
# import matplotlib.pyplot as pt
# # Triplot meshes triangles. Takes x_y mesh plots
# pt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)
# pt.show()

#adjlist = [[5, 6], [2,7], [1,9,10], [5,9], [7,10], [0, 3, 9, 8], [0,10], [9, 1, 10, 4], [5], [5,3,10,2,7], [6,9,2,7,4]]
#edgelist = adjlist_to_edgelist(adjlist)

test = nx.erdos_renyi_graph(15,0.2)
adjlist = test.adjacency_list()

# Pymetis nested dissection
n = 1
type_solver = 'isakov'
solver_address = r"C:\Users\John\PycharmProjects\public-isakov\public-isakov\build\visual_studio_2013\bin\an_ss_ge_fi_vdeg_omp.exe"

# FIRST IDENTIFY EDGE SEPARATORS
optimal_solution = recursive_bisection(n,adjlist, solver_type=type_solver, isakov_address=solver_address)
opt_list = qmosesnodelist_to_pymetisoutput(optimal_solution['nodes'])

edgelist = adjlist_to_edgelist(adjlist)
# removing duplicate edges
for edge in edgelist:
    if (edge[1],edge[0]) in edgelist:
        edgelist.remove((edge[1],edge[0]))

# FINDING EDGE SEPARATORS
# Identifying the Edges between partitions
edge_separators = []
possible_node_separators = []
for part_1 in range(2**n):
    for part_2 in range(2**n):
        if part_1 != part_2:
            for node_A in optimal_solution['nodes'][part_1][0]:
                for node_B in optimal_solution['nodes'][part_2][0]:
                    if (node_A, node_B) in edgelist:
                        edge_separators.append((node_A, node_B))
                        if node_A not in possible_node_separators:
                            possible_node_separators.append(node_A)
                        if node_B not in possible_node_separators:
                            possible_node_separators.append(node_B)
                    elif (node_A, node_B) in edgelist:
                        edge_separators.append((node_B, node_A))
                        if node_A not in possible_node_separators:
                            possible_node_separators.append(node_A)
                        if node_B not in possible_node_separators:
                            possible_node_separators.append(node_B)

print 'NODES A', optimal_solution['nodes'][0]
print 'NODES B', optimal_solution['nodes'][1]
print 'Edge Separators', edge_separators

# The node separators are found from the node cover of the edges. This is the
# vertex cover.
# THE VERTEX COVER CAN BE FOUND FROM THE COMPLEMENT OF THE MAXIMAL INDEPENDENT
# SET
G = nx.Graph()
G.add_edges_from(edge_separators)
max_ind_set = nx.maximal_independent_set(G)

N = nx.Graph()
N.add_nodes_from(nx.maximal_independent_set(G))
N.add_edges_from(edge_separators)
S = [x for x in possible_node_separators if x not in max_ind_set]

print 'Possible Node Separators', possible_node_separators
print 'Maximal Independent Set', nx.maximal_independent_set(G)
print 'Node separators:', S

# REMOVE NODE SEPARATORS FROM Q and P
Q = [x for x in optimal_solution['nodes'][0][0] if x not in S]
P = [x for x in optimal_solution['nodes'][1][0] if x not in S]

graph = nx.Graph()
graph.add_edges_from(edgelist)
graph.add_nodes_from(range(len(adjlist)))
adj_matrix = nx.adjacency_matrix(graph)
original_matrix = sps.csr_matrix(adj_matrix)

Q_subgraph = graph.subgraph(Q)
P_subgraph = graph.subgraph(P)
Q_adjmatrix = nx.adjacency_matrix(Q_subgraph)
P_adjmatrix = nx.adjacency_matrix(P_subgraph)

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

# NOW TO NEED FIND THE ORDERING OF P and Q using minimum tree
(Q_nonzero,Q_idx_loc,Q_jdx_loc) = adjmatrix_to_cvxformat(Q,Q_adjmatrix)
(P_nonzero,P_idx_loc,P_jdx_loc) = adjmatrix_to_cvxformat(P,P_adjmatrix)
Q_spmatrix = spmatrix(Q_nonzero,Q_idx_loc,Q_jdx_loc)
P_spmatrix = spmatrix(P_nonzero,P_idx_loc,P_jdx_loc)
Q_order_final = [x for x in amd.order(Q_spmatrix)]
P_order_final = [x for x in amd.order(P_spmatrix)]

# Change numbering of orderings in original_adjlist
# S becomes N - len(S) + Element In S
S_new = [(len(adjlist) - len(S) + x) for x in range(len(S))]
Q_new = Q_order_final
P_new = [x+len(Q) for x in P_order_final]

print 'Q_new: %s P_new: %s S_new: %s' % (Q_new, P_new, S_new)
print 'Q_old: %s P_old: %s S_old: %s' % (Q, P, S)

NodeList_old = Q + P + S
NodeList_new = Q_new + P_new + S_new

print NodeList_old
print NodeList_new

def update_adjlist(A_old,A_new,A_adjlist):
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