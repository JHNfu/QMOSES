####################################
# MESH PARTITIONING HAMILTONIAN
# Generates symbolic function for Mesh Partitioning Hamiltonian.
# John Furness, September 2015
####################################
from __future__ import division
from __future__ import absolute_import
import meshpy.triangle as triangle
import numpy as np
import numpy.linalg as la
from six.moves import range
from sympy import *
import networkx as nx
__author__ = 'John Furness'
#####################################
# PROBLEM VARIABLES
#####################################
no_part = 3

#####################################
# GENERATING A MESH
#####################################
def needs_refinement(vertices, area):
    bary = np.sum(np.array(vertices), axis=0)/3
    max_area = 0.8 + (la.norm(bary, np.inf)-1)*0.6
    return bool(area > max_area)

# Points define the shape
points_outer = [(-2,2),(2,2),(2,-2),(-2,-2)]
points_inner = [(-1,1),(1,1),(1,-1),(-1,-1)]
points = []
points.extend(points_inner)
points.extend(points_outer)

#Facets define the line between two points, denoted by the number of node as
# outlined in the points vector
facets_outer = [(4,5),(5,6),(6,7),(7,4)]
facets_inner = [(0,1),(1,2),(2,3),(3,0)]

facets = []
facets.extend(facets_inner)
facets.extend(facets_outer)

meshinfo = triangle.MeshInfo()
# Holes define which shape (series of faces creating a closed shape)
meshinfo.set_holes([(0,0)])
meshinfo.set_points(points)
meshinfo.set_facets(facets)

mesh = triangle.build(meshinfo, refinement_func=needs_refinement)
#mesh = triangle.build(meshinfo)

def round_trip_connect(start, end):
    return [(i, i+1) for i in range(start, end)] + [(end, start)]

mesh_points = np.array(mesh.points)
mesh_tris = np.array(mesh.elements)

mesh_points_list = mesh_points.tolist()
mesh_tris_list = mesh_tris.tolist()
adjlist = meshpytrielements_to_adjlist(mesh_tris_list)
print 'Final Adjacency List', adjlist
no_vert = len(adjlist)

# FOR PLOTTING MESH
# import matplotlib.pyplot as pt
# Triplot meshes triangles. Takes x_y mesh plots
# pt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)

# DETERMINING POSITION OF CENTRE OF ELEMENTS FOR NETWORKX
element_centres = []
position = dict()
for element_no, element in enumerate(mesh_tris_list):
    element_x = 0
    element_y = 0
    for node in element:
        element_x += mesh_points_list[node][0]
        element_y += mesh_points_list[node][1]
    element_centres.append([element_x/3, element_y/3])
    position[element_no] = np.array([element_x/3, element_y/3])

# SETTING UP GRAPH OBJECT FOR NETWORKX
edgelist = adjlist_to_edgelist(adjlist)
graph = nx.Graph()
graph.add_edges_from(edgelist)
graph.add_nodes_from(range(len(adjlist)))

# PRINTING DATA
print 'Number of nodes:', no_vert
print 'Number of partition:', no_part
print 'Number of edges:', len(edgelist)
print 'Therefore need %s qubits.\n' % (no_vert*no_part)

#####################################
# SETTING UP HAMILTONIAN USING SYMPY
#####################################
# DEFINING SYMBOLIC FUNCTIONS
x = symbols('x(1:%d\,1:%d)' % (((no_vert + 1) , (no_part + 1))))
H1 = Symbol('H1')
H2 = Symbol('H2')
H3 = Symbol('H3')
A = Symbol('A')
B = Symbol('B')
C = Symbol('C')

# TERM 1: PENALTY TERM TO ENSURE EVERY NODE IS GIVEN A PARTITION
H1 = 0
term = 1
for v in range(no_vert):
    for i in range(no_part):
        v += 1 # FIX THIS!
        i += 1 # FIX THIS!
        term -= x[((v - 1) * no_part) + i - 1]
        v -= 1 # FIX THIS!
        i -= 1 # FIX THIS!
    term = term**2
    H1 += term
    term = 1

# TERM 2: PENALTY TERM TO ENSURE SAME NO. OF NODES IN EACH PARTITION
H2 = 0
term = no_vert/no_part
for v in range(no_part):
    for i in range(no_vert):
        term -= x[((i - 1) * no_part) + v - 1]
    term = term**2
    H2 += term
    term = no_vert/no_part

# Combining Penalty Terms
HA = A*H1 + B*H2

# TERM 3: OPTIMISER TERM If nodes sharing an edge are in different partitions,
# a penalty is given. Vice vera, if nodes in the same partition and share
# an edge, the function is minimised.
H3 = 0
for edge in edgelist:
    #print edge
    for idx in range(no_part):
        for jdx in range(idx+1,no_part):
                H3 += x[((edge[0]*no_part) + idx)]*x[((edge[1]*no_part) + jdx)]

# Combining Optimiser Terms
HB = C*H3
# Combining to Final Hamiltonian
H = HA + HB

# For substituting coefficients
# H = H.subs('A', A_subs)
# H = H.subs('B', B_subs)
# H = H.subs('C', C_subs)

# Expanding Hamiltonian to eventually yield coefficients
H = H.expand(H)

# Accounting for Degeneracy, assume first node is in partition 1
H = H.subs(x[0], 1)
for s in range(no_part):
    H = H.subs(x[s],0)
H = H.expand(H)

#####################################
# CREATING Q in QUBO
#####################################
Q = dict()
for idx in range(len(x)):
    for jdx in range(len(x)):
        Q[(idx, jdx)] = 0
# ALSO FINDS MAXIMUM COEFFICIENT for NORMALISATION
print 'Printing Q:'
max_coeff = 0
for idx in range(len(x)):
    print ' '
    for jdx in range(len(x)):
        if idx != jdx:
            #Q[(idx, jdx)] = round(H.coeff(x[idx]*x[jdx]),2)
            Q[(idx, jdx)] = H.coeff(x[idx]*x[jdx])
            #maxtest = round(H.coeff(x[idx]*x[jdx]),2)
            maxtest = H.coeff(x[idx]*x[jdx])
            print maxtest, '\t',
        else:
            E = H.coeff(1*x[idx])
            #Q[(idx,jdx)] = round((H.coeff(x[idx]**2) - E.coeff(-1)),2)
            Q[(idx,jdx)] = (H.coeff(x[idx]**2) - E.coeff(-1))
            #maxtest = round((H.coeff(x[idx]**2) - E.coeff(-1)),2)
            maxtest = (H.coeff(x[idx]**2) - E.coeff(-1))
            print maxtest, '\t',

        #if abs(maxtest) > max_coeff:
        #    max_coeff = maxtest

#print max_coeff