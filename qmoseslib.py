##########################################################
# Q-MOSES - A D-WAVE Mesh Partitioner
##########################################################
# This is the staff for you to use, Moses.
# Not developed with a good wind and good timing, or even
# divine intervention, but rather an understanding of
# fluid dynamics, applied mathematics, quantum mechanics,
# & computational complexity theory.
##########################################################
# Developed as part of the Honours Project by John Furness
# With generous in-kind support provided by QxBranch and Shoal Group
# August 2015
##########################################################
# For algorithm explanations and use cases, see attached documentation.
# Written to PEP8 Python standard.
# Requires:
# numpy, networkx, matplotlib, pyplot, dwave_sapi
##########################################################
# The MIT License (MIT)
# Copyright (c) 2015 John Furness
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies
# or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
# AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##########################################################

import networkx as nx
import numpy as nm
import dwave_sapi
from dwave_sapi import local_connection, find_embedding

# GRAPH PARTITION FUNCTION
# Imports a variety of graph/mesh/adjacency files, you give the number of partitions,
# whether to use local or remote solver,
# Partitions into any number of partitions
# sub functions include Lucas Hamiltonian for 2^N partitions
# Your own Hamiltonian for prime numbers etc etc

# Function Splits Adjacency List Into Two Parts Based on the Hamiltonian
# Requires solver information

def partition(adjlist,solver):
    edgelist = adjlist_to_edgelist(adjlist)
    graph = nx.Graph()
    graph.add_edges_from(edgelist)

    max_degree = max(graph.degree().values())
    test_size = len(graph.nodes())
    num_vertices = len(adjlist)

    # Calculating A, assuming B = 1, ensuring a minimum of 1
    A = max(min(2*max_degree, num_vertices) / 8.0, 1)

    # Number of Qubits Required
    num_qubits = test_size - 1 # Minus 1 for degeneracy

    # Bias for each qubit
    h = [2*A] * num_qubits

    num_reads = 1000 # number of reads from D-Wave
    annealing_time  = 20

    # Collecting information on the available Chimera graph and number of qubits of the solver
    qubits = dwave_sapi.get_hardware_adjacency(solver)
    q_size = solver.properties["num_qubits"]

    J = dict()
    # Set up the default couplings to have a value of 2
    for i in range(num_qubits):
        for j in range(num_qubits):
            J[(i, j)] = 2*A

    # Set up the existing connections to have a value of 1.5
    for i,j in edgelist:
        # Deal with degeneracy by ignoring the last spin and setting it
        # up as a bias instead.
        if i == num_qubits:
            h[j] = 2*A-0.5 # For degeneracy
        elif j == num_qubits:
            h[i] = 2*A-0.5 # For degeneracy
        else:
            J[(i, j)] = 2*A-0.5

    # Fill in the Zeros
    for i in range(len(h)):
        for j in range(len(h)):
            if i >= j:
                J[(i, j)] = 0.0

    # Using the D-Wave API heuristic to find an embedding on the Chimera graph for the problem
    embeddings = dwave_sapi.find_embedding(J, len(h), qubits, q_size)

    # Sending problem to the solver
    embedding_solver = dwave_sapi.EmbeddingSolver(solver, embeddings)
    answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
    sols = answers['solutions']

    map(lambda sols: sols.append(1), answers['solutions'])
    opt_sol = sols[0]

    part_A_nodes = []
    part_B_nodes = []
    for idx, s in enumerate(opt_sol):
        if s == 1:
            part_A_nodes.append(idx)
        elif s == -1:
            part_B_nodes.append(idx)
        else:
            print 'ERROR: geez Louise - ' \
                  'You got something other than -1 or 1 in your sol.'

    # BUILDING EDGE LIST FOR COLOURED GRAPH
    graph_A = graph.subgraph(part_A_nodes)
    graph_B = graph.subgraph(part_B_nodes)

    return opt_sol

# A RECURSIVE BISECTION SCHEME FOR 2^N PARTITIONS
# Based on previously defined function, partition.
# INPUT recursive_bisection(n,adjlist,solver)
# N:        - 2^N partitions
# adjlist   - adjacency list of original graph
# solver    - D-Wave Solver
# OUTPUT list of nodes containing partition number
# def recursive_bisection(N,adjlist,solver):
#
#     edgelist = adjlist_to_edgelist(adjlist)
#     graph = nx.Graph()
#     graph.add_edges_from(edgelist)
#
#     if N == 1:
#         opt_sol = partition(adjlist,solver)
#     elif N > 2:
#         i = 2
#         for i in N:
#             opt_sol = partition(adjlist,solver)
#
#             part_A_nodes = []
#             part_B_nodes = []
#             for idx, s in enumerate(opt_sol):
#                 if s == 1:
#                     part_A_nodes.append(idx)
#                 elif s == -1:
#                     part_B_nodes.append(idx)
#                 else:
#                     print 'ERROR: geez Louise - ' \
#                           'You got something other than -1 or 1 in your sol.'
#
#             # BUILDING EDGE LIST FOR COLOURED GRAPH
#             graph_A = graph.subgraph(part_A_nodes)
#             graph_B = graph.subgraph(part_B_nodes)
#
#             adjlist_A = graph_A.adjacency_list()
#             adjlist_B = graph_B.adjacency_list()
#
#             opt_solA = partition(adjlist_A,solver)
#             opt_solB = partition(adjlist_B,solver)
#
#
#
#     else:
#         print 'ERROR: Enter value for N of one or greater'
#
#
#
#
#     return opt_sol

def adjlist_to_edgelist(adjlist):
    j = 0
    edgelist = []
    for y in adjlist:
        for i in y:
            edgelist.append((j, i))
        j += 1
    return edgelist


# FILL-REDUCING ORDERINGS FUNCTION FOR SPARSE MATRICES
# Using graph partitioning, or graph colouring methods, even clique methods

# MESSAGE SCHEDULING GRAPH COLOURING PROBLEM

# QUANTUM LATTICE GAS for FUN DEMONSTRATION perhaps