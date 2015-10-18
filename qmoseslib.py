from __future__ import division
##############################################################################
# +++++++++++++++++++ Q-MOSES - A D-WAVE Mesh Partitioner +++++++++++++++++++
##############################################################################
# This is the staff for you to use, Moses.
# Not developed with a good wind and good timing, or even
# divine intervention, but rather an understanding of
# fluid dynamics, applied mathematics, quantum mechanics,
# & computational complexity theory.
##############################################################################
# Developed as part of the Honours Project by John Furness
# With generous in-kind support provided by QxBranch and Shoal Group
# August 2015
##############################################################################
# For algorithm explanations and use cases, see attached documentation.
# Written to PEP8 Python standard.
# Requires:
# numpy, networkx, matplotlib, pyplot, dwave_sapi
##############################################################################
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
##############################################################################

import scipy.linalg as sl
import scipy.sparse as sp
from random import randint
import networkx as nx
from networkx.algorithms import bipartite
import numpy as np
import nxmetis
import dwave_sapi
import time
from dwave_sapi import local_connection, find_embedding
from collections import defaultdict
import copy as copy
from sympy import *
import isakovlib

# tristan checked
##############################################################################
# VARIOUS FUNCTIONS FOR MESH PARTITIONING
##############################################################################


def mosespartitioner_crossgrade(no_part, adjlist, load = None, solver_type = None, dwave_solver = None, isakov_address = None):
    """

    In hours of attempts to ensure, mosespartitioner was written correctly,
    the following was written to improve readability of the code.

    Each term is expanded into a Q and then added together at the end.

    This includes changes accounting for degeneracy (ie. that one element is
    allocated to partition one).

    Ultimately, the original function and this one yield the same solutions
    for the terms; ensuring balance of nodes, and term ensure each node is
    allocated a partition.

    The difference occurs in ensuring the edge boundary
    is minimised. When accounting for degeneracy there are some diagonal terms
    that don't occur in the original function.

    The additional change is that the LoadVectors are appropriately applied
    to the their respectively partitions, ie. for node 2 (x2) for three
    partitions: X2,1 is affect by LoadVector1, and X2,2, is affected by
    LoadVector2 etc

    NOTE ON COEFFICIENTS:
    A : Attributed to term that ensures element balance
    B : Attributed to term that ensures one in every element
    C : Attributed to term that minimises edgeboundary

    Continuation of previous description..
    Moses Partitioner: implementation of the flagship Hamiltonian for this
    library of impeccable functions.
    Moses partitioner can do a  lot, you know? Partition into a positive
    integer number of parts, which is real handy. Additionally it can be used
    for Load Balancing on heterogeneous systems.

    Requires (Num_Partitions - 1) * Num_Elements fully connected qubits.

    :param no_part: number of partitions
    :param adjlist: adjacency list
    :param load: load vector consisting of Num_Partitions elements and summing
    to one
    :param solver_type: string or 'isakov' or 'dwave'
    :param dwave_solver: solver class, local connection etc
    :param isakov_address: solver address e.g., '= r"C:..........'
    :return: well, quite a lot really.
    """

    no_vert = len(adjlist)
    no_qubits = no_vert*no_part

    print 'Number of Elements', no_vert

    #### adjust load vector so that one partition yields an extra node
    if no_vert % no_part == 1 and load == None:

        # adjust load vector so that one partition yields an extra node

        LoadVector = [i * (no_vert+1)/(no_vert*no_part) for i in [1]*(no_part-1)]
        LoadVector.append(1 - sum(LoadVector))

    elif load == None:
        LoadVector = [i * (1/no_part) for i in [1]*no_part]
        print 'Assuming homogenous system...'

    elif len(load) == no_part and sum(load) == 1:
        print 'Partitioning heterogenous system...'
        LoadVector = load
        print 'Load Vector:', LoadVector

    else:
        print 'Moses Partitioner Error: Load Vector is incorrect type.' \
              '\t Load Vector should be a list of length(no. of partitions) and sum to 1.'

    print 'LoadVector:', LoadVector

    # MESH PARTITIONING HAMILTONIAN HARD-CODE
    edgelist = adjlist_to_edgelist(adjlist)
    no_edges = len(edgelist)
    print 'Number of edges', no_edges
    graph = nx.Graph()
    graph.add_edges_from(edgelist)
    graph.add_nodes_from(range(len(adjlist)))

    # NOTE ON COEFFICIENTS:
    # A : Attributed to term that ensures element balance
    # B : Attributed to term that ensures one in every element
    # C : Attributed to term that minimises edgeboundary

    B = 2 # Ensure everyone gets a partition
    A = 1 # Correct share of nodes
    C = 1 # Minimise edge boundary

    print 'Moses Partitioner Coefficients: A = %s, B = %s, C = %s' % (A, B, C)

    # SETTING UP Q DICT
    Q1 = dict()
    for idx in range(no_qubits):
        for jdx in range(no_qubits):
            Q1[(idx, jdx)] = 0

    # TERM 1: ENSURE CORRECT NODE SHARE IN LOAD BALANCE

    Q1 = dict()
    for idx in range(no_qubits):
        for jdx in range(no_qubits):
            Q1[(idx, jdx)] = 0

    # The off diagonals
    for vert_ydx in range(no_vert):
        for partdx in range(no_part):
            #print partdx
            for vert_xdx in range(no_vert):
                jdx = (vert_ydx * no_part) + partdx
                idx = (vert_xdx * no_part) + partdx
                if idx > jdx:
                    Q1[(jdx,idx)] = 2*A

    # The diagonals
    for vertdx in range(no_vert):
        for partdx in range(no_part):

            idx = (vertdx*no_part) + partdx

            if partdx == 0:
                # Energy on first partition is higher due to degeneracy
                Q1[(idx,idx)] = 1*A - 2*no_vert*LoadVector[partdx]*A + 2*A # 1A for the square terms
            else:
                # Energy on remaining partitions is lower due to degeneracy
                Q1[(idx,idx)] = (1*A) - 2*no_vert*LoadVector[partdx]*A # 1A for the square terms

    # print "TERM 1 EVEN NODE SHARE:"
    # for idx in range(no_qubits):
    #     print ' '
    #     for jdx in range(no_qubits):
    #         print Q1[(idx, jdx)],

    print "\n"

    # TERM 2: ENSURE EVERY NODE GETS A PARTITION
    Q2 = dict()
    for idx in range(no_qubits):
        for jdx in range(no_qubits):
            Q2[(idx, jdx)] = 0

    # The diagonals
    for idx in range(no_qubits):
        Q2[(idx, idx)] = -2*B + 1*B # 1*B is for (x1)^2 terms

    # might be possible solution without itertools
    from itertools import combinations
    for vertdx in range(no_vert):
        for Q2_key in list(combinations(range(vertdx*no_part, (vertdx+1)*no_part), 2)):
                Q2[Q2_key] = 2*B

    # print "Printing Q2: everynode gets a partition"
    # for idx in range(no_qubits):
    #     print ''
    #     for jdx in range(no_qubits):
    #         print Q2[(idx, jdx)],

    # TERM 3: ENSURE EDGE BOUNDARY IS MINIMISED
    Q3 = dict()
    for idx in range(no_qubits):
        for jdx in range(no_qubits):
            Q3[(idx, jdx)] = 0


    combo = list(combinations(range(0, no_part), 2))
    print combo
    for i, j in edgelist:
        #print "i %s, j %s" % (i, j)
        for idx, jdx in combo:
            if i == 0:
                Q3_jdx = jdx + (j*no_part)
                Q3[(Q3_jdx, Q3_jdx)] = 1*C
            elif j == 0:
                Q3_idx = idx + (i*no_part)
                Q3[(Q3_idx, Q3_idx)] = 1*C

    for i, j in edgelist:
        print i, j
        #print "i %s, j %s" % (i, j)
        if j > i:
            for part_i in range(no_part):
                for part_j in range(part_i+1,no_part):
                    Q3_idx = (i*no_part) + part_i
                    Q3_jdx = (j*no_part) + part_j
                    Q3[(Q3_idx, Q3_jdx)] = 1

            for part_i in range(no_part):
                for part_j in range(part_i+1,no_part):
                    Q3_jdx = (j*no_part) + part_i
                    Q3_idx = (i*no_part) + part_j
                    Q3[(Q3_idx, Q3_jdx)] = 1*C
    # print edgelist
    # print "Printing Q3"
    # for idx in range(no_qubits):
    #     print ''
    #     for jdx in range(no_qubits):
    #         print Q3[(idx, jdx)],

    # BRINGIN' IT ALL TOGETHER
    Q_final = dict()
    for idx in range(no_qubits):
        print ''
        for jdx in range(no_qubits):
            Q_final[(idx, jdx)] = Q1[(idx, jdx)] + Q2[(idx, jdx)] + Q3[(idx, jdx)]
            print Q_final[(idx, jdx)],

    # REMOVE DEGENERATE SPIN
    Q_deg = dict()
    for idx in range(no_part, no_vert*no_part):
        # print ''
        for jdx in range(no_part, no_vert*no_part):
            # print '%.f' % (Q[(idx,jdx)]),
            Q_deg[(idx - no_part,jdx - no_part)] = (Q_final[(idx,jdx)])

    # print "\n\nQ with degeneracy:"
    # for idx in range(no_qubits-no_part):
    #     print ''
    #     for jdx in range(no_qubits-no_part):
    #         print Q_deg[(idx, jdx)],

    (h, J, offset) = dwave_sapi.qubo_to_ising(Q_deg)

    num_reads = 1000

    if solver_type == 'dwave':
        qubits = dwave_sapi.get_hardware_adjacency(dwave_solver)
        q_size = dwave_solver.properties["num_qubits"]

        # Using the D-Wave API heuristic to find an embedding on the Chimera graph for the problem
        embeddings = dwave_sapi.find_embedding(J, len(h), qubits, q_size)

        # Sending problem to the solver
        embedding_solver = dwave_sapi.EmbeddingSolver(dwave_solver, embeddings)
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        sols = answers['solutions']

        t0 = time.time()
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        t1 = time.time()
        answers['timing'] = t1 - t0

        answers['num_reads'] = num_reads

        # adding offset from Q to h, J conversion back to energies
        answers['energies'] = [energy + offset for energy in answers['energies']]
        answers['hvalues'] = h
        answers['jvalues'] = J
        answers['offset'] = offset

    elif solver_type == 'isakov':

        answers = isakovlib.solve_Ising(h, J, offset = offset)

        # adding degeneracy
        negative = '-'
        deg_plusminus = '+' + negative*(no_part-1)
        answers['plusminus'] = [deg_plusminus + sol for sol in answers['plusminus']]

    else:
        print "Moses Partitioner Error: Solver type not recognised. Try dwave or isakov"

    # adding degeneracy back in
    deg = [1] + [-1]*(no_part-1)
    answers['solutions'] = [deg + sol for sol in answers['solutions']]

    # Creating Pymetis style node output based on optimal solution
    opt_sol = answers['solutions'][0]
    print "Optimal Solution", opt_sol

    final_solution = []
    incorrect = 0
    for nodes in range(no_vert):
        colour_sols = []
        for idx in range(nodes*no_part,(nodes*no_part)+no_part):
            colour_sols.append(opt_sol[idx])

        if sum(colour_sols) == (no_part*-1)+2:
            for idx in range(nodes*no_part,(nodes*no_part)+no_part):
                if opt_sol[idx] == 1:
                    final_solution.append(idx - nodes*no_part + 1)
                    # print final_solution
        else:
            final_solution.append(0)
            incorrect += 1

    print "Final_Solution", final_solution

    answers['nodelist'] = final_solution
    answers['adjlist'] = adjlist
    answers['num_parts'] = no_part
    answers['load_vector'] = LoadVector
    answers['solver'] = solver_type
    answers['Q'] = Q
    answers['Qdeg'] = Q_deg
    answers['A'] = A
    answers['B'] = B
    answers['C'] = C

    return answers


def mosespartitioner(no_part, adjlist, load = None, solver_type = None, dwave_solver = None, isakov_address = None):
    """

    Moses Partitioner: implementation of the flagship Hamiltonian for this
    library of impeccable functions.
    Moses partitioner can do a  lot, you know? Partition into a positive
    integer number of parts, which is real handy. Additionally it can be used
    for Load Balancing on heterogeneous systems.

    Requires (Num_Partitions - 1) * Num_Elements fully connected qubits.

    :param no_part: number of partitions
    :param adjlist: adjacency list
    :param load: load vector consisting of Num_Partitions elements and summing
    to one
    :param solver_type: string or 'isakov' or 'dwave'
    :param dwave_solver: solver class, local connection etc
    :param isakov_address: solver address e.g., '= r"C:..........'
    :return: well, quite a lot really.
    """

    no_vert = len(adjlist)

    print 'Number of Elements', no_vert

    #### adjust load vector so that one partition yields an extra node
    if no_vert % no_part == 1 and load == None:

        # adjust load vector so that one partition yields an extra node

        LoadVector = [i * (no_vert+1)/(no_vert*no_part) for i in [1]*(no_part-1)]
        LoadVector.append(1 - sum(LoadVector))

    elif load == None:
        LoadVector = [i * (1/no_part) for i in [1]*no_part]
        print 'Assuming homogenous system...'

    elif len(load) == no_part and sum(load) == 1:
        print 'Partitioning heterogenous system...'
        LoadVector = load
        print 'Load Vector:', LoadVector

    else:
        print 'Moses Partitioner Error: Load Vector is incorrect type.' \
              '\t Load Vector should be a list of length(no. of partitions) and sum to 1.'

    print 'LoadVector:', LoadVector

    # MESH PARTITIONING HAMILTONIAN HARD-CODE
    edgelist = adjlist_to_edgelist(adjlist)
    no_edges = len(edgelist)
    print 'Number of edges', no_edges
    graph = nx.Graph()
    graph.add_edges_from(edgelist)
    graph.add_nodes_from(range(len(adjlist)))

    # for larger mesh
    A = 0.5*no_edges # Ensure everyone gets a partition
    B = 1 # Correct share of nodes
    C = 1 # Minimise edge boundary

    # A: Ensure everyone gets a partition
    # B: Correct share of nodes
    # C: Minimise edge boundary

    # Following coefficients worked superbly on small graphs
    # A = 0.5*no_edges # ensures each has a partition
    # B = 1 # ensures equal nodes
    # C = 1 # ensures edge bondary minimised.

    # # Following coefficients worked superbly on small graphs
    # A = 0.5*no_edges # ensures each has a partition
    # B = 1 # ensures equal nodes
    # C = 0.1*no_vert  # ensures edge bondary minimised.

    # # Following coefficients worked superbly on small graphs
    # A = 0.5*no_edges # ensures each has a partition
    # B = 1 # ensures equal nodes
    # C = 0.3*no_vert  # ensures edge bondary minimised.

    print 'Moses Partitioner Coefficients: A = %s, B = %s, C = %s' % (A, B, C)

    Q_HC = dict()
    for idx in range(no_part*no_vert):
        for jdx in range(no_part*no_vert):
            Q_HC[(idx, jdx)] = 0

    for idx in range(no_vert):
        for jdx in range(no_vert):
            for print_idx in range(idx*no_part,idx*no_part + no_part):
                # print ' '
                for print_jdx in range(jdx*no_part,jdx*no_part + no_part):
                    if (idx, jdx) in edgelist:
                        Q_HC[(print_idx, print_jdx)] = C
                    # print print_idx, print_jdx, idx*no_part, jdx*no_part
                    if (print_idx - idx*no_part) == (print_jdx - jdx*no_part):
                        Q_HC[(print_idx, print_jdx)] = 2*B
                    elif print_idx - idx*no_part == print_jdx:
                        Q_HC[(print_idx, print_jdx)] = 2*B

    # WRITING EACH NODE * NO_PART
    for node_dx in range(no_vert):
        for idx in range(node_dx*no_part,node_dx*no_part + no_part):
            # print ' '
            #print idx, jdx
            for part, jdx in enumerate(range(node_dx*no_part,node_dx*no_part + no_part)):
                if idx == jdx and part == (no_part - 1):
                    #Q_HC[(idx,jdx)] = -A - (2*no_vert*LoadVector[0]*B/no_part) + B
                    Q_HC[(idx,jdx)] = -A - (2*no_vert*LoadVector[0]*B) + B
                    # print Q_HC[(idx,jdx)],
                elif idx == jdx and part != (no_part - 1):
                    #Q_HC[(idx,jdx)] = -A - (2*no_vert*LoadVector[part]*B/no_part) + B
                    Q_HC[(idx,jdx)] = -A - (2*no_vert*LoadVector[part]*B) + B
                else:
                    Q_HC[(idx,jdx)] = 2*A
                    #print 'else', Q_HC[(idx,jdx)]
                    # print Q_HC[(idx,jdx)],

    # A LIL FIXER FOR DEGENERACY
    for node_dx in range(no_vert):
        for idx in range(node_dx*no_part,node_dx*no_part + 1):
            # print ' '
            for jdx in range(node_dx*no_part,node_dx*no_part + 1):
                if idx == jdx:
                    # Q_HC[(idx,jdx)] = -A - ((len(edgelist)/no_part)-1)*B
                    Q_HC[(idx,jdx)] += 2*B
                    # print Q_HC[(idx,jdx)],
                    #print 'D', idx, jdx, Q_HC[(idx,jdx)]

    # print 'Printing Q HARDCODE:'
    # for idx in range(no_part*no_vert):
    #     print ' '
    #     for jdx in range(no_part*no_vert):
    #         print Q_HC[(idx,jdx)],

    # Setting up Q accounting for degeneracy
    #print '\nQ:'
    Qdeg = dict()
    for idx in range(no_part, no_vert*no_part):
        # print ''
        for jdx in range(no_part, no_vert*no_part):
            # print '%.f' % (Q[(idx,jdx)]),
            Qdeg[(idx - no_part,jdx - no_part)] = (Q_HC[(idx,jdx)])

    (h, J, offset) = dwave_sapi.qubo_to_ising(Qdeg)

    num_reads = 1000

    if solver_type == 'dwave':
        qubits = dwave_sapi.get_hardware_adjacency(dwave_solver)
        q_size = dwave_solver.properties["num_qubits"]

        # Using the D-Wave API heuristic to find an embedding on the Chimera graph for the problem
        embeddings = dwave_sapi.find_embedding(J, len(h), qubits, q_size)

        # Sending problem to the solver
        embedding_solver = dwave_sapi.EmbeddingSolver(dwave_solver, embeddings)
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        sols = answers['solutions']

        t0 = time.time()
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        t1 = time.time()
        answers['timing'] = t1 - t0

        answers['num_reads'] = num_reads

        # adding offset from Q to h, J conversion back to energies
        answers['energies'] = [energy + offset for energy in answers['energies']]
        answers['hvalues'] = h
        answers['jvalues'] = J
        answers['offset'] = offset

    elif solver_type == 'isakov':

        answers = isakovlib.solve_Ising(h, J, offset = offset)

        # adding degeneracy
        negative = '-'
        deg_plusminus = '+' + negative*(no_part-1)
        answers['plusminus'] = [deg_plusminus + sol for sol in answers['plusminus']]

    else:
        print "Moses Partitioner Error: Solver type not recognised. Try dwave or isakov"

    # adding degeneracy back in
    deg = [1] + [-1]*(no_part-1)
    answers['solutions'] = [deg + sol for sol in answers['solutions']]

    # Creating Pymetis style node output based on optimal solution
    opt_sol = answers['solutions'][0]
    print "Optimal Solution", opt_sol

    final_solution = []
    incorrect = 0
    for nodes in range(no_vert):
        colour_sols = []
        for idx in range(nodes*no_part,(nodes*no_part)+no_part):
            colour_sols.append(opt_sol[idx])

        if sum(colour_sols) == (no_part*-1)+2:
            for idx in range(nodes*no_part,(nodes*no_part)+no_part):
                if opt_sol[idx] == 1:
                    final_solution.append(idx - nodes*no_part + 1)
                    # print final_solution
        else:
            final_solution.append(0)
            incorrect += 1

    print "Final_Solution", final_solution

    answers['nodelist'] = final_solution
    answers['adjlist'] = adjlist
    answers['num_parts'] = no_part
    answers['load_vector'] = LoadVector
    answers['solver'] = solver_type
    answers['Q'] = Q_HC
    answers['Qdeg'] = Qdeg
    answers['A'] = A
    answers['B'] = B
    answers['C'] = C

    return answers


def partition(adjlist, solver_type = None, dwave_solver = None, isakov_address = None):
    """

    Partition function based on Lucas Hamiltonian. This one can only divide
    into 2 parts, but requires (Num_Elements - 1) fully connected qubits.

    :param adjlist:
    :param solver_type:
    :param dwave_solver:
    :param isakov_address:
    :return:
    """

    edgelist = adjlist_to_edgelist(adjlist)
    graph = nx.Graph()
    graph.add_edges_from(edgelist)
    nodes = range(len(adjlist))
    graph.add_nodes_from(nodes)

    max_degree = max(graph.degree().values())
    num_vertices = len(adjlist)
    test_size = num_vertices

    # Calculating A, assuming B = 1, ensuring a minimum of 1
    A = max(min(2*max_degree, num_vertices) / 8.0, 1)

    # Number of Qubits Required
    num_qubits = test_size - 1 # Minus 1 for degeneracy

    # Bias for each qubit
    h = [2*A] * num_qubits
    num_reads = 10000 # number of reads from D-Wave
    annealing_time  = 20

    J = dict()
    # Set up the default couplings to have a value of 2
    for i in range(num_qubits):
        for j in range(num_qubits):
            J[(i, j)] = 2*A

    print "Edgelist", edgelist
    print 'Length h', len(h)
    print 'h', h

    # Set up the existing connections to have a value of 1.5
    for i,j in edgelist:
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

    # expected energy
    X = num_vertices + len(edgelist)*0.5

    answers = dict()
    errors = []
    if solver_type == 'dwave':


        print "\nSolving Partition with D-Wave Local Solver"

        qubits = dwave_sapi.get_hardware_adjacency(dwave_solver)
        q_size = dwave_solver.properties["num_qubits"]

        # Using the D-Wave API heuristic to find an embedding on the Chimera graph for the problem
        embeddings = dwave_sapi.find_embedding(J, len(h), qubits, q_size)
        print "Embeddings", embeddings
        # Sending problem to the solver
        embedding_solver = dwave_sapi.EmbeddingSolver(dwave_solver, embeddings)
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)


        MaxChain = max(len(L) for L in embeddings)
        print 'The Maximum Chain for the Embedding is %s long.' % (MaxChain)
        t0 = time.time()
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        t1 = time.time()
        answers['timing'] = t1 - t0
        answers['maxchain'] = MaxChain
        sols = answers['solutions']

    elif solver_type == 'isakov':
        answers = isakovlib.solve_Ising(h, J)

        # adding degeneracy
        negative = '-'
        deg_plusminus = '+'
        answers['plusminus'] = [deg_plusminus + sol for sol in answers['plusminus']]

    else:
        print "PARTITION ERROR: Solver type not recognised. Try DWave or Isakov"

    # adding degeneracy back in
    deg = [1]
    answers['solutions'] = [deg + sol for sol in answers['solutions']]

    # optimal solution
    opt_sol = answers['solutions'][0]
    answers['opt_sol'] = opt_sol

    # finding length of edge boundary
    answers['edge_length'] = length_edgeboundary(adjlist, opt_sol)

    # creating node list
    answers['node_list'] = split_nodelist(opt_sol, adjlist)

    answers['num_qubits'] = num_qubits
    answers['num_reads'] = num_reads
    answers['expected_energy'] = X

    answers['h'] = h
    answers['J'] = J

    answers['adjlist'] = adjlist

    return answers

# needs serious serious serious improvement
def recursive_bisection(n, adjlist, solver_type = None, dwave_solver = None, isakov_address = None):
    """
    recursive_bisection(n, adjlist, solver_type = None, dwave_solver = None, isakov_address = None):
    A recursive bisection algorithm based on function partition for 2^N
    partitions.

    Some may claim 'recursive' also pertains to the repetitive cursing that
    takes places in understanding the code.

    It works just trust me!

    :param n: 2^n partitions
    :param adjlist: adjacency list of original graph
    :param solver_type: 'isakov' or 'dwave'
    :param dwave_solver:
    :param isakov_address:
    :return:
    """
    global output
    global nodes_orig1

    print 'Partitioning', len(adjlist), 'elements...'

    if pow(2,n) > len(adjlist):
        print 'Recursive Bisection Error: n is too large. The number of' \
              ' partitions required is greater than the number of nodes.'
        print '\t%s partitions is greater than %s nodes' % (pow(2, n), len(adjlist))
        print '\tSee Nuclear Fission, or "Splitting the Atom".'

        exit()
    elif n == 1:
        results = partition(adjlist, solver_type = solver_type, dwave_solver = dwave_solver, isakov_address = isakov_address)
        opt_sol1 = results['opt_sol']

        print opt_sol1

        [adjlist_A_orig, adjlist_B_orig] = split_adjlist(opt_sol1,adjlist)

        try:
            nodes_orig1
        except NameError:
            [nodes_A, nodes_B] = split_nodelist(opt_sol1,adjlist)
            adjlist_A = adjlist_A_orig
            adjlist_B = adjlist_B_orig
        else:
            [nodes_A, nodes_B] = split_nodelist_fromnodes(opt_sol1,nodes_orig1)
            adjlist_A = enlarge_adjlist(nodes_orig1,adjlist_A_orig)
            adjlist_B = enlarge_adjlist(nodes_orig1,adjlist_B_orig)

        print "Number of Qubits:", results['num_qubits']

        try:
            output
        except NameError:

            output = dict()
            output['optimal'] = opt_sol1

            output['nodes'] = []
            output['nodes'].append(nodes_A)
            output['nodes'].append(nodes_B)

            output['edge_length'] = []
            output['edge_length'].append(results['edge_length'])

            output['adjacency'] = []
            output['adjacency'] = adjlist_A
            output['adjacency'] = adjlist_B

            output['num_qubits'] = []
            output['num_qubits'].append(results['num_qubits'])

        else:
            #output['optimal'].append([opt_sol1])
            output['nodes'].append(nodes_A)
            output['nodes'].append(nodes_B)
            print "EDGE LENGTH", output['edge_length']
            print "EDGE LENGTH TYPE", type(output['edge_length'])
            output['edge_length'].append(results['edge_length'])
            output['adjacency'].append(adjlist_A)
            output['adjacency'].append(adjlist_B)
            output['num_qubits'].append(results['num_qubits'])

    elif n > 1:
        print "Entered n>1 with n =", n

        try:
            nodes_orig1
        except NameError:
            print "Recursive Bisection: Round One"
        else:
            adjlist_en = enlarge_adjlist(nodes_orig1,adjlist)

        results = partition(adjlist, solver_type = solver_type,
                            dwave_solver = dwave_solver,
                            isakov_address = isakov_address)

        opt_sol = results['opt_sol']
        print "Number of Qubits", results['num_qubits']

        try:
            nodes_orig1
        except NameError:
            [part_A_adj, part_B_adj] = split_adjlist(opt_sol, adjlist)
            [part_A_nodes, part_B_nodes] = split_nodelist(opt_sol, adjlist)
        else:
            [part_A_adj, part_B_adj] = split_adjlist_fromnodes(opt_sol,adjlist_en,nodes_orig1)
            [part_A_nodes, part_B_nodes] = split_nodelist_fromnodes(opt_sol, nodes_orig1)

        reduce_adjlist(part_A_nodes, part_A_adj)
        reduce_adjlist(part_B_nodes, part_B_adj)

        nodes_orig1 = part_A_nodes

        recursive_bisection(n-1, part_A_adj, solver_type = solver_type, dwave_solver = dwave_solver, isakov_address = isakov_address)

        nodes_orig1 = part_B_nodes

        recursive_bisection(n-1, part_B_adj, solver_type = solver_type, dwave_solver = dwave_solver, isakov_address = isakov_address)

    else:
        print "Recursive Bisection Error: I'm sorry, wrong input, ya goose."
        exit()

    output['adjlist'] = adjlist

    return output


##############################################################################
# FUNCTIONS FOR NESTED DISSECTION FOR FILL-REDUCING ON SPARSE MATRICS
##############################################################################


def find_edgeseparators(node_list,edgelist):

    """

    edgeseparator = find_edgeseparators(node_list,edgelist)

    :param node_list: list of lists, with each list containing
    the nodes in each partition
    :param edgelist: all the edges of the graph
    :return: edge_separators: the separators between the nodes in a list
    :return: possible_node_separators: a list of all nodes on the edge separators

    """

    # removing duplicate edges
    for edge in edgelist:
        if (edge[1],edge[0]) in edgelist:
            edgelist.remove((edge[1],edge[0]))

    n = 1

    # If you are reading this, forgive the original coder. This is horrible.
    # Please update append to the following comment.
    # total_time_wasted_trying_to_understand = 0 hours 10 minutes
    edge_separators = []
    possible_node_separators = []
    for part_1 in range(2**n):
        for part_2 in range(2**n):
            if part_1 != part_2:
                for node_A in node_list[part_1]:
                    for node_B in node_list[part_2]:
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

    return edge_separators, possible_node_separators


def vertexcover(adjlist, Node_A = None, Node_B = None):
    '''

    Vertex cover generates the h, J values for vertex cover Hamiltonian in
    Lucas.
    This was converted from a QUBO to an Ising spin problem by hand

    :param input: input can be an edge list of adjlist
    :return:
    '''

    edgelist = adjlist_to_edgelist(adjlist)

    N = len(adjlist)
    num_edges = len(edgelist)

    # coefficients
    # A = (N)*0.5
    A = 3
    B = 1.5
    C = 0.25*1/len(edgelist)

    # determining optimiser term (B term)
    h = []
    for idx in range(N):
        h.append(0.5*B)

    # determining penalty term (A term)
    J = dict()
    for idx in range(N):
        for jdx in range(N):
            J[(idx,jdx)] = 0


    for edge in edgelist:
        h[edge[0]] += -0.25*A
        h[edge[1]] += -0.25*A
        J[(edge[0],edge[1])] = 0.25*A

    # fixer for problem of many subgraphs
    for idx in range(N):
        for jdx in range(N):
            if idx != jdx:
                if (idx, jdx) in J.keys():
                    J[(idx, jdx)] = J[(idx, jdx)] + 0.25
                else:
                    J[(idx, jdx)] = 0.25

    for idx in range(N):
        h[idx] = h[idx] + 0.25


    # Extra term to pick solution to ensure equal nodes from list A and List B
    if Node_A != None and Node_B != None:
        for nidx in Node_A:
            for njdx in Node_A:
                if nidx != njdx:
                    J[(nidx, njdx)] = J[(nidx, njdx)] + 2*C

        for nidx in Node_B:
            for njdx in Node_B:
                if nidx != njdx:
                    J[(nidx, njdx)] = J[(nidx, njdx)] + 2*C

        for nidx in Node_A:
            for njdx in Node_B:
                if nidx != njdx:
                    J[(nidx, njdx)] = J[(nidx, njdx)] - 2*C

    # offset energy
    offset = 0.5*B*N + (num_edges*0.25)

    return h, J, offset


def original_vertexcover(adjlist):
    '''

    Vertex cover generates the h, J values for vertex cover Hamiltonian in
    Lucas.
    This was converted from a QUBO to an Ising spin problem by hand

    :param input: input can be an edge list of adjlist
    :return:
    '''

    edgelist = adjlist_to_edgelist(adjlist)

    N = len(adjlist)
    num_edges = len(edgelist)

    # coefficients
    A = 2
    B = 1

    # determining optimiser term (B term)
    h = []
    for idx in range(N):
        h.append(0.5*B)

    # determining penalty term (A term)
    J = dict()
    for edge in edgelist:
        h[edge[0]] += -0.25*A
        h[edge[1]] += -0.25*A
        J[(edge[0],edge[1])] = 0.25*A

    # offset energy
    offset = 0.5*B*N + (num_edges*0.25)

    return h, J, offset


def nesteddissection(matrix, solver_type = 'isakov', isakov_address = None):
    '''

    Nested Dissection computes an ordering for fill reduction. This is based on
    graph partitioning Hamiltonian and an amended vertex cover Hamiltonian
    from Lucas.

    Additionally prepares a benchmark against nxmetis Nested Dissection

    :param matrix: takes a numpy array
    :param solver_type: choose solver type
    :param isakov_address: choose isakov address
    :return: answers dictionary

    '''

    print 'Is matrix positive definite?',
    if is_pos_def(matrix) is true:
        print 'Yes'
    else:
        print 'Nested Dissection Error: matrix not positive definite.'
        exit()

    G = nx.from_numpy_matrix(matrix)
    adjlist = G.adjacency_list()

    adjlist = remove_diagonals(adjlist)

    ######################################################################
    # Step 0: Benchmark against import nxmetis Nested Dissection
    ######################################################################

    AMD_ordering = nxmetis.node_nested_dissection(G)

    ######################################################################
    # Step 1: Identify edge separators with Graph Partitioning
    ######################################################################

    #answers = partition(adjlist, solver_type='dwave', dwave_solver = solver)
    answers = partition(adjlist, solver_type = 'isakov', isakov_address = isakov_address)

    node_list = answers['node_list']
    opt_list = qmosesnodelist_to_pymetisoutput(node_list)

    edgelist = adjlist_to_edgelist(adjlist)
    edge_separators, possible_node_separators = find_edgeseparators(node_list,edgelist)

    ######################################################################
    # Step 2: Identify node separators with Vertex Cover
    ######################################################################

    # Finding node separator from edge separator using networkx function
    G_sub = nx.Graph()
    G_sub.add_edges_from(edge_separators)
    print 'Is the graph of edge separators bipartite?', bipartite.is_bipartite(G_sub)
    matching = nx.bipartite.maximum_matching(G_sub)
    vertex_cover_nx = list(nx.bipartite.to_vertex_cover(G_sub, matching))
    Snx = vertex_cover_nx

    # finding sides of bipartite graph
    X, Y = bipartite.sets(G_sub)
    NodesA = list(X)
    NodesB = list(Y)

    # Finding node separator using Vertex Cover Hamiltonian
    edge_sep_adjlist = edgelist_to_adjlist(edge_separators)
    h, J, offset = vertexcover(edge_sep_adjlist, Node_A = NodesA, Node_B = NodesB)

    answer = isakovlib.solve_Ising(h,J,offset = offset, solver_address = isakov_address)

    # adding offset from Q to h, J conversion back to energies
    answer['energies'] = [energy + offset for energy in answer['energies']]

    vertex_cover_qc = sol_to_vertexcover(answer['solutions'][0])
    S = vertex_cover_qc

    ####################################################################
    # Step 3: Now to rearrange the bloody thing
    ####################################################################

    # REMOVE NODE SEPARATORS S TO CREATE Q and P
    Qnodes = [x for x in node_list[0] if x not in S]
    Pnodes = [x for x in node_list[1] if x not in S]

    Q_new = [idx for idx in range(len(Qnodes))]
    P_new = [idx+len(Qnodes) for idx in range(len(Pnodes))]

    # Ensures S_new is in correct order
    S_new = [(len(adjlist) - len(S) + x) for x in range(len(S))]

    RowOrdering = Qnodes + Pnodes + S
    NodeList_new = Q_new + P_new + S_new

    def ordering_to_perm_matrix(p):
        '''

        Based on a list of row orderings, produces a numpy array
        as the permutation matrix

        :param p: List of row orderings
        :return: numpy array permutation matrix
        '''

        perm_list = []
        size = len(p)
        for row in p:
            perm_row = [0]*size
            perm_row[row] = 1
            perm_list.append(perm_row)

        perm_matrix = np.array(perm_list)
        return perm_matrix

    def reorder_matrix(perm_matrix, M):
        '''

        Reorders matrix based on permutation matrix
        P^T * A * P

        :param perm_matrix:
        :param M:
        :return:
        '''
        res = np.dot(np.transpose(perm_matrix), M)
        reordered_matrix = np.dot(res,perm_matrix)
        return reordered_matrix

    ####################################################################
    # Comparing Results
    ####################################################################

    dw_perm_matrix = ordering_to_perm_matrix(RowOrdering)
    dw_matrix = reorder_matrix(dw_perm_matrix, matrix)
    print 'Is Nested Dissection reordering positive definite?', is_pos_def(dw_matrix)

    amd_perm_matrix = ordering_to_perm_matrix(AMD_ordering)
    amd_matrix = reorder_matrix(amd_perm_matrix, matrix)
    print 'Is Approx. Min. Degree reordering positive definite?', is_pos_def(amd_matrix)

    Lorig = sl.cholesky(matrix) # defaults on upper triangular
    nz_count_orig = np.count_nonzero(Lorig)

    Ldw = sl.cholesky(dw_matrix)
    nz_count_dw = np.count_nonzero(Ldw)

    Lamd = sl.cholesky(amd_matrix)
    nz_count_amd = np.count_nonzero(Lamd)

    print 'No. Non-Zero Elements (Lower the better)'
    print 'No Reordering:', nz_count_orig
    print 'Dw Reordering:', nz_count_dw
    print 'AMD Reorderin:', nz_count_amd

    ####################################################################
    # Printing/saving results results
    ####################################################################

    print 'AMD Ordering:', AMD_ordering
    print 'DW  Ordering:', RowOrdering

    answer['Q'] = Qnodes
    answer['P'] = Pnodes
    answer['S'] = S
    answer['edge_separators'] = edge_separators
    answer['ordering_amd'] = AMD_ordering
    answer['ordering_dw'] = RowOrdering
    answer['adjlist'] = adjlist
    answer['original_matrix'] = matrix
    answer['matrix_amd'] = amd_matrix
    answer['matrix_dw'] = dw_matrix
    answer['nz_count_orig'] = nz_count_orig
    answer['nz_count_dw'] = nz_count_dw
    answer['nz_count_amd'] = nz_count_amd

    return answer


def remove_diagonals(adjlist):
    '''

    Removes nodes that link to themselves.

    This only confuses partition function.

    :param adjlist:
    :return:
    '''

    for idx, node in enumerate(adjlist):
        for adj in node:
            if adj == idx:
                node.remove(adj)

    return adjlist


def sol_to_vertexcover(opt_sol):
    '''

    Simply converts solution from vertex cover into a list of node of the
    vertex cover.

    :param opt_sol:
    :return:
    '''
    vertex_cover = []
    for idx, node in enumerate(opt_sol):
        if node == 1:
            vertex_cover.append(idx)

    return vertex_cover

# needs improvement
def adjmatrix_to_cvxformat(node_list,adjmatrix):
    '''

    Takes a <class 'scipy.sparse.csr.csr_matrix'> and returns the format
    required to make a spmatrix in the cvxopt module.

    :param node_list:
    :param adjmatrix:
    :return:
    '''
    nonzero = []
    idx_loc = []
    jdx_loc = []
    for idx in range(len(node_list)):
        for jdx in range(len(node_list)):
            if adjmatrix[(idx,jdx)] != 0:
                nonzero.append(adjmatrix[(idx,jdx)])
                idx_loc.append(idx)
                jdx_loc.append(jdx)
    return(nonzero, idx_loc, jdx_loc)


def sparseBanded(rank, bandwidth=1):
    '''

    Function generates a random sparse symmetric positive definite banded
    matrix.

    Needs some work - advise to use Rogues Wathen
    See: https://github.com/macd/rogues

    :param size:
    :return:
    '''

    max_val = 5

    center = [randint(8, max_val*(bandwidth+1)) for r in xrange(rank)]

    sides = []
    bandwidths = [0]
    for side in range(bandwidth):
        sides.append([randint(1, max_val) for r in xrange(rank)])
        bandwidths.append(side+1)

    data = np.vstack((center, sides))
    diags = np.array(bandwidths)

    S = sp.spdiags(data, diags, rank, rank).toarray()
    S = (S + rank*sp.identity(rank))
    S = (S + S.transpose())/2  # ensures matrix is symmetric

    return S


def sparseSym(rank, density=0.01, format='csr', dtype=None, random_state=None):

    '''

    Function generates a random sparse symmetric positive definite matrix by
    first generating a sparse matrix using scipy.sparse, followed by ensuring
    it's symmetrical by adding its transpose.

    This is followed by ensuring the matrix is diagonally dominant and hence
    that it is symmetric positive definite by adding an identity matrix
    multiplied by the size.

    :param rank: Size of matrix
    :param density:
    :param format:
    :param dtype:
    :param random_state:
    :return: S3: returns numpy array
    '''

    density = density / (2.0 - 1.0/rank)
    S = sp.rand(rank, rank, density=density, format=format, dtype=dtype, random_state=random_state)
    S = (S + S.transpose())/2  # ensures matrix is symmetric
    S = (S + rank*sp.identity(rank))
    S = S.A # converts it to numpy array

    return S


def is_pos_def(x):
    '''

    Checks if matrix is positive definite by checking if all
    eigenvalues are greater than zero

    :param x:
    :return:
    '''
    return np.all(np.linalg.eigvals(x) > 0)


###############################################################################
# UNDER THE HOOD: The following make other functions work
###############################################################################


def length_edgeboundary(adjlist,opt_sol):

    """

    no_edges = length_edgeboundary(adjlist,opt_sol):

    Counts the number of edges cut between two partitions.
    Requires networkx not very efficient

    :param adjlist:
    :param opt_sol:
    :return: no_edges: integer number of edges
    """

    # Converting Adjacency List to a Graph
    global_edgelist = adjlist_to_edgelist(adjlist)
    graph = nx.Graph()
    graph.add_edges_from(global_edgelist)

    [part_A_nodes, part_B_nodes] = split_nodelist(opt_sol,adjlist)

    graph_part_A = graph.subgraph(part_A_nodes)
    graph_part_B = graph.subgraph(part_B_nodes)
    edgelist_A = nx.edges(graph_part_A)
    edgelist_B = nx.edges(graph_part_B)

    no_edges = (len(global_edgelist)/2) - (len(edgelist_A) + len(edgelist_B))

    return no_edges


def adjlist_to_edgelist(adjlist):
    """

    edgelist = adjlist_to_edgelist(adjlist: Pretty self explanatory. Takes
    an adjacency list and converts it to an edgelist, how useful!

    The functions inverse is 'edgelist_to_adjlist', naturally!

    :param adjlist: adjacency list
    :return: edgelist: list of edges

    """

    edgelist = []
    for j, y in enumerate(adjlist):
        for i in y:
            edgelist.append((j, i))

    return edgelist


def edgelist_to_adjlist(edgelist):
    """

    edgelist = adjlist_to_edgelist(adjlist: Pretty self explanatory. Takes
    an adjacency list and converts it to an edgelist, how useful!

    The functions inverse is 'edgelist_to_adjlist', naturally!

    :param: edgelist: list of edges
    :return: adjlist: adjacency list

    """

    max_iter = 0
    for edge in edgelist:
        if max(edge) > max_iter:
            max_iter = max(edge)

    adjlist = [[] for i in range(max_iter+1)]
    for edge in edgelist:
        adjlist[edge[0]].append(edge[1])

    return adjlist


def split_adjlist(opt_sol,adjlist):

    """

    split_adjlist(opt_sol,adjlist): Returns two adjacency lists based on the
    output from the QC.

    :param opt_sol: a list of [-1, 1, 1, -1... etc
    :param adjlist: adjacency list
    :return: two adjacency lists

    """

    # Converting Adjacency List to a Graph
    edgelist = adjlist_to_edgelist(adjlist)
    graph = nx.Graph()
    graph.add_edges_from(edgelist)

    # Splitting Nodes from Optimal Solution
    part_A_nodes = []
    part_B_nodes = []
    for idx, s in enumerate(opt_sol):
        if s == 1:
            part_A_nodes.append(idx)
        elif s == -1 or s == 0:
            part_B_nodes.append(idx)
        else:
            print 'ERROR in split_adjlist: geez Louise - You got something ' \
                  'other than -1 or ' \
                  '1 in your sol.'

    # Converting into Graph files from subgraphs
    graphA = graph.subgraph(part_A_nodes)
    graphB = graph.subgraph(part_B_nodes)

    # Getting adjlist from each graph
    adjlist_A = graphA.adjacency_list()
    adjlist_B = graphB.adjacency_list()

    return [adjlist_A, adjlist_B]


def split_nodelist(opt_sol, adjlist):

    """

    split_nodelist(opt_sol,adjlist) spits out a list of lists with first list
    the nodes in partition A and the second is the ndoes in Paritition B

    :param opt_sol:
    :param adjlist:
    :return:

    """

    # Splitting Nodes from Optimal Solution
    part_A_nodes = []
    part_B_nodes = []
    for idx, s in enumerate(opt_sol):
        if s == 1:
            part_A_nodes.append(idx)
        elif s == -1 or s == 0:
            part_B_nodes.append(idx)
        else:
            print 'ERROR in split_nodelist: geez Louise - You got something ' \
                  'other than -1 or ' \
                  '1 in your sol.'

    return [part_A_nodes, part_B_nodes]

# needs improvements
def split_adjlist_fromnodes(opt_sol,adjlist,nodes_orig):

    """

    Splits the partition into two adjacency lists based on a list of original
    nodes. This was created due to the adjacency list reduction function, the
    change the element the in the adjacency list to the least.

    :param opt_sol:
    :param adjlist:
    :param nodes_orig:
    :return:

    """

    # Converting Adjacency List to a Graph
    edgelist = adjlist_to_edgelist(adjlist)
    graph = nx.Graph()
    graph.add_edges_from(edgelist)

    # Splitting Nodes from Optimal Solution
    part_A_nodes = []
    part_B_nodes = []
    for idx, s in enumerate(opt_sol):
        if s == 1:
            part_A_nodes.append(nodes_orig[idx])
        elif s == -1 or s == 0:
            part_B_nodes.append(nodes_orig[idx])
        else:
            print 'ERROR in split_adjlist: geez Louise - You got something ' \
                  'other than -1 or ' \
                  '1 in your sol.'

    # Converting into Graph files from subgraphs
    graphA = graph.subgraph(part_A_nodes)
    graphB = graph.subgraph(part_B_nodes)

    # Getting adjlist from each graph
    adjlist_A = graphA.adjacency_list()
    adjlist_B = graphB.adjacency_list()

    return [adjlist_A, adjlist_B]


def split_nodelist_fromnodes(opt_sol,nodes_orig):
    """

    Splits the partition into two nodes lists based on a list of original
    nodes. This was created due to the adjacency list reduction function, the
    change the element the in the adjacency list to the least.

    :param opt_sol:
    :param adjlist:
    :param nodes_orig:
    :return:

    """

    part_AA_nodes = []
    part_BB_nodes = []
    for idx, s in enumerate(opt_sol):
        if s == 1:
            part_AA_nodes.append(nodes_orig[idx])
        elif s == -1 or s == 0:
            part_BB_nodes.append(nodes_orig[idx])
        else:
            print 'Error: Function partition contains strange output.'

    return [part_AA_nodes, part_BB_nodes]

# needs improvement
def reduce_adjlist(nodes,adjlist):

    """

    reduce_adjlist(nodes,adjlist): from an adjacency list of nodes with
    strange labels return a reduced adjacency list containing nodes from
    0 to number of nodes.

    Tread with caution.

    Inverse function is enlarge_adjlist.

    i.e. input: [[5,13],[0,14],[0],[5]]
         output: [[1,2],[0,3],[0],[1]]

    :param nodes:
    :param adjlist:
    :return:
    """

    j = 0
    i = 0
    k = 0
    adjlistReturn = adjlist
    for p in nodes:
        k = 0
        for x in adjlist:
            for n in x:
                if p == n:
                    adjlistReturn[k][j] = i
                    j += 1
                else:
                    j += 1
            j = 0
            k += 1
        i += 1
    return adjlistReturn

# needs improvement
def enlarge_adjlist(nodes,adjlist):

    """

    enlarge_adjlist(nodes,adjlist): from an adjacency list of nodes with
    strange labels return a reduced adjacency list containing nodes from
    0 to number of nodes.

    Tread with caution.

    Inverse function is reduce_adjlist.

    i.e. output: [[1,2],[0,3],[0],[1]]
         input: [[5,13],[0,14],[0],[5]]


    :param nodes:
    :param adjlist:
    :return:
    """

    j = 0
    i = 0
    k = 0
    # Because Python only creates a reference between variables
    # Hence used copy module
    adjlistReturn = copy.deepcopy(adjlist)
    for p in nodes:
        # print 'All the', k, 's should be', p
        k = 0
        for x in adjlist:
            # print '\t', x
            for n in x:
                # print '\t', n
                if n == i:
                    # print 'ADDED:', n, "== Position", i
                    # print 'Position', i, 'is a', p
                    # print '=>', adjlist[k][j], "==", p
                    adjlistReturn[k][j] = p
                    j += 1
                else:
                    j += 1
            j = 0
            k += 1
        i += 1
        # print 'Adjlist after that round: \n', adjlistReturn
        # print 'Original:\n', adjlist
        # print "\n NEW ROUND \n \n"
    return adjlistReturn


def qmosesnodelist_to_pymetisoutput(nodes_list):

    """
    qmosesnodelist_to_pymetisoutput(nodes_list):
    Converts a list of list of nodes in each partition into a similar style
    output to what PyMETIS does.

    e.g. [[0,3,4],[1,2]] to [1,2,2,1,1]

    :param nodes_list:
    :return:
    """

    max_value = 0
    for part in nodes_list:
        max_test = max(part)
        if max_test > max_value:
            max_value = max_test

    output = [0]*(max_value + 1)
    for part_num, part in enumerate(nodes_list):
        for element in part:
            output[element] = part_num + 1

    return output

# needs improvement
def meshpytrielements_to_adjlist(meshpy_elements):

    """

    meshpytrielements_to_adjlist(meshpy_elements):
    Converts MeshPy Triangle 2D or Tetrahedral 3D elements list or array to
    an adjacency file.

    :param meshpy_elements:
    :return:

    """

    if type(meshpy_elements) != list:
        meshpy_elements = meshpy_elements.tolist()
    elif type(meshpy_elements) == list:
        pass
    else:
        print 'Weird Input Try again.'
        exit()

    # Dimensions: checking if triangular, or tetrahedral
    dimensions = len(meshpy_elements[0])-1

    # Creates an empty adjacency list
    adjlist = [[] for x in xrange(len(meshpy_elements))]
    globallistno = 0
    for list1 in meshpy_elements:
        count = [[] for x in xrange(len(meshpy_elements))]
        #print 'Entering Global List No.:', globallistno
        for element1 in list1:
            sublistno = 0
            #print 'Element One:', element1
            for list_onwards in meshpy_elements:
                #print '\tExploring Sub-List No.:', sublistno
                for element2 in list_onwards:
                    #print '\t\tElement Two:', element2
                    if element2 == element1:
                        #print '\t\t\tYou BLOODY RIPPER: at', globallistno, sublistno
                        if globallistno == sublistno:
                            continue
                        elif globallistno != sublistno:
                            if count[sublistno] >= [1]:
                                count[sublistno][0] += 1
                            else:
                                count[sublistno].append(1)
                            #print '\t\t\t Count at', sublistno, 'is now:', count[sublistno]
                            if count[sublistno][0] >= dimensions:
                                adjlist[globallistno].append(sublistno)
                                #print '\t\t\t Added count'
                sublistno += 1
        globallistno += 1

    return adjlist


def normalise_weightings(h, J, offset):
    '''

    normalise_weightings finds the maximum coupling weight and then normalises
    all the weightings to that value.

    e.g. h = [3, 1.5, 1, 0.5, 2], h_new = [1, 0.5, 0.3333, 0.1667, 0.6667]

    :param h: H weightings in a list
    :param J: J coupling strengths in a dictionary
    :param offset: offset energy
    :return: Returns the altered input
    '''

    # finding max in h
    maximum = max(h)

    # finding max in J
    for key, value in J.items():
        if value > maximum:
            maximum = value

    # normalising h values
    h = [x/maximum for x in h]

    # normalising J values
    for key, value in J.items():
        J[key] = value / maximum

    factor = 1 / maximum
    print 'Normalising weightings by factor:', factor

    # not certain this is technically correct
    offset *= factor

    return h, J, offset

##############################################################################
# FUNCTIONS FOR RESULTS ANALYSES
##############################################################################

# needs significant improvement
# it's likely this doesn't even output the correct results
# also need to ensure it works for moses partition whereby partitions are
# not actually allocated
def edge_resultsanalysis(list_of_lists_of_nodes,edgelist, num_parts):

    """
    (edge_results, total_edges) =
    edge_resultsanalysis(list_of_lists_of_nodes,edgelist, num_parts)

    Function for analysing the condition of minimising the edge boundary in a
    mesh/graph partitioning problem.

    Args:
        list of lists of nodes: which nodes are in each partition.
        edgelist: edgelist of undirectional graph.
        num_parts: number of partitions

    Returns:
        edges_between: returns a dict for the num of edges between partitions
        total_edges: returns the total number of edges cut in the partitioning
        edges_each: returns list of

    """

    # removing duplicate edges for algorithm
    for edge in edgelist:
        if (edge[1],edge[0]) in edgelist:
            edgelist.remove((edge[1],edge[0]))

    # determines total_edges and edges_between
    edges_between = dict()   # (0,1) = 28, (0,2) = 20 (1,0) =
    total_edges = 0
    for part_1 in range(num_parts):
        for part_2 in range(num_parts):
            no_edges = 0
            if part_1 != part_2:
                for node_A in list_of_lists_of_nodes[part_1]:
                    for node_B in list_of_lists_of_nodes[part_2]:
                        if (node_A, node_B) in edgelist:
                            no_edges += 1
                        elif (node_B, node_A) in edgelist:
                            no_edges += 1
                edges_between[(part_1,part_2)] = no_edges
                total_edges += no_edges

    # for some reason, after testing its half what the solution says
    total_edges = total_edges*0.5

    edges_each = []
    for results in edges_between.keys():
        sum = 0
        avail_connec = [t for t in edges_between.keys() if t[0] == results[0]]
        for relevant in avail_connec:
            # sums key value for relevant boundary edges
            sum += edges_between[relevant]
            #print 'relevant: %s, idx: %s, results: %s' % (relevant, idx, results)
        #print 'Y', results
        edges_each.append(sum)
    edges_each = edges_each[:num_parts]

    return (edges_between, edges_each, total_edges)


def num_nodes_per_part(list_of_lists_of_nodes,num_parts):

    """
    (result, variance) = num_nodes_per_part(list_of_lists_of_nodes,num_parts)

    Function for analysing the balance of nodes in a partition in a
    mesh/graph partitioning problem.

    Args:
        list of lists of nodes: which nodes are in each partition.
        num_parts: number of partitions

    Returns:
        result_list: a list where each element is the number of nodes in that partition.
        variance: returns the variance of the result. Zero is ideal.

    """

    def var(any_list):
        mean = sum(any_list)/len(any_list)
        result = 0
        for i in any_list:
            result += (mean - i) ** 2
        return result/2

    result_list = []
    for part_1 in range(num_parts):
        result_list.append(len(list_of_lists_of_nodes[part_1]))

    variance = var(result_list)

    return [result_list, variance]


def energy_from_solution(h, J, opt_sol, offset = None):

    """
    total_energy = energy_from_solution(h, J, offset = None, opt_sol):

    Function for taking a solution generated by another algorithm and
    calculating the energy based on the hvalues and Jvalues.

    Args:
        h: h value weightings for Ising spin in a list
        J: J value couplings for Ising spin in a dictionary.
        opt_sol: a list of spins, e.g. [-1, 1, 1, -1...
        offset: offset energy if applicable

    Returns:
        total_energy

    """

    hvalue_energy = sum([a*b for a, b in zip(opt_sol, h)])

    # calculating J value energy
    Jvalue_energy = 0
    for idx in range(len(h)):
        for jdx in range(len(h)):
            #print idx, jdx
            if (idx,jdx) in J:
                Jvalue_energy += J[(idx, jdx)]*opt_sol[idx]*opt_sol[jdx]

    if offset == None:
        total_energy = Jvalue_energy + hvalue_energy
    else:
        total_energy = Jvalue_energy + hvalue_energy + offset

    return total_energy

##############################################################################
# MISC FUNCTIONS - ATTEMPTS AT OTHER THINGS
##############################################################################

def colourgraph(no_col,adjlist):
    '''

    Colour Graph: it's not a matter of colouring between the lines, it's a
    matter of colouring the nodes on the lines (or edges).

    Hardcoded implementation of the graph coloring Hamiltonian found in Lucas

    Does not consider permutation symmetry and likely there will be degeneracy
    in the solution.

    If an energy of -1*no_vert is not acquired, there is likely no solution
    for that number of colours. e.g. 20 nodes, opt sol should yield -20

    :param no_col:
    :param adjlist:
    :return:
    '''

    no_vert = len(adjlist)  # Number of vertices
    no_part = no_col

    no_qubits = (no_vert*no_part)

    edgelist = adjlist_to_edgelist(adjlist)

    # Coefficients
    A = 1
    B = 1

    # TERM 1: The term A ensures every vertex has one colour, taken from Moses
    # Partitioner

    Q1 = dict()
    for idx in range(no_qubits):
        for jdx in range(no_qubits):
            Q1[(idx, jdx)] = 0

    # The diagonals
    for idx in range(no_qubits):
        Q1[(idx, idx)] = -2*A + 1*A # 1*B is for (x1)^2 terms

    # might be possible solution without itertools
    from itertools import combinations
    for vertdx in range(no_vert):
        for Q1_key in list(combinations(range(vertdx*no_part, (vertdx+1)*no_part), 2)):
                Q1[Q1_key] = 2*A

    # print "Printing QA or Q1, whatever ya wanna call it"
    #
    # for idx in range(no_qubits):
    #     print ''
    #     for jdx in range(no_qubits):
    #         print Q1[(idx,jdx)],

    # TERM 2: Energy penalty each time an connects two vertices of the same colour

    Q2 = dict()
    for idx in range(no_qubits):
        for jdx in range(no_qubits):
            Q2[(idx, jdx)] = 0

    for i, j in edgelist:
        if j > i:
            for part in range(no_part):

                idx = i*no_part + part
                jdx = j*no_part + part

                Q2[(idx, jdx)] += 1*B

    # print "\nPrinting QB or Q2, the edges one"

    # for idx in range(no_qubits):
    #     print ''
    #     for jdx in range(no_qubits):
    #         print Q2[(idx,jdx)],

    # THE FINAL Q
    Q_final = dict()
    for idx in range(no_qubits):
        print ''
        for jdx in range(no_qubits):
            Q_final[(idx, jdx)] = Q1[(idx, jdx)] + Q2[(idx, jdx)]
            print Q_final[(idx, jdx)],

    return Q_final
