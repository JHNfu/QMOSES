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

import networkx as nx
import numpy as nm
import dwave_sapi
from dwave_sapi import local_connection, find_embedding
from collections import defaultdict
import copy as copy
from sympy import *
import isakovlib


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


def vertexcover(adjlist):
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
    A = (N)*0.5
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


##############################################################################
# VARIOUS FUNCTIONS FOR MESH PARTITIONING
##############################################################################


def mosespartitioner(no_part,adjlist, load = None, solver_type = None, dwave_solver = None, isakov_address = None):
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

    if load == None:
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
    graph = nx.Graph()
    graph.add_edges_from(edgelist)
    graph.add_nodes_from(range(len(adjlist)))

    A = 2*no_vert
    B = no_vert*0.1 # KEEP AT ONE
    C = 0.5*no_vert

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

        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        answers['num_reads'] = num_reads

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

    # adding offset from Q to h, J conversion back to energies
    answers['energies'] = [energy + offset for energy in answers['energies']]

    # Creating Pymetis style node output based on optimal solution
    opt_sol = answers['solutions'][0]
    final_solution = []
    incorrect = 0
    for nodes in range(no_vert):
        colour_sols = []
        for idx in range(nodes*no_part,(nodes*no_part)+no_part):
            colour_sols.append(opt_sol[idx])

        if sum(colour_sols) == -1:
            for idx in range(nodes*no_part,(nodes*no_part)+no_part):
                if opt_sol[idx] == 1:
                    final_solution.append(idx - nodes*no_part + 1)
                    # print final_solution
        else:
            final_solution.append(0)
            incorrect += 1

    answers['nodelist'] = final_solution
    answers['adjlist'] = adjlist
    answers['num_parts'] = no_part
    answers['load_vector'] = LoadVector
    answers['solver'] = solver_type
    answers['Q'] = Q_HC
    answers['Qdeg'] = Qdeg

    return answers


def partition(adjlist,solver_type = None, dwave_solver = None, isakov_address = None):

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

    num_reads = 1000 # number of reads from D-Wave
    annealing_time  = 20

    J = dict()
    # Set up the default couplings to have a value of 2
    for i in range(num_qubits):
        for j in range(num_qubits):
            J[(i, j)] = 2*A

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

    answers = dict()
    if solver_type == 'dwave':
        print "\nSolving Partition with D-Wave Local Solver"
        ## PRINT ERROR IF SOLVER ISN'T DEFINED
        qubits = dwave_sapi.get_hardware_adjacency(dwave_solver)
        q_size = dwave_solver.properties["num_qubits"]

        # Using the D-Wave API heuristic to find an embedding on the Chimera graph for the problem
        embeddings = dwave_sapi.find_embedding(J, len(h), qubits, q_size)

        # Sending problem to the solver
        embedding_solver = dwave_sapi.EmbeddingSolver(dwave_solver, embeddings)
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        sols = answers['solutions']

        map(lambda sols: sols.append(1), answers['solutions'])
        opt_sol = sols[0]

        output['opt_sol'] = opt_sol

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
    answers['edge length'] = length_edgeboundary(adjlist, opt_sol)

    # creating node list
    answers['node_list'] = split_nodelist(opt_sol, adjlist)

    answers['num_qubits'] = num_qubits

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
            output = defaultdict(list)
            output['optimal'].append([opt_sol1])
            output['nodes'].append([nodes_A])
            output['nodes'].append([nodes_B])
            output['edge length'].append(results['edge length'])
            output['adjacency'].append(adjlist_A)
            output['adjacency'].append(adjlist_B)
            output['num_qubits'].append(results['num_qubits'])
        else:
            output['optimal'].append([opt_sol1])
            output['nodes'].append([nodes_A])
            output['nodes'].append([nodes_B])
            output['edge length'].append(results['edge length'])
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

    return output


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
        return result

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
            if (idx,jdx) in J and idx > jdx:
                # print J[(i, j)],
                # print 'hlen', len(h),
                # print idx, jdx, opt_sol[idx*jdx], J[(idx, jdx)]
                Jvalue_energy += J[(idx, jdx)]*opt_sol[idx]*opt_sol[jdx]

    if offset == None:
        total_energy = Jvalue_energy + hvalue_energy
    else:
        total_energy = Jvalue_energy + hvalue_energy + offset

    return total_energy


##############################################################################
# MISC FUNCTIONS - ATTEMPTS AT OTHER THINGS
##############################################################################

# needs serious improvement and tests
# based on sympy, not ideal
def colourgraph(no_col,adjlist,solver_type = None, dwave_solver = None, isakov_address = None):
    """

    Colour Graph: it's not a matter of colouring between the lines, it's a
    matter of colouring the nodes on the lines (or edges).

    sympy implementation of the graph coloring Hamiltonian found in Lucas

    :param no_col: number of colours
    :param adjlist:
    :param solver_type:
    :param dwave_solver:
    :param isakov_address:
    :return:
    """
    no_vert = len(adjlist)  # Number of vertices

    Aye = 1
    Bee = 1

    A = Symbol('A')
    B = Symbol('B')

    x = symbols('x(1:%d\,1:%d)' % (((no_vert + 1) , (no_col + 1))))

    H1 = 0
    term = 1
    for v in range(no_vert):
        for i in range(no_col):
            v += 1 # FIX THIS!
            i += 1 # FIX THIS!
            term -= x[((v - 1) * no_col) + i - 1]
            v -= 1 # FIX THIS!
            i -= 1 # FIX THIS!
        term = term**2
        H1 += term
        term = 1
    H1 = A * H1

    edgelist = adjlist_to_edgelist(adjlist)

    # removing duplicate edges
    for edge in edgelist:
        if (edge[1],edge[0]) in edgelist:
            edgelist.remove((edge[1],edge[0]))

    H2 = 0
    for edge in edgelist:
        for colour in range(no_col):
            H2 += x[((edge[0] - 1) * no_col) + colour - 1]*x[((edge[1] - 1) * no_col) + colour - 1]
    H2 = (B)*H2 # Dividing by 2 as edgelist defines edges twice

    #print 'H1:', H1
    #print 'H2:', H2
    H = H1 + H2

    H = H.subs('A', Aye)
    H = H.subs('B', Bee)

    # print 'Hamiltonian (before subs):', H

    # Accounting for Degeneracy

    H = H.subs(x[0], 1)
    for s in range(no_col):
        H = H.subs(x[s],0)

    # print 'Hamiltonian (after subs):', H

    H = H.expand(H)

    #print 'Final Hamiltonian:', H

    Q = dict()
    for idx in range(len(x)):
        for jdx in range(len(x)):
            Q[(idx, jdx)] = 0

    for idx in range(len(x)):
        for jdx in range(len(x)):
            if idx != jdx:
                Q[(idx, jdx)] = round(H.coeff(x[idx]).coeff(x[jdx]),2)
                Q[(idx, jdx)] = round(H.coeff(x[idx]*x[jdx]),2)
            else:
                E = H.coeff(1*x[idx])
                #Q[(idx,jdx)] = round(H.coeff(x[idx]**2 + H.coeff(x[idx])),2)
                Q[(idx,jdx)] = round((H.coeff(x[idx]**2) - E.coeff(-1)),2)

    # Creating a Q accounting for degeneracy
    #print '\nQ:'
    Qdeg = dict()
    for idx in range(no_col, no_vert*no_col):
        #print ''
        for jdx in range(no_col, no_vert*no_col):
            #print '%.f' % (Q[(idx,jdx)]),
            Qdeg[(idx - no_col,jdx - no_col)] = (Q[(idx,jdx)])

    (h, J, offset) = dwave_sapi.qubo_to_ising(Qdeg)

    #print 'H vaules', h
    #print 'J values', J

    num_reads = 1000

    print "Using %s to solve Graph Colouring Problem..." % (solver_type)
    if solver_type == 'dwave':
        ## PRINT ERROR IF SOLVER ISN'T DEFINED
        qubits = dwave_sapi.get_hardware_adjacency(dwave_solver)
        q_size = dwave_solver.properties["num_qubits"]

        # Using the D-Wave API heuristic to find an embedding on the Chimera graph for the problem
        embeddings = dwave_sapi.find_embedding(J, len(h), qubits, q_size)

        # Sending problem to the solver
        embedding_solver = dwave_sapi.EmbeddingSolver(dwave_solver, embeddings)
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)

        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)

    elif solver_type == 'isakov':

        answers = isakovlib.solve_Ising(h, J, offset)

        # adding degeneracy
        negative = '-'
        deg_plusminus = '+' + negative*(no_col-1)
        answers['plusminus'] = [deg_plusminus + sol for sol in answers['plusminus']]

    else:
        print "COLOUR GRAPH ERROR: Solver type not recognised. Try dwave or isakov"

    # adding degeneracy back in
    deg = [1] + [-1]*(no_col-1)
    answers['solutions'] = [deg + sol for sol in answers['solutions']]

    # adding offset from Q to h, J conversion back to energies
    answers['energies'] = [energy + offset for energy in answers['energies']]

    opt_sol = answers['solutions'][0]

    final_solution = []
    incorrect = 0
    for nodes in range(no_vert):
        colour_sols = []
        for idx in range(nodes*no_col,(nodes*no_col)+no_col):
            colour_sols.append(opt_sol[idx])
            # if opt_sol[idx] == 1:
            #     nice_output.append(idx - nodes*no_col + 1)

        if 1 not in colour_sols:
            final_solution.append(0)
            incorrect += 1

        else:
            for idx in range(nodes*no_col,(nodes*no_col)+no_col):
                if opt_sol[idx] == 1:
                    final_solution.append(idx - nodes*no_col + 1)

    if incorrect > 0:
        print 'Some nodes were not worthy of a colour.\n' \
              'ie. There is not a perfect solution to this problem or ' \
              'the solver did not find one.'

    return final_solution





# FILL-REDUCING ORDERINGS FUNCTION FOR SPARSE MATRICES
# Using graph partitioning, or graph colouring methods, even clique methods

# MESSAGE SCHEDULING GRAPH COLOURING PROBLEM

# QUANTUM LATTICE GAS for FUN DEMONSTRATION perhaps