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
from collections import defaultdict
import copy as copy
import subprocess
import os
from sympy import *

# GRAPH COLOURING FUNCTION
def colourgraph(no_col,adjlist,solver_type = None, dwave_solver = None, isakov_address = None):
    no_vert = len(adjlist)  # Number of vertices

    print 'FUNCTION COLOURGRAPH Beginning...'
    print 'Number of nodes:', no_vert
    print 'Number of colours:', no_col
    print 'Therefore need %s qubits.\n' % (no_vert*no_col)

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

    print 'Final Hamiltonian:', H

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

    if solver_type == 'dwave':
        "Using dwave solver for Graph Colouring problem..."
        ## PRINT ERROR IF SOLVER ISN'T DEFINED
        qubits = dwave_sapi.get_hardware_adjacency(dwave_solver)
        q_size = dwave_solver.properties["num_qubits"]

        # Using the D-Wave API heuristic to find an embedding on the Chimera graph for the problem
        embeddings = dwave_sapi.find_embedding(J, len(h), qubits, q_size)

        # Sending problem to the solver
        embedding_solver = dwave_sapi.EmbeddingSolver(dwave_solver, embeddings)
        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        sols = answers['solutions']

        answers = embedding_solver.solve_ising(h, J, num_reads = num_reads)
        sols = answers['solutions']
        print answers['energies']
        # Adding Degeneracy Back In
        Degenerate = [1] + [-1]*(no_col-1)
        opt_sol = Degenerate + sols[0]

    elif solver_type == 'isakov':
        output = defaultdict(list)
        print "Using Isakov to solve Graph Colouring problem..."
        # PRINT ERROR IF SOLVER ADDRESS ISN'T DEFINED
        # PRINT LATTICE FILE

        if os.path.isfile("ColourGraph.lattice") == 1:
            os.remove("ColourGraph.lattice")

        file = open("ColourGraph.lattice", "w")
        file.write("ColourGraph Function File")

        # Printing h values to Lattice File
        for idx, element in enumerate(h):
            file.write("\n%s %s %s" % (idx, idx, element))

        # Printing J values to Lattice File
        #print 'J values:',
        for idx in range(len(h)):
            #print ' '
            for jdx in range(len(h)):
                if (idx, jdx) in J:
                    file.write("\n%s %s %s" % (idx, jdx, J[(idx, jdx)]))
                # if i >= j:
                #     continue
                # else:
                #     file.write("\n%s %s %s" % (i, j, J[(i,j)]))
                #     #file.write("\n%s %s %s" % (j, i, J_dict[(i,j)]))

        file.close()

        lattice_file = "ColourGraph.lattice"

        args = [isakov_address, "-l", lattice_file, "-s", "100", "-r", "1000"]
        #subprocess.call(args)# Uncomment to run again and print output
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        #popen.wait()
        isakov_output = popen.stdout.read()
        isakov_output = isakov_output.splitlines()

        solution_listoflist = []
        for s in isakov_output:
            solution_listoflist.append(s.split())

        for s in solution_listoflist:
            output['energy'].append(float(s[0]))
            output['occurences'].append(int(float(s[1])))
            output['plusorminus'].append(s[2])

        for element in output['plusorminus']:
            conversion = []
            for spin in element:
                if spin == '+':
                    conversion.append(1)
                elif spin == '-':
                    conversion.append(-1)
                else:
                    print "It's output something other than +, -"
            output['solutions'].append(conversion)

        sols = output['solutions']
        print 'SOLS', sols
        print J

        # Adding Degeneracy Back In
        Degenerate = [1] + [-1]*(no_col-1)
        opt_sol = Degenerate + sols[0]

    else:
        print "COLOUR GRAPH ERROR: Solver type not recognised. Try dwave or isakov"

    print 'Convert to nice format from opt_sol =', opt_sol
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



# GRAPH PARTITION FUNCTION
# Imports a variety of graph/mesh/adjacency files, you give the number of partitions,
# whether to use local or remote solver,
# Partitions into any number of partitions
# sub functions include Lucas Hamiltonian for 2^N partitions
# Your own Hamiltonian for prime numbers etc etc

# Function Splits Adjacency List Into Two Parts Based on the Hamiltonian
# Requires solver information

def partition(adjlist,solver_type = None, dwave_solver = None, isakov_address = None):
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

    output = defaultdict(list)
    output['num_qubits'] = num_qubits

    if solver_type == 'dwave':
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

        print "Optimum Solution:", opt_sol
        print "Type:", type(opt_sol)

        print "HERE TIS!", output

    elif solver_type == 'isakov':
        print "fine! I'll do Isakov"
        # PRINT ERROR IF SOLVER ADDRESS ISN'T DEFINED
        # PRINT LATTICE FILE

        if os.path.isfile("PartitionLattice.lattice") == 1:
            os.remove("PartitionLattice.lattice")

        file = open("PartitionLattice.lattice", "w")
        file.write("Partition Function File")

        # Printing h values to Lattice File
        for idx, element in enumerate(h):
            file.write("\n%s %s %s" % (idx, idx, element))

        # Printing J values to Lattice File
        #print 'J values:',
        for i in range(len(h)):
            #print ' '
            for j in range(len(h)):
                if i >= j:
                    continue
                else:
                    file.write("\n%s %s %s" % (i, j, J[(i,j)]))
                    #file.write("\n%s %s %s" % (j, i, J_dict[(i,j)]))

        file.close()

        lattice_file = "PartitionLattice.lattice"

        args = [isakov_address, "-l", lattice_file, "-s", "200", "-r", "1000"]
        #subprocess.call(args)
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        #popen.wait()
        isakov_output = popen.stdout.read()
        isakov_output = isakov_output.splitlines()

        solution_listoflist = []
        for s in isakov_output:
            solution_listoflist.append(s.split())

        for s in solution_listoflist:
            output['energy'].append(float(s[0]))
            output['occurences'].append(int(float(s[1])))
            output['plusorminus'].append(s[2])

        for element in output['plusorminus']:
            conversion = []
            for spin in element:
                if spin == '+':
                    conversion.append(1)
                elif spin == '-':
                    conversion.append(-1)
                else:
                    print "It's output something other than +, -"
            output['solutions'].append(conversion)


        # YOU'RE NOT ADDING THE EXTRA QUBIT YA TWAT
        sols = output['solutions']
        map(lambda sols: sols.append(1), output['solutions'])
        opt_sol = sols[0]
        output['opt_sol'] = opt_sol

    else:
        print "PARTITION ERROR: Solver type not recognised. Try DWave or Isakov"

    # Finding length of edge boundary
    output['edge length'] = length_edgeboundary(output,opt_sol)

    return output

# A RECURSIVE BISECTION SCHEME FOR 2^N PARTITIONS
# Based on previously defined function, partition.
# INPUT recursive_bisection(n,adjlist,solver)
# N:        - 2^N partitions
# adjlist   - adjacency list of original graph
# solver    - D-Wave Solver
# OUTPUT list of nodes containing partition number

def recursive_bisection(n, adjlist, solver_type = None, dwave_solver = None, isakov_address = None):
    global output
    global nodes_orig1

    print 'Partitioning', len(adjlist), 'elements...'

    if pow(2,n) > len(adjlist):
        print pow(2,n), '>', adjlist
        print 'Recursive Bisection Error: n is too large. The number of' \
              ' partitions required is greater than the number of nodes.' \
              ' \n  See Nuclear Fission, or "Splitting the Atom".'
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
        #print 'Adjlist A', adjlist_A
        #print "Child Nodes A:", nodes_A
        #print "Child Nodes B:", nodes_B

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
        #print "Adjlist:", adjlist
        try:
            nodes_orig1
        except NameError:
            print "Round One"
        else:
            adjlist_en = enlarge_adjlist(nodes_orig1,adjlist)
            #print "Adjlist Orig:", adjlist_en
            #print "Nodes Orig:", nodes_orig1

        #print "GOT TO PARTITION"
        results = partition(adjlist, solver_type = solver_type,
                            dwave_solver = dwave_solver,
                            isakov_address = isakov_address)

        opt_sol = results['opt_sol']
        print "Number of Qubits", results['num_qubits']
        #print "GOT THROUGH PARTITION"

        try:
            nodes_orig1
        except NameError:
            [part_A_adj, part_B_adj] = split_adjlist(opt_sol, adjlist)
            [part_A_nodes, part_B_nodes] = split_nodelist(opt_sol, adjlist)
            #print "GOT TO ZERO"
            #print "A Nodes", part_A_nodes
            #print "B Nodes", part_B_nodes
        else:
            ### IN HERE, NON-REDUCE ADJLIST
            #adjlist = enlarge_adjlist(nodes_orig1,adjlist)
            [part_A_adj, part_B_adj] = split_adjlist_fromnodes(opt_sol,adjlist_en,nodes_orig1)
            [part_A_nodes, part_B_nodes] = split_nodelist_fromnodes(opt_sol, nodes_orig1)
            #print "GOT TO TWO"
            #print "A Nodes", part_A_nodes
            #print "B Nodes", part_B_nodes

        #print "Adjlist A, wo reduction:", part_A_adj

        reduce_adjlist(part_A_nodes, part_A_adj)
        reduce_adjlist(part_B_nodes, part_B_adj)

        nodes_orig1 = part_A_nodes
        #print "Nodes Original", nodes_orig1
        #print "Adjlist A, with reduction::", part_A_adj
        #print "n", n
        recursive_bisection(n-1, part_A_adj, solver_type = solver_type,
                            dwave_solver = dwave_solver,
                            isakov_address = isakov_address)

        nodes_orig1 = part_B_nodes
        #print "Nodes Original", nodes_orig1
        #print "Adjlist B:", part_B_adj
        #print "n", n
        recursive_bisection(n-1, part_B_adj,  solver_type = solver_type,
                            dwave_solver = dwave_solver,
                            isakov_address = isakov_address)

    else:
        print "Recursive Bisection Error: I'm sorry, wrong input, ya goose."
        exit()

    return output

def length_edgeboundary(adjlist,opt_sol):
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
    j = 0
    edgelist = []
    for y in adjlist:
        for i in y:
            edgelist.append((j, i))
        j += 1
    return edgelist

# SPLIT
# Returns two adjacency lists given the partition and one adjlist

def split_adjlist(opt_sol,adjlist):

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

def split_nodelist(opt_sol,adjlist):

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
            print 'ERROR in split_nodelist: geez Louise - You got something ' \
                  'other than -1 or ' \
                  '1 in your sol.'
    return [part_A_nodes, part_B_nodes]

def split_adjlist_fromnodes(opt_sol,adjlist,nodes_orig):
    """Splits the partition into two adjacency files based on a lit of original nodes"""
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


def reduce_adjlist(nodes,adjlist):
    j = 0
    i = 0
    k = 0
    adjlistReturn = adjlist
    for p in nodes:
        k = 0
        for x in adjlist:
            for n in x:
                if p == n:
                    # print 'Adjlist K,J:', adjlist_A[k][j]
                    # print 'Success: p:', p,  'n:', n, 'x:', x, 'i:', i,'k:', k, 'j:', j
                    adjlistReturn[k][j] = i
                    j += 1
                else:
                    # print 'WRONG: p:', p,  'n:', n, 'x:', x, 'i:', i,'k:', k, 'j:', j
                    j += 1
            j = 0
            k += 1
        i += 1
    return adjlistReturn

def enlarge_adjlist(nodes,adjlist):
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

# Converts output from recursive bisection of qmoses to the pymetis output
def qmosesnodelist_to_pymetisoutput(nodes_list):
    my_list = [l[0] for l in nodes_list]
    list1 = range(max(map(max, my_list))+1)
    for idx, j in enumerate(my_list):
        for i in j:
            list1[i] = idx
    return list1

# Converts MeshPy Triangle 2D or Tetrahedral 3D elements list or array to
# adjacency list for partitioning
# Output: adjacency list
def meshpytrielements_to_adjlist(meshpy_elements):
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

# FILL-REDUCING ORDERINGS FUNCTION FOR SPARSE MATRICES
# Using graph partitioning, or graph colouring methods, even clique methods

# MESSAGE SCHEDULING GRAPH COLOURING PROBLEM

# QUANTUM LATTICE GAS for FUN DEMONSTRATION perhaps