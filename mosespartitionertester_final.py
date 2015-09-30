__author__ = 'John'
##############################################################################
# Moses Partitioner Tester Final
# September 2015
##############################################################################
# This code is used to run mulitple tests of the modes partitioner Hamiltonian
# as defined in the qmoseslib.py. This includes post-processing of energy
# information, balancing of nodes, and minimising edges.
# Data is stored in a results file in a dictionary with each test with its
# own individual key, dumped using the pickle library.
##############################################################################
import cPickle as pickle
from qmoseslib import *
import isakovlib as isakov
import networkx as nx
import time
import pymetis as pm

# Solvers
type_solver = 'isakov'
solver_address = r"C:\Users\John\PycharmProjects\public-isakov\public-isakov\build\visual_studio_2013\bin\an_ss_ge_fi_vdeg_omp.exe"

results = dict()
timestr = time.strftime("%m%d_%H%M%S")
results['starttime'] = timestr
pickle.dump(results, open("Results_%s.results" % (timestr), "wb"))

min_nodes = 25
max_nodes = 31
min_parts = 2
max_parts = 4

# Load Balancing Vector
LoadBalancing = dict()
LoadBalancing[2] = [[0.25, 0.75], [0.75, 0.25], [0.66, 0.33], [0.33, 0.66]]
LoadBalancing[3] = [[0.25, 0.25, 0.5], [0.5, 0.25, 0.25], [0.25, 0.5, 0.25]]
LoadBalancing[4] = [[0.33, 0.33, 0.16, 0.16], [0.16, 0.33, 0.33, 0.16]]

# Initalising Results Dictionary
for num_nodes in range(min_nodes,max_nodes+1):
    for num_parts in range(min_parts,max_parts+1):
        results[(num_nodes,num_parts)] = dict()

##############################################################################
# RUNNING TESTS
##############################################################################
totaltests = (max_nodes-min_nodes+1)*(max_parts-min_parts+1)
test_num = 1
for num_nodes in range(min_nodes,max_nodes+1):
    for num_parts in range(min_parts,max_parts+1):

        req_qubits = (num_parts*num_nodes)-num_parts

        graph = nx.erdos_renyi_graph(num_nodes,0.5)
        adjlist = graph.adjacency_list()

        print '==============================================================='
        print 'RUNNING TEST %s OF %s TOTAL TESTS' % (test_num, totaltests)
        print '==============================================================='
        print 'TEST PARAMETERS:'
        print 'Elements in Mesh: %s, Number of Partitions: %s' % (num_nodes, num_parts)
        print 'Number of required fully connected qubits:', req_qubits
        print 'Running the %s solver' % type_solver
        print '==============================================================='
        print 'BEGIN TESTING....'

        answers = mosespartitioner(num_parts,adjlist, solver_type = 'isakov')
        print 'Results received in time %.2f seconds' % (answers['timing'])

        # storing results in dict
        results[(num_nodes,num_parts)] = answers

        print '==============================================================='
        if test_num == totaltests:
            print 'TESTING COMPLETE!, now for some postprocessing fun...'
        else:
            test_num += 1
            print 'You beauty! Now onto test number', test_num,'...'
        print '==============================================================='
        print '\n\n'

pickle.dump(results, open("Results_%s.results" % (timestr), "wb"))


##############################################################################
# POST-PROCESSING
# Using functions from qmoseslib.py
##############################################################################

results = pickle.load(open("Results_%s.results" % (timestr), "rb"))

test_num = 1
for num_nodes in range(min_nodes,max_nodes+1):
    for num_parts in range(min_parts,max_parts+1):

        solver = results[(num_nodes,num_parts)]['solver']

        print '==============================================================='
        print 'POST PROCESSING %s OF %s TOTAL TESTS' % (test_num, totaltests)
        print '==============================================================='
        print 'TEST PARAMETERS:'
        print 'Elements in Mesh: %s, Number of Partitions: %s' % (num_nodes, num_parts)
        print 'This test was using the %s solver' % solver
        print '==============================================================='
        print 'INITIATE POST PROCESSING...'

        nodes_listlist = []
        for idx in range(num_parts):
            nodes_listlist.append(results[(num_nodes,num_parts)]['nodelist'])

        # minimised edge condition analysis
        print 'Analysis of minimised edge boundary...'
        adjlist = results[(num_nodes,num_parts)]['adjlist']
        edgelist = adjlist_to_edgelist(adjlist)
        (edges_between, edges_each, total_edges) = edge_resultsanalysis(nodes_listlist, edgelist, num_parts)

        results[(num_nodes,num_parts)]['edges_between'] = edges_between
        results[(num_nodes,num_parts)]['edges_each'] = edges_each
        results[(num_nodes,num_parts)]['total_edges'] = total_edges

        # equal sharing of nodes condition analysis
        print 'Analysis of nodes in each partition...'
        nodeanalysis = num_nodes_per_part(nodes_listlist, num_parts)
        results[(num_nodes,num_parts)]['nodeanalysis'] = nodeanalysis

        # COMPARING AGAINST PYMETIS
        print 'Now comparing against METIS...'
        adjlist_pm = results[(num_nodes,num_parts)]['adjlist']
        pymetis_output = [i + 1 for i in pm.part_graph(num_parts,adjlist_pm)[1]]

        # rearranges pymetis first element to be in partition one
        # "A rose by any other name would smell as sweet"
        # Translation: "A partition by any other value would effectively be the same"

        def replace(list, X, Y):
          i = 0
          for v in list:
             if v == X:
                list.pop(i)
                list.insert(i, Y)
             i += 1

        if pymetis_output[0] != 1:
            A = pymetis_output[0]
            replace(pymetis_output,A,'x')
            replace(pymetis_output,1,A)
            replace(pymetis_output,'x',1)

        nodes_listlist_pm = [[] for x in xrange(num_parts)]
        for idx, element in enumerate(pymetis_output):
            nodes_listlist_pm[element-1].append(idx)

        results[(num_nodes,num_parts)]['nodelist_pm'] = nodes_listlist_pm

        # minimised edge condition analysis
        (edges_between, edges_each, total_edges) = edge_resultsanalysis(nodes_listlist_pm, edgelist, num_parts)
        results[(num_nodes,num_parts)]['edges_between_pm'] = edges_between
        results[(num_nodes,num_parts)]['edges_each_pm'] = edges_each
        results[(num_nodes,num_parts)]['total_edges_pm'] = total_edges

        # equal sharing of nodes condition analysis
        node_analysis = num_nodes_per_part(nodes_listlist_pm,num_parts)
        results[(num_nodes,num_parts)]['nodeanalysis_pm'] = node_analysis

        # calculating Pymetis energy and solution
        # pymetis output to mesh partitioner spin output
        spin_output_pm = []
        for element in pymetis_output:
            spin = [-1]*num_parts
            spin[element-1] = 1
            spin_output_pm.extend(spin)

        # calculating energy
        hvalues = results[(num_nodes,num_parts)]['hvalues']
        J = results[(num_nodes,num_parts)]['jvalues']
        offset = results[(num_nodes,num_parts)]['offset']

        # removing degeneracy
        del spin_output_pm[0:num_parts]

        energy_pm = energy_from_solution(hvalues, J, spin_output_pm, offset=offset)
        results[(num_nodes,num_parts)]['energies_pm'] = energy_pm

        print '==============================================================='
        print '++++++++++++++++++++  INITIAL COMPARISONS  ++++++++++++++++++++'
        print '==============================================================='
        energy_moses = results[(num_nodes,num_parts)]['energies'][0]
        print 'Energy of METIS is %s and MOSES is %s' % (energy_pm, energy_moses)

        print '\n++++++++++++++++++++++++ NODE ANALYSIS ++++++++++++++++++++++++'
        var_pm = results[(num_nodes,num_parts)]['nodeanalysis_pm'][1]
        var_moses = results[(num_nodes,num_parts)]['nodeanalysis'][1]
        print 'Nodes variance of METIS is %s and MOSES is %s' % (var_pm, var_moses)
        print 'Node Lists:',
        print 'METIS:', results[(num_nodes,num_parts)]['nodeanalysis_pm'][0],
        print 'MOSES:', results[(num_nodes,num_parts)]['nodeanalysis'][0]

        print '\n++++++++++++++++++++++++ EDGE ANALYSIS ++++++++++++++++++++++++'
        totedge_moses = results[(num_nodes,num_parts)]['total_edges']
        totedge_pm = results[(num_nodes,num_parts)]['total_edges_pm']
        print 'Total edge boundary of METIS is %s and MOSES is %s' % (totedge_pm, totedge_moses)
        print '==============================================================='

        if test_num == totaltests:
            print 'POST-PROCESSING COMPLETE! File can now be analysed.'
        else:
            test_num += 1
            print 'You beauty! Now post processing test number', test_num,'...'
        print '==============================================================='
        print '\n\n'

pickle.dump(results, open("ResultsAnalysed_%s.results" % (timestr), "wb"))