__author__ = 'John Furness'
##############################################################################
# +++++++++++++++++++++++++++++ METIS  BENCHMARK +++++++++++++++++++++++++++++
##############################################################################
# John Furness
# October 2015
##############################################################################
# This code runs a very simple benchmark of the timing of METIS using the
# Erdos Renyi graph generator. METIS will be used as a graph partitioner
# for varying connectivity and graph sizes.
# The main goal will be to reduce the edge boundary between partitions,
# this test will only consider 3 and four partitions.
##############################################################################

from qmoseslib import *
import pymetis as pm
import networkx as nx
import time
import cPickle as pickle

print '==============================================================='
print '++++++++++++++++  METIS BENCHMARKER TESTER  +++++++++++++++++++'
print '==============================================================='

###### OUTLINING CAMPAIGN

timestr = time.strftime("%m%d_%H%M")
filename = "METIS_benchmark_%s.result" % timestr
# results = pickle.load(open(filename, "rb"))

min_vert = 10
max_vert = 1000
step_vert = 20

min_part = 2
max_part = 4
step_part = 1

edge_probabilities = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

#number of tests each random problem
num_tests = 50

results = dict()

results['start_time'] = timestr

results['min_vert'] = min_vert
results['max_vert'] = max_vert
results['step_vert'] = step_vert

results['min_part'] = min_part
results['max_part'] = max_part
results['step_part'] = step_part

results['edge_probabilities'] = edge_probabilities
results['num_tests'] = num_tests

pickle.dump(results, open(filename, "wb"))

total_tests = len(range(min_part, max_part, step_part))*len(range(min_vert, max_vert, step_vert))*len(edge_probabilities)*num_tests
count = 1
for no_part in range(min_part, max_part, step_part):
    for no_vert in range(min_vert, max_vert, step_vert):
        for edge_prob in edge_probabilities:

            rand_graph = nx.erdos_renyi_graph(no_vert,edge_prob)
            adjlist = rand_graph.adjacency_list()
            edgelist = adjlist_to_edgelist(adjlist)

            # generating random solution for measure of quality
            rand_solution = generate_RandomSol(no_part, no_vert)

            rand_solution_list = [[] for i in range(no_part+1)]
            for idx, result in enumerate(rand_solution):
                rand_solution_list[result-1].append(idx)

            # analysis of the quality of random results
            (rand_node_res, rand_node_var) = num_nodes_per_part(rand_solution_list, no_part)
            (rand_edges_between, rand_edges_each, rand_total_edges) = edge_resultsanalysis(rand_solution_list, edgelist, no_part)

            test_results = dict()

            test_results['rand_total_edges'] = rand_total_edges
            test_results['rand_edges_each'] = rand_edges_each
            test_results['rand_edges_between'] = rand_edges_between
            test_results['rand_node_res'] = rand_node_res
            test_results['rand_node_var'] = rand_node_var

            for test in range(num_tests):

                # for loading previous result
                # if results[(no_part, no_vert, edge_prob, test)] != {}:
                #     print "Test %s %s %s %s has been completed, continuing... " % (no_part, no_vert, edge_prob, test)
                #     count += 1
                #     continue


                print '==============================================================='
                print 'RUNNING TEST %s OF %s TOTAL TESTS' % (count, total_tests)
                print '==============================================================='
                print 'Num Parts %s \t No Vert %s \t Edge Prob %s' % (no_part, no_vert, edge_prob)
                print '==============================================================='

                print "RUNNING METIS TEST...",

                t0_m = time.time()
                metis_sol = pm.part_graph(no_part,adjlist)[1]
                t1_m = time.time()
                metis_solution = [idx + 1 for idx in metis_sol]
                t_final = t1_m - t0_m

                print "Done!"
                print "Timing: %s" % t_final

                metis_solution_list = [[] for i in range(no_part+1)]
                for idx, result in enumerate(metis_solution):
                    metis_solution_list[result-1].append(idx)

                # analysis of quality of metis solution

                (metis_node_res, metis_node_var) = num_nodes_per_part(metis_solution_list, no_part)
                (metis_edges_between, metis_edges_each, metis_total_edges) = edge_resultsanalysis(metis_solution_list, edgelist, no_part)

                quality_factor = 1 - metis_total_edges/rand_total_edges
                print "Quality of Solution is %s" % quality_factor

                # saving results

                test_results['quality_factor'] = quality_factor

                test_results['metis_total_edges'] = metis_total_edges
                test_results['metis_edges_each'] = metis_edges_each
                test_results['metis_edges_between'] = metis_edges_between
                test_results['metis_node_res'] = metis_node_res
                test_results['metis_node_var'] = metis_node_var

                test_results['timing'] = t_final

                results[(no_part, no_vert, edge_prob, test)] = test_results

                count += 1

        print "Dumping results...",
        pickle.dump(results, open(filename, "wb"))
        print "Dumped!"

print "Dumping final results...",
pickle.dump(results, open(filename), "wb")
print "Done!"

print "\n TESTING COMPLETE!"





