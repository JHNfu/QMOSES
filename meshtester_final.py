__author__ = 'John Furness (furnessjw@gmail.com)'
import pymetis as pm
import matplotlib.pyplot as plt
from qmoseslib import *

# INITIALISING ALL RELEVANT PARAMETERS, SOLVERS AND WRITE FILES
# Parameters
two_pow_n_parts = 1
degree = 3
nodes = 100

# Solvers
type_solver = 'isakov'
solver_address = r"C:\Users\John\PycharmProjects\public-isakov\public-isakov\build\visual_studio_2013\bin\an_ss_ge_fi_vdeg_omp.exe"

# Write File
import time
timestr = time.strftime("%Y%m%d_%H%M%S")

results_file = open("MeshResults_%s.results" % (timestr), "w")
print >>results_file, ("Results Mesh Tester - Test: %s" % (timestr))
print >>results_file, ("Run\tPartitions\tNodes\tDegree\tSolver"
                       "\tType\tTotalEdges\tNodeList\tVariance"
                       "\tType\tTotalEdges\tNodeList\tVariance")

run_num = 1
for two_pow_n_parts in range(1,3):

    for degree in range(3,6):
        n = two_pow_n_parts
        m = 2**n

        print "================================================================"
        print "RUNNING CONDITIONS for run number %s....." % (run_num)
        print "Number of Partitions:", 2**n,
        print "\t Number of Nodes:", nodes
        print "\t Degree of Nodes:", degree,
        print "\t Solver Name:", type_solver
        print "================================================================"
        # the 0 <= d < n inequality must be satisfied
        # n * d must be even
        G = nx.random_regular_graph(degree,100)
        adjlist1 = G.adjacency_list()

        opt_sol = recursive_bisection(n,adjlist1, solver_type=type_solver, isakov_address=solver_address)

        nodes_list2 = opt_sol['nodes']
        list2 = qmosesnodelist_to_pymetisoutput(nodes_list2)

        # PRINT NODELIST FOR EACH PARTITION
        nodes_list = []
        for idx in range(2**n):
            nodes_list.append(opt_sol['nodes'][idx][0])

        mesh = nx.Graph()
        edgelist1 = adjlist_to_edgelist(adjlist1)
        mesh.add_edges_from(edgelist1)
        mesh.add_nodes_from(range(len(adjlist1)))

        # removing duplicate edges
        for edge in edgelist1:
            if (edge[1],edge[0]) in edgelist1:
                edgelist1.remove((edge[1],edge[0]))

        def num_per_part(list_of_lists_of_nodes,num_parts):
            num_nodes = []
            for part_1 in range(num_parts):
                num_nodes.append(len(list_of_lists_of_nodes[part_1]))
            return num_nodes

        num_nodes = num_per_part(nodes_list,2**n)
        (edge_results, total_edges) = num_edges_per_partition(nodes_list, edgelist1, 2**n)
        #PYMETIS OUTPUT TO NODE LIST OF LISTS
        node_list_pymetis = [[] for x in xrange(2**n)]
        pymetis_output = pm.part_graph(m,adjlist1)[1]
        for idx in range(len(adjlist1)):
            for part in range(2**n):
                if pymetis_output[idx] == part:
                    node_list_pymetis[part].append(idx)

        ## CREATING RESULTS FOR PYMETIS
        (edge_results_pymetis, total_edges_pymetis) = num_edges_per_partition(node_list_pymetis, edgelist1, 2**n)
        num_nodes_pymetis = num_per_part(node_list_pymetis,2**n)

        def var(any_list):
            mean = sum(any_list)/len(any_list)
            var = 0
            for i in any_list:
                var += (mean - i) ** 2
            return var

        print "\nQMOSES RECURSIVE BISECTION RESULTS:"
        print "Minimising Edge Boundary:"
        print '\t Edges between partitions:', edge_results
        print '\t Total number of edges between partitions:', total_edges/2
        print '\t TOTAL EDGES:', len(edgelist1)
        print 'Same Nodes Per Partition:'
        print '\t Number of Nodes in Each Partition:', num_nodes
        print '\t Variance:', var(num_nodes)

        print "\nPYMETIS RESULTS:"
        print "Minimising Edge Boundary:"
        print '\t Edges between partitions:', edge_results_pymetis
        print '\t Total number of edges between partitions:', total_edges_pymetis/2
        print '\t TOTAL EDGES:', len(edgelist1)
        print 'Same Nodes Per Partition:'
        print '\t Number of Nodes in Each Partition:', num_nodes_pymetis
        print '\t Variance:', var(num_nodes_pymetis)

        print "Writing results to data file...",

        # Running Conditions
        print >>results_file, "%s\t%s\t%s\t%s\t%s" % (run_num, 2**n, nodes, degree, type_solver),
        # QMOSES Results
        print >>results_file, "\tQMOSES\t%s\t%s\t%s" % (total_edges/2, num_nodes, var(num_nodes)),
        # PYMETIS RESULTS
        print >>results_file, "\tPYMETIS\t%s\t%s\t%s" % (total_edges_pymetis/2, num_nodes_pymetis, var(num_nodes_pymetis))

        print "DONE!\n\nNow for run %s ...." % (run_num + 1)
        run_num += 1

print "\nTEST STOPPED"

