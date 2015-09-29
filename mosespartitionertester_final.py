__author__ = 'John'
import cPickle as pickle
from qmoseslib import *
import isakovlib as isakov
import networkx as nx
# Solvers
type_solver = 'isakov'
solver_address = r"C:\Users\John\PycharmProjects\public-isakov\public-isakov\build\visual_studio_2013\bin\an_ss_ge_fi_vdeg_omp.exe"

total_nodes = 100

for num_nodes in range(10,total_nodes):
    for num_parts in range(2,4):

        graph = nx.erdos_renyi_graph(num_nodes,0.5)
        adjlist = graph.adjacency_list()

        answers = mosespartitioner(num_parts,adjlist, solver_type = 'isakov')

