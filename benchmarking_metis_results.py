__author__ = 'John Furness'
##############################################################################
# +++++++++++++++++++++++++++++ METIS  BENCHMARK +++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++ RESULTS ANALYSIS +++++++++++++++++++++++++++++
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
import numpy as np
import scipy.stats
import networkx as nx
import time
import matplotlib.pyplot as plt
import cPickle as pickle

print '==============================================================='
print '++++++++++++++++  METIS BENCHMARK RESULTS  +++++++++++++++++++'
print '==============================================================='

###### OUTLINING CAMPAIGN

filename = "METIS_benchmark_1024_0245.result"
results = pickle.load(open(filename, "rb"))

timestr = results['start_time']

min_vert = results['min_vert']
max_vert = results['max_vert']
step_vert = results['step_vert']

min_part = results['min_part']
max_part = results['max_part']
step_part = results['step_part']

edge_probabilities = results['edge_probabilities']
num_tests = results['num_tests']

total_tests = len(range(min_part, max_part, step_part))*len(range(min_vert, max_vert, step_vert))*len(edge_probabilities)*num_tests
count = 1

print results.keys()

result_analysis = dict()
for no_part in range(min_part, max_part, step_part):
    for no_vert in range(min_vert, max_vert, step_vert):

        Average_Quality = []
        Average_Timing = []

        for edge_prob in edge_probabilities:

            for test in range(num_tests):

                print '==============================================================='
                print 'ANALYSING TEST %s OF %s TOTAL TESTS' % (count, total_tests)
                count += 1
                print '==============================================================='
                print 'Num Parts %s \t No Vert %s \t Edge Prob %s' % (no_part, no_vert, edge_prob)
                print '==============================================================='
                print "Timing", results[(no_part, no_vert, edge_prob, test)]['timing']


                Average_Timing.append(results[(no_part, no_vert, edge_prob, test)]['timing'])
                Average_Quality.append(results[(no_part, no_vert, edge_prob, test)]['quality_factor'])

        Mean_Timing = sum(Average_Timing)/float(len(Average_Timing))
        Std_Timing = np.std(Average_Timing)
        SEM_Timing = scipy.stats.sem(Average_Timing)

        print "Mean_Timing", Mean_Timing
        print "Timing Std", Std_Timing
        print "SEM Quality", SEM_Timing

        Mean_Quality = sum(Average_Quality)/float(len(Average_Quality))
        Std_Quality = np.std(Average_Quality)
        SEM_Quality = scipy.stats.sem(Average_Quality)

        print "Mean Quality", Mean_Quality
        print "Std Quality", Std_Quality
        print "SEM Quality", SEM_Quality

        result_analysis[(no_part, no_vert)] = dict()
        result_analysis[(no_part, no_vert)]['Mean_Timing'] = Mean_Timing
        result_analysis[(no_part, no_vert)]['Std_Timing'] = Std_Timing
        result_analysis[(no_part, no_vert)]['SEM_Timing'] = SEM_Timing

        result_analysis[(no_part, no_vert)]['Mean_Quality'] = Mean_Quality
        result_analysis[(no_part, no_vert)]['Std_Quality'] = Std_Quality
        result_analysis[(no_part, no_vert)]['SEM_Timing'] = SEM_Quality

plt.figure()
plt.title("METIS: Timing Versus Problem Size")
plt.xlabel("Number of Vertices")
plt.ylabel("Time (s)")

for no_part in range(min_part, max_part, step_part):

    plot_x = []
    plot_y = []

    plot_y_error = []

    for no_vert in range(min_vert, max_vert, step_vert):

        plot_x.append(no_vert)

        Mean_Timing = result_analysis[(no_part, no_vert)]['Mean_Timing']
        plot_y.append(Mean_Timing)

        Std_Timing = result_analysis[(no_part, no_vert)]['Std_Timing']
        SEM_Timing = result_analysis[(no_part, no_vert)]['SEM_Timing']
        plot_y_error.append(SEM_Timing)

        Mean_Quality = result_analysis[(no_part, no_vert)]['Mean_Quality']
        Std_Quality = result_analysis[(no_part, no_vert)]['Std_Quality']

    plt.plot(plot_x, plot_y, '-s', label="Parts = %s" % no_part)
    plt.errorbar(plot_x, plot_y, yerr = plot_y_error)
plt.legend(loc ='upper left')
plt.show()

# plt.figure()
# plt.title("METIS: Timing Versus Problem Size")
# plt.errorbar(x, y, xerr=0.2, yerr=0.4)


