__author__ = 'John Furness'
##############################################################################
# Nested Dissection Method Tester
# September 2015
##############################################################################

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