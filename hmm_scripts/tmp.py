#!/usr/bin/env python

from pomegranate import *

dists = [NormalDistribution(5, 1), NormalDistribution(1, 7)]
trans_mat = numpy.array([[0.7, 0.3],
                         [0.2, 0.8])

starts = numpy.array([1,0])
ends   = numpy.array([1,0])

sim_model = HiddenMarkovModel.from_matrix(trans_mat, dists, starts, ends)

