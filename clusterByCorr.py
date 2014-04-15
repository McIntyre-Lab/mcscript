#!/usr/bin/env python


import numpy as np
import scipy as sc



fname = '/home/jfear/tmp/b52_test.csv';


myCorr = numpy.loadtxt(fname,delimiter=',')


myList = []
for index, row in enumerate(myCorr > .8):
    myList.append(set(np.where(row)[0].tolist()))

results = []
for mySet in myList:









