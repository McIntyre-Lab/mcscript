#!/usr/bin/env python
import numpy as np

def isoCount (myVars):
    if hasattr(myVars, '__iter__'):
        cnt = len(myVars)
    else:
        cnt = 1
    return cnt

def flatten_list(myVars, return_list=0):
    """ Isoforms are being treated together as a group. To delineate them,
    isoforms are grouped together as tuples. For the calis statement, xvar and
    yvar need to each be a string. This function flattens these tuples into a
    single list and creates a string. """

    result = []
    for element in myVars: 
        if hasattr(element, '__iter__'):
            # if tuple flatten and append
            result.extend(element)
        else:
            # if not a tuple just append
            result.append(element)
    if return_list:
        return result
    else:
        return(' '.join(result))
