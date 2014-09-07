#!/usr/bin/env python
def flatten_list(myVars):
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
    return(' '.join(result))
