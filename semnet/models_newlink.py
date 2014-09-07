#!/usr/bin/env python
import logging
from os.path import splitext
import numpy as np
from semnet.writer import createOutput
from semnet import utils

################################################################################
# Add additional links to existing genes in the network
################################################################################
def add_additional_links(path, newgene, args):
    """ Function to take a gene that is already in the model and add additional links. """
    model_type = "Adding additional links"

    # Initialize default values
    path.reinit()

    # Create a flattend list
    _yvar = utils.flatten_list(path.yvar).split(' ')
    _xvar = utils.flatten_list(path.xvar).split(' ')

    yind = []
    for iso in newgene.name:
        for index, y in enumerate(_yvar):
            if iso in y:
                yind.append(index)

    xind = []
    for iso in newgene.name:
        for index, x in enumerate(_xvar):
            if iso in x:
                xind.append(index)
    
    if yind:
        # Find which columns in the beta matrix
        start = yind[0]
        end = start + len(yind)
        _bCols = set(range(start, end))

        # Iterate through rows with 0's and change to 1
        _bRows = np.where(path.beta[:,start] == 0)[0]
        for row in _bRows:
            if row not in _bCols:    # For identifiability you cannot act on yourself
                path.beta[row, start:end] = 1
                createOutput(path, model_type, args)
                path.beta[row, start:end] = 0
                path.count_increment()

    elif xind:
        # Find which column in the gamma matrix
        start = xind[0]
        end = start + len(xind)

        # Iterate through rows with 0's and change to 1
        _gRows = np.where(path.gamma[:,start] == 0)[0]
        for row in _gRows:
            path.gamma[row, start:end] = 1
            createOutput(path, model_type, args)
            path.gamma[row, start:end] = 0
            path.count_increment()
    else:
        message="This gene {} does not appear in the baseline model, so new links cannot be added".format(newgene.name)
        logging.warn(message)
        return
