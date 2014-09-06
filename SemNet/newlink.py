#!/usr/bin/env python
import logging
import numpy as np

def add_additional_links(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Function to take a gene that is already in the model and add additional links. """
    model_type = "Adding additional links"

    # Bind global instance of myCount variable.
    global myCount

    _yvar_ = flatten_list(yvar).split(' ')
    _xvar_ = flatten_list(xvar).split(' ')

    yind = []
    for iso in newGene:
        for index, y in enumerate(_yvar_):
            if iso in y:
                yind.append(index)

    xind = []
    for iso in newGene:
        for index, x in enumerate(_xvar_):
            if iso in x:
                xind.append(index)
    
    if yind:
        # Copy matrices 
        _BETA_ = BETA.copy()

        # Find which columns in the beta matrix
        start = yind[0]
        end = start + len(yind)
        bCols = set(range(start, end))

        # Iterate through rows with 0's and change to 1
        bRows = np.where(_BETA_[:,start] == 0)[0]
        for bRow in bRows:
            if bRow not in bCols:    # For identifiability you cannot act on yourself
                _BETA_[bRow, start:end] = 1
                generate_sas(yvar, xvar, _BETA_, GAMMA, PHI, args)
                output_model_to_log(yvar, xvar, _BETA_, GAMMA, PHI, model_type)
                _BETA_[bRow, start:end] = 0
                myCount +=1

    elif xind:
        # Copy matrices 
        _GAMMA_ = GAMMA.copy()

        # Find which column in the gamma matrix
        start = xind[0]
        end = start + len(xind)

        # Iterate through rows with 0's and change to 1
        gRows = np.where(_GAMMA_[:,start] == 0)[0]
        for gRow in gRows:
            _GAMMA_[gRow, start:end] = 1
            generate_sas(yvar, xvar, BETA, _GAMMA_, PHI, args)
            output_model_to_log(yvar, xvar, BETA, _GAMMA_, PHI, model_type)
            _GAMMA_[gRow, start:end] = 0
            myCount +=1
    else:
        message="This gene {} does not appear in the baseline model, so new links cannot be added".format(newGene)
        logging.warn(message)
        return
