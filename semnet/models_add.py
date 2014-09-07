#!/usr/bin/env python
import logging
from os.path import splitext
import numpy as np
from semnet.writer import createOutput
from semnet import utils

################################################################################
# Add genes to the network
################################################################################

# Add genes downstream of nodes
def add_genes_ds_beta(path, newgene, args):
    """ Iterate through beta matrix and place new genes downstream of each beta. """
    model_type = "Adding Genes Downstream of Beta"

    # Initialize default values
    path.reinit()

    # Expand beta to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_beta()

    path.yvar.append(newgene.name)
    # Iterate through BETA matrix and output all possible models
    for i in range(0, path.bCol):
        path.beta[path.bRow:, i] = 1
        createOutput(path, model_type, args)
        path.beta[path.bRow:, i] = 0
        path.count_increment()

def add_genes_ds_gamma(path, newgene, args):
    """ Iterate through gamma matrix and place new genes downstream of each gamma. """
    model_type = "Adding Genes Downstream of Gamma"

    # Initialize default values
    path.reinit()

    # Expand gamma to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_gamma()

    path.yvar.append(newgene.name)
    # Iterate through GAMMA matrix and output all possible models
    index = 0
    for gene in path.xvar:
        if hasattr(gene, '__iter__'):
            xvarCount = len(gene)
        else:
            xvarCount = 1

        start = index
        end = start + xvarCount
        index = end

        path.gamma[path.bRow:, start:end] = 1
        createOutput(path, model_type, args)
        path.gamma[path.bRow:, start:end] = 0
        path.count_increment()

# Add genes upstream of nodes
def add_genes_above_beta(path, newgene, args):
    """ Iterate through gamma matrix and place new genes upstream of betas. """
    model_type = "Adding Genes upstream of Beta"

    # Initialize default values
    path.reinit()
    path.xvar.append(newgene.name)

    index = 0
    for gene in path.yvar:
        # determine if gene has multiple isoforms
        if hasattr(gene,'__iter__'):
            yvarCount = len(gene)
        else:
            yvarCount = 1

        yvarStart = index
        yvarEnd = yvarStart + yvarCount
        coords = (yvarStart, yvarEnd)
        index = yvarEnd

        # Initialize default values
        path.reinit()

        # Edit variable list
        path.xvar.append(newgene.name)

        # Expand matrices GAMMA and PHI matrices
        for iso in range(0, newgene.count):
            # Adding new exogenous genes
            path.expand_gamma(axis=1)
            path.expand_phi()

        # Fill in 1's for GAMMA and PHI
        _gRow, _gCol = np.shape(path.gamma)
        path.gamma[yvarStart:yvarEnd, (_gCol - newgene.count):] = 1

        _pRow, _pCol = np.shape(path.phi)
        path.phi[(_pRow - newgene.count):, (_pCol - newgene.count):] = 1

        # Build SAS output
        createOutput(path, model_type, args)
        path.count_increment()

def add_genes_above_gamma(path, newgene, args):
    """ Iterate through gamma matrix and place new genes upstream of gammas. """
    model_type = "Adding Genes upstream of Gamma"

    # Initialize default values
    path.reinit()

    index = 0
    for xloc, gene in enumerate(path.xvar):
        # determine if gene has multiple isoforms
        if hasattr(gene,'__iter__'):
            xvarCount = len(gene)
        else:
            xvarCount = 1

        xvarStart = index
        xvarEnd = xvarStart + xvarCount
        coords = (xvarStart, xvarEnd)
        index = xvarEnd

        # Initialize default values again
        path.reinit()

        # Expand matrices
        for iso in range(0, newgene.count):
            # Adding new exogenous genes
            path.expand_gamma(axis=1)
            path.expand_phi()

        for iso in range(0, xvarCount):
            path.expand_beta()
            path.expand_gamma()

        # Fill in 1's for GAMMA and PHI
        _gRow, _gCol = np.shape(path.gamma)
        path.gamma[(_gRow - xvarCount):, (_gCol - newgene.count):] = 1

        _pRow, _pCol = np.shape(path.phi)
        path.phi[(_pRow - newgene.count):, (_pCol - newgene.count):] = 1

        # Fill in 1's for BETA. Need to use locations information from GAMMA.
        # To do this, slice out the current column from GAMMA identify
        # locations of 1's and add to end of BETA 
        for gammaCol in range(*coords):
            gSlice = path.gamma[:, gammaCol]
            for nPos in np.flatnonzero(gSlice):
                path.beta[nPos, path.bCol:] = 1

        # Delete appropriate columns/rows from GAMMA and PHI
        path.del_col_gamma(coords)
        path.del_row_col_phi(coords)

        # Edit variable list
        del path.xvar[xloc]
        path.xvar.append(newgene.name)

        path.yvar.append(gene)

        # Build SAS output
        createOutput(path, model_type, args)
        path.count_increment()

# Add genes between connected nodes
def add_genes_bt_beta(path, newgene, args):
    """ Iterate through beta matrix and add genes between two betas. """
    model_type = "Adding Genes between Betas"

    # Initialize default values again
    path.reinit()

    # Add newGene to yvar because they are new endogenous genes
    path.yvar.append(newgene.name)

    # Expand beta to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_beta()

    rows, cols = np.nonzero(path.beta)
    for row, col in zip(rows, cols):
        path.beta[row, col] = 0
        path.beta[row, path.bCol:] = 1
        path.beta[path.bRow:, col] = 1
        createOutput(path, model_type, args)
        path.beta[row, col] = 1
        path.beta[row, path.bCol:] = 0
        path.beta[path.bRow:, col] = 0
        path.count_increment()

def add_genes_bt_gamma_beta(path, newgene, args):
    """ Iterate through gamma matrix and place new genes between gamma and ds beta. """
    model_type = "Adding Genes between Gamma and Beta"

    # Initialize default values again
    path.reinit()

    # Add newGene to yvar because they are new endogenous genes
    path.yvar.append(newgene.name)

    # Expand gamma and beta to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_beta()
        path.expand_gamma()

    # Iterate through each row of GAMMA matrix and place new gene in between
    # gamma and beta
    for row_index in range(0, path.gRow):
        col_index = 0
        path.beta[row_index, path.bCol:] = 1
        for gene in path.xvar:
            if hasattr(gene,'__iter__'):
                xvarCount = len(gene)
            else:
                xvarCount = 1
            start = col_index
            end = col_index + xvarCount
            col_index = end

            if np.all(path.gamma[row_index, start:end] == 1):
                path.gamma[path.gRow:, start:end] = 1
                path.gamma[row_index, start:end] = 0
                createOutput(path, model_type, args)
                path.gamma[path.gRow:, start:end] = 0
                path.gamma[row_index, start:end] = 1
                path.count_increment()
            else:
                logging.debug("The gene: {0} does not act on {1}".format(gene,path.yvar[row_index]))
        path.beta[row_index, path.bCol:] = 0
