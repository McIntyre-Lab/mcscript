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

    # Initialize default values
    path.reinit()

    # Expand gamma and beta to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_beta()
        path.expand_gamma()

    # Make copy of oiriginal for iterating
    _yvar = list(path.yvar)

    # Append new gene to yvar
    path.yvar.append(newgene.name)

    # Iterate through BETA matrix and output all possible models
    index = 0
    for gene in _yvar:
        # How many isoforms does the current gene have
        isoCnt = utils.isoCount(gene)

        # Create slice index of bottom 'newgene.count' rows of beta and
        # alterating columns of beta
        bSlice = np.s_[-newgene.count:, index:index+isoCnt]

        model_type = "Adding Genes Downstream of Beta:\n {0} -> {1}".format(gene, newgene.name)

        path.beta[bSlice] = 1
        createOutput(path, model_type, args)
        path.beta[bSlice] = 0
        path.count_increment()
        index += isoCnt

def add_genes_ds_gamma(path, newgene, args):
    """ Iterate through gamma matrix and place new genes downstream of each gamma. """

    # Initialize default values
    path.reinit()

    # Expand gamma and beta to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_beta()
        path.expand_gamma()

    # Append new gene to yvar
    path.yvar.append(newgene.name)

    # Iterate through GAMMA matrix and output all possible models
    index = 0
    for gene in path.xvar:
        # How many isoforms does the current gene have
        isoCnt = utils.isoCount(gene)

        # Create slice index of bottom 'newgene.count' rows of gamma and
        # alternating columns of gamma
        gSlice = np.s_[-newgene.count:, index:index+isoCnt]

        model_type = "Adding Genes Downstream of Gamma:\n {0} -> {1}".format(gene, newgene.name)

        path.gamma[gSlice] = 1
        createOutput(path, model_type, args)
        path.gamma[gSlice] = 0
        path.count_increment()
        index += isoCnt

# Add genes upstream of nodes
def add_genes_above_beta(path, newgene, args):
    """ Iterate through gamma matrix and place new genes upstream of betas. """

    # Initialize default values
    path.reinit()

    # Expand matrices GAMMA and PHI matrices
    for iso in range(0, newgene.count):
        path.expand_gamma(axis=1) # Add column not row
        path.expand_phi()

    # Append new gene to yvar
    path.xvar.append(newgene.name)

    index = 0
    for gene in path.yvar:
        # How many isoforms does the current gene have
        isoCnt = utils.isoCount(gene)

        # Create slice index of right most 'newgene.count' columns of gamma and
        # alterating rows of gamma
        gSlice = np.s_[index:index+isoCnt, -newgene.count:]

        model_type = "Adding Genes Upstream of beta:\n {1} -> {0}".format(gene, newgene.name)

        path.gamma[gSlice] = 1
        createOutput(path, model_type, args)
        path.gamma[gSlice] = 0
        path.count_increment()
        index += isoCnt

def add_genes_above_gamma(path, newgene, args):
    """ Iterate through gamma matrix and place new genes upstream of gammas. """

    # Initialize default values
    path.reinit()

    # Make copy of oiriginal for iterating
    _xvar = list(path.xvar)

    for gene in _xvar:
        # How many isoforms does the current gene have
        isoCnt = utils.isoCount(gene)

        # Convert the current exogenous gene to an enogenous gene
        path.convert_ExogToEndog(gene)

        # Append new gene to xvar
        path.xvar.append(newgene.name)

        # Expand gamma and phi matrices for my new gene
        for iso in range(0, newgene.count):
            path.expand_gamma(axis=1) # Add column not row
            path.expand_phi()

        # Create slice index of right most 'newgene.count' columns of gamma and
        # the bottom 'isoCnt' rows of gamma.
        gSlice = np.s_[-isoCnt:, -newgene.count:]

        model_type = "Adding Genes Upstream of Gamma:\n {1} -> {0}".format(gene, newgene.name)

        path.gamma[gSlice] = 1
        createOutput(path, model_type, args)
        path.count_increment()

        # Reset I need to move different exogenous genes 
        path.reinit()

# Add genes between connected nodes
def add_genes_bt_beta(path, newgene, args):
    """ Iterate through beta matrix and add genes between two betas. """

    # Initialize default values again
    path.reinit()

    # Expand gamma and beta to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_beta()
        path.expand_gamma()

    # Make copy of oiriginal for iterating
    _yvar = list(path.yvar)

    # Append new gene to yvar
    path.yvar.append(newgene.name)

    # Iterate through BETA matrix and add genes between betas
    cIndex = 0
    for col in _yvar:
        # How many isoforms does the current column have
        colCnt = utils.isoCount(col)

        rIndex = 0
        for row in _yvar:
            # A gene cannot act on itself, so skip when row=col
            if col == row:
                continue

            # How many isoforms does the current row have
            rowCnt = utils.isoCount(row)

            # Zero slice, is the current location in Beta, if there are 1s
            # there I will turn them 0
            zSlice = np.s_[rIndex:rIndex+rowCnt, cIndex:cIndex+colCnt]

            # Target slice, is the bottom 'newgene.count' rows, that will get turned 1
            # showing that current 'col' is acting on the newgene
            tSlice = np.s_[-newgene.count:, cIndex:cIndex+colCnt]

            # Source slice, is the right most 'newgene.count' cols, that will get turned 1
            # showing that newgene is acting on the current row
            sSlice = np.s_[rIndex:rIndex+rowCnt, -newgene.count:]

            # If there was a previous interaction (i.e. 1) then make models
            if np.all(path.beta[zSlice]):
                model_type = "Adding Genes Between Betas:\n {0} -> {2} -> {1}".format(col, row, newgene.name)

                path.beta[zSlice] = 0
                path.beta[tSlice] = 1
                path.beta[sSlice] = 1
                createOutput(path, model_type, args)
                path.beta[zSlice] = 1
                path.beta[tSlice] = 0
                path.beta[sSlice] = 0
                path.count_increment()

            rIndex += rowCnt
        cIndex += colCnt

def add_genes_bt_gamma_beta(path, newgene, args):
    """ Iterate through gamma matrix and place new genes between gamma and ds beta. """

    # Initialize default values again
    path.reinit()

    # Expand gamma and beta to account for new isoforms
    for iso in range(0, newgene.count):
        path.expand_beta()
        path.expand_gamma()

    # Make copy of oiriginal for iterating
    _yvar = list(path.yvar)

    # Add newGene to yvar because they are new endogenous genes
    path.yvar.append(newgene.name)

    # Iterate through each row of GAMMA matrix and place new gene in between
    # gamma and beta
    cIndex = 0
    for col in path.xvar:
        # How many isoforms does the current column have
        colCnt = utils.isoCount(col)

        rIndex = 0
        for row in _yvar:
            # How many isoforms does the current row have
            rowCnt = utils.isoCount(row)

            # Zero slice, is the current location in Gamma, if there are 1s
            # there I will turn them 0
            zSlice = np.s_[rIndex:rIndex+rowCnt, cIndex:cIndex+colCnt]

            # Target slice, is the bottom 'newgene.count' rows in GAMMA, that
            # will get turned 1 showing that current 'col' is acting on the
            # newgene
            tSlice = np.s_[-newgene.count:, cIndex:cIndex+colCnt]

            # Source slice, is the right most 'newgene.count' cols in BETA,
            # that will get turned 1 showing that newgene is acting on the
            # current row
            sSlice = np.s_[rIndex:rIndex+rowCnt, -newgene.count:]

            # If there was a previous interaction (i.e. 1) then make models
            if np.all(path.gamma[zSlice]):
                model_type = "Adding Genes Between Gamma and Beta:\n {0} -> {2} -> {1}".format(col, row, newgene.name)

                path.gamma[zSlice] = 0
                path.gamma[tSlice] = 1
                path.beta[sSlice] = 1
                createOutput(path, model_type, args)
                path.gamma[zSlice] = 1
                path.gamma[tSlice] = 0
                path.beta[sSlice] = 0
                path.count_increment()
            rIndex += rowCnt
        cIndex += colCnt
