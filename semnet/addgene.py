#!/usr/bin/env python
import logging
import numpy as np

# Add genes downstream of nodes
def add_genes_ds_beta(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through beta matrix and place new genes downstream of each beta. """
    model_type = "Adding Genes Downstream of Beta"

    # Bind global instance of myCount variable.
    global myCount

    # Add newGene to yvar because they are new endogenous genes
    _yvar_ = list(yvar)
    _yvar_.append(newGene)

    # Copy matrices 
    _BETA_ = BETA.copy()

    bRow, bCol = np.shape(_BETA_)

    # Expand beta to account for new isoforms
    for iso in range(0, newCount):
        _BETA_ = expand_beta(_BETA_)

    # Iterate through BETA matrix and output all possible models
    for i in range(0, bCol):
        _BETA_[bRow:, i] = 1
        generate_sas(_yvar_, xvar, _BETA_, GAMMA, PHI, args)
        output_model_to_log(_yvar_, xvar, _BETA_, GAMMA, PHI, model_type)
        myCount +=1
        _BETA_[bRow:, i] = 0

def add_genes_ds_gamma(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through gamma matrix and place new genes downstream of each gamma. """
    model_type = "Adding Genes Downstream of Gamma"

    # Bind global instance of myCount variable.
    global myCount

    # Add newGene to yvar because they are new endogenous genes
    _yvar_ = list(yvar)
    _yvar_.append(newGene)

    # Copy matrices 
    _BETA_ = BETA.copy()
    _GAMMA_ = GAMMA.copy()

    bRow, bCol = np.shape(_BETA_)

    # Expand gamma to account for new isoforms
    for iso in range(0, newCount):
        _GAMMA_ = expand_gamma(_GAMMA_)

    # Iterate through GAMMA matrix and output all possible models
    index = 0
    for gene in xvar:
        if hasattr(gene, '__iter__'):
            xvarCount = len(gene)
        else:
            xvarCount = 1

        start = index
        end = start + xvarCount
        index = end

        _GAMMA_[bRow:, start:end] = 1
        generate_sas(_yvar_, xvar, _BETA_, _GAMMA_, PHI, args)
        output_model_to_log(_yvar_, xvar, _BETA_, _GAMMA_, PHI, model_type)
        _GAMMA_[bRow:, start:end] = 0
        myCount +=1

# Add genes upstream of nodes
def add_genes_above_gamma(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through gamma matrix and place new genes upstream of gammas. """
    model_type = "Adding Genes upstream of Gamma"

    # Bind global instance of myCount variable.
    global myCount

    index = 0
    for xloc, gene in enumerate(xvar):
        # determine if gene has multiple isoforms
        if hasattr(gene,'__iter__'):
            xvarCount = len(gene)
        else:
            xvarCount = 1

        xvarStart = index
        xvarEnd = xvarStart + xvarCount
        index = xvarEnd

        # Copy matrices 
        _BETA_ = BETA.copy()
        _GAMMA_ = GAMMA.copy()
        _PHI_ = PHI.copy()

        bRow, bCol = np.shape(_BETA_)

        # Expand matrices
        for iso in range(0, newCount):
            # Adding new exogenous genes
            _GAMMA_ = expand_gamma(_GAMMA_, axis=1)
            _PHI_ = expand_phi(_PHI_)

        for iso in range(0, xvarCount):
            _BETA_ = expand_beta(_BETA_)
            _GAMMA_ = expand_gamma(_GAMMA_)

        # Fill in 1's for GAMMA and PHI
        gRow, gCol = np.shape(_GAMMA_)
        _GAMMA_[(gRow - xvarCount):, (gCol - newCount):] = 1

        pRow, pCol = np.shape(_PHI_)
        _PHI_[(pRow - newCount):, (pCol - newCount):] = 1

        # Fill in 1's for BETA. Need to use locations information from GAMMA.
        # To do this, slice out the current column from GAMMA identify
        # locations of 1's and add to end of BETA 
        for gammaCol in range(xvarStart, xvarEnd):
            gSlice = _GAMMA_[:, gammaCol]
            for nPos in np.flatnonzero(gSlice):
                _BETA_[nPos, bCol:] = 1

        # Delete appropriate columns/rows from GAMMA and PHI
        _GAMMA_ = np.delete(_GAMMA_, range(xvarStart, xvarEnd), 1)
        _PHI_ = np.delete(_PHI_, range(xvarStart, xvarEnd), 0)
        _PHI_ = np.delete(_PHI_, range(xvarStart, xvarEnd), 1)

        # Edit variable list
        _xvar_ = list(xvar)
        del _xvar_[xloc]
        _xvar_.append(newGene)

        _yvar_ = list(yvar)
        _yvar_.append(gene)

        # Build SAS output
        generate_sas(_yvar_, _xvar_, _BETA_, _GAMMA_, _PHI_, args)
        output_model_to_log(_yvar_, _xvar_, _BETA_, _GAMMA_, _PHI_, model_type)
        myCount +=1

def add_genes_above_beta(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through gamma matrix and place new genes upstream of betas. """
    model_type = "Adding Genes upstream of Beta"

    # Bind global instance of myCount variable.
    global myCount

    # Edit variable list
    _xvar_ = list(xvar)
    _xvar_.append(newGene)

    index = 0
    for gene in yvar:
        # determine if gene has multiple isoforms
        if hasattr(gene,'__iter__'):
            yvarCount = len(gene)
        else:
            yvarCount = 1

        yvarStart = index
        yvarEnd = yvarStart + yvarCount
        index = yvarEnd

        # Copy matrices 
        _BETA_ = BETA.copy()
        _GAMMA_ = GAMMA.copy()
        _PHI_ = PHI.copy()

        # Expand matrices GAMMA and PHI matrices
        for iso in range(0, newCount):
            # Adding new exogenous genes
            _GAMMA_ = expand_gamma(_GAMMA_, axis=1)
            _PHI_ = expand_phi(_PHI_)

        # Fill in 1's for GAMMA and PHI
        gRow, gCol = np.shape(_GAMMA_)
        _GAMMA_[yvarStart:yvarEnd, (gCol - newCount):] = 1

        pRow, pCol = np.shape(_PHI_)
        _PHI_[(pRow - newCount):, (pCol - newCount):] = 1

        # Build SAS output
        generate_sas(yvar, _xvar_, _BETA_, _GAMMA_, _PHI_, args)
        output_model_to_log(yvar, _xvar_, _BETA_, _GAMMA_, _PHI_, model_type)
        myCount +=1

# Add genes between connected nodes
def add_genes_bt_gamma_beta(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through gamma matrix and place new genes between gamma and ds beta. """
    model_type = "Adding Genes between Gamma and Beta"

    # Bind global instance of myCount variable.
    global myCount

    # Add newGene to yvar because they are new endogenous genes
    _yvar_ = list(yvar)
    _yvar_.append(newGene)

    # Copy matrices 
    _BETA_ = BETA.copy()
    _GAMMA_ = GAMMA.copy()

    bRow, bCol = np.shape(_BETA_)
    gRow, gCol = np.shape(_GAMMA_)

    # Expand gamma and beta to account for new isoforms
    for iso in range(0, newCount):
        _BETA_ = expand_beta(_BETA_)
        _GAMMA_ = expand_gamma(_GAMMA_)

    # Iterate through each row of GAMMA matrix and place new gene in between
    # gamma and beta
    for row_index in range(0, gRow):
        col_index = 0
        _BETA_[row_index, bCol:] = 1
        for gene in xvar:
            if hasattr(gene,'__iter__'):
                xvarCount = len(gene)
            else:
                xvarCount = 1
            start = col_index
            end = col_index + xvarCount
            col_index = end
            if np.all(_GAMMA_[row_index, start:end] == 1):
                _GAMMA_[gRow:, start:end] = 1
                _GAMMA_[row_index, start:end] = 0
                generate_sas(_yvar_, xvar, _BETA_, _GAMMA_, PHI, args)
                output_model_to_log(_yvar_, xvar, _BETA_, _GAMMA_, PHI, model_type)
                _GAMMA_[gRow:, start:end] = 0
                _GAMMA_[row_index, start:end] = 1
                myCount +=1
            else:
                logging.debug("The gene: {0} does not act on {1}".format(gene,yvar[row_index]))
        _BETA_[row_index, bCol:] = 0

def add_genes_bt_beta(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through beta matrix and add genes between two betas. """
    model_type = "Adding Genes between Betas"

    # Bind global instance of myCount variable.
    global myCount

    # Add newGene to yvar because they are new endogenous genes
    _yvar_ = list(yvar)
    _yvar_.append(newGene)

    # Copy matrices 
    _BETA_ = BETA.copy()

    bRow, bCol = np.shape(_BETA_)

    # Expand beta to account for new isoforms
    for iso in range(0, newCount):
        _BETA_ = expand_beta(_BETA_)

    rows, cols = np.nonzero(_BETA_)
    for row, col in zip(rows, cols):
        _BETA_[row, col] = 0
        _BETA_[row, bCol:] = 1
        _BETA_[bRow:, col] = 1
        generate_sas(_yvar_, xvar, _BETA_, GAMMA, PHI, args)
        output_model_to_log(_yvar_, xvar, _BETA_, GAMMA, PHI, model_type)
        _BETA_[row, col] = 1
        _BETA_[row, bCol:] = 0
        _BETA_[bRow:, col] = 0
        myCount +=1
