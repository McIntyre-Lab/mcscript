#!/usr/bin/env python
import logging
from os.path import splitext
import numpy as np

# Initialize model counter.
myCount = 1

# Functions to expand matrices appropriately
def expand_beta(BETA):
    bRow = np.shape(BETA)[0]
    BETA = np.insert(BETA, bRow, values=0, axis=1)
    BETA = np.insert(BETA, bRow, values=0, axis=0)
    return(BETA)

def expand_gamma(GAMMA,axis=0):
    gPos = np.shape(GAMMA)[axis]
    GAMMA = np.insert(GAMMA, gPos, values=0, axis=axis)
    return(GAMMA)

def expand_phi(PHI):
    pRow = np.shape(PHI)[0]
    PHI = np.insert(PHI, pRow, values=0, axis=0)
    PHI = np.insert(PHI, pRow, values=0, axis=1)
    return(PHI)

# Functions to Add Gene to different locations in the network
def run_baseline(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through beta matrix and place new genes downstream of each beta. """
    model_type = "Baseline"

    # Bind global instance of myCount variable.
    global myCount

    generate_sas(yvar, xvar, BETA, GAMMA, PHI, args, fbase=1)
    output_model_to_log(yvar, xvar, BETA, GAMMA, PHI, model_type)

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

# Functions to Build SAS Output
def build_matrix_output(myMatrix, prefix='', find_zero=0):
    """ Construct proc calis matrix statements from full matrices. Iterate over
        the full matrix and identify the important coordinates and construct a
        string with the correct syntax.
        
        myMatrix is a matrix of {1,0} either Beta, Gamma, or Phi from SEM.

        prefix {beta,gamma,''} is the prefix for the parameter labels, leave blank if none.

        find_zero {0,1} if 0 then identify non-zero indices (for beta and gamma matrices), if 1 identify 0 indices (for covariance matrix phi).

    """

    # pull out the relevant coordinates
    if find_zero:
        # run this if you are looking at the covariance matrix
        rows, cols = np.where(myMatrix == 0)
    else:
        # run this if you are looking at the gamma or beta matrix
        rows, cols = np.where(myMatrix > 0)

    # Build the matrix statements
    myList = list()
    for index,(row,col) in enumerate(zip(rows, cols)):
        if prefix:
            myList.append("[{0},{1}] = {2}{3}".format(row+1,col+1,prefix,index+1))
        else:
            myList.append("[{0},{1}] = 0".format(row+1,col+1))
    return(",\n".join(myList))

def build_calis(yvar, xvar, BETA, GAMMA, PHI, args, fbase):
    """ Write the complete proc calis statement for sas """

    # Bind global instance of myCount variable.
    global myCount

    if fbase:
        modelNum = 'baseline'
    else:
        modelNum = myCount

    # create output file name
    _oname_ = "{0}_{1}.sas".format(splitext(args.oname)[0], modelNum)

    # write everything out
    with open(_oname_,'w') as OUT:
        OUT.write("""\
                        libname SEM "{0}";

                        proc calis data=SEM.{1} method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
                            lismod
                                yvar = {2}, 
                                xvar = {3}
                                ;
                                matrix _BETA_ {4}
                                ;
                                matrix _GAMMA_ {5}
                                ;
                                matrix _PHI_ {6}
                                ;
                            run;

                        data SEM.gene_{7}_model_{8};
                            length gene $12.;
                            length model $14.;
                            retain gene model;
                            set fitstat;
                            where IndexCode = 312;
                            gene = "{7}";
                            model = "Model {8}";
                            rename FitValue = BIC;
                            keep gene model FitValue;
                            run;
                        """.format(args.lname, args.mname, yvar, xvar, BETA, GAMMA, PHI, args.gname, modelNum))

def flatten_list(myList):
    """ Isoforms are being treated together as a group. To delineate them,
    isoforms are grouped together as tuples. For the calis statement, xvar and
    yvar need to each be a string. This function flattens these tuples into a
    single list and creates a string. """

    result = []
    for element in myList: 
        if hasattr(element, '__iter__'):
            # if tuple flatten and append
            result.extend(element)
        else:
            # if not a tuple just append
            result.append(element)
    return(' '.join(result))

def generate_sas(yvar, xvar, BETA, GAMMA, PHI, args, fbase=0):
    """ Build the SAS output """
    # Flatten variable list into a string
    _yvar_ = flatten_list(yvar)
    _xvar_ = flatten_list(xvar)

    # Create SAS proc calis matrix statements
    _BETA_ = build_matrix_output(BETA,'beta')
    _GAMMA_ = build_matrix_output(GAMMA,'gamma')
    il1 = np.triu_indices_from(PHI)
    PHI[il1] = 2
    _PHI_ = build_matrix_output(PHI,find_zero=1)

    # Create SAS proc calis statement and write it out
    build_calis(_yvar_, _xvar_, _BETA_, _GAMMA_, _PHI_, args, fbase)

# Function to Build log
def output_model_to_log(yvar, xvar, BETA, GAMMA, PHI, model_type):
    """ Create output for the log """
    # Bind global instance of myCount variable.
    global myCount

    # Construct Output
    if model_type == 'Baseline':
        # Set model number to 0 for the log
        message = "\nModel number {0}\nModel type: {1}\nY-variables: {2}\nX-variables: {3}\n\nBeta:\n{4}\n\nGAMMA:\n{5}\n\nPHI:\n{6}\n\n".format(0, model_type, yvar, xvar, BETA, GAMMA, PHI)
    else:
        message = "\nModel number {0}\nModel type: {1}\nY-variables: {2}\nX-variables: {3}\n\nBeta:\n{4}\n\nGAMMA:\n{5}\n\nPHI:\n{6}\n\n".format(myCount, model_type, yvar, xvar, BETA, GAMMA, PHI)
    logging.info(message)
