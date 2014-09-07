from collections import defaultdict
import numpy as np
import semnet

def rowSplitter(row, separator):
    """ Split the row of a path file based on a give separator. """
    rowSplit = [i.strip(' ') for i in row.strip().split(separator)]
    left = isoSplitter(rowSplit[0])
    right = isoSplitter(rowSplit[1])
    return(left, right)

def isoSplitter(lst, separator='|'):
    """ Split isoforms that are concatenated with '|' """
    genes = lst.split(' ')
    for index, gene in enumerate(genes):
        if '|' in gene:
            genes[index] = tuple(gene.split('|'))
    return genes

def parsePathFile(fname):
    """ Iterate over a 'Path' file. Determine which genes covary and which
    genes are sources and targets.

    gene 1 <-> gene2: indicates two genes covary

    gene 1 -> gene2 gene3: indicates gene1 regulates gene2 and gene3
    gene 1 <- gene2 gene3: indicates gene2 and gene3 regulate gene1
    """
    try:
        paths = defaultdict(list)
        covs = []
        with open(fname) as FH:
            for row in FH:
                if '<->' in row:
                   left, right = rowSplitter(row, '<->')
                   covs.append(tuple(left+right))
                elif '->' in row:
                   left, right = rowSplitter(row, '->')
                   for i in left:
                       paths[i].extend(right)
                elif '<-' in row:
                   left, right = rowSplitter(row, '<-')
                   for i in right:
                       paths[i].extend(left)
        return(paths, covs)
    except IOError:
        print "The file {0} does not exist.".format(fname)

def identifyEndogenous(paths):
    """ Endogenous variables are those that have an arrow going into them. So
    they are the values of the paths dictionary. I need to create a list of
    these endogenous variables.
    """
    yvar = set()
    for genes in paths.values():
        for gene in genes:
            yvar.add(gene)
    return sorted(yvar)

def identifyExogenous(yvar, paths):
    """ Exogenous variables are those that do not have an arrow going into
    them. So they are the keys in the paths dictionary that are not also values
    (ie not in yvar).
    """
    xvar = set()
    for key in paths:
        if key not in yvar:
            xvar.add(key)
    return sorted(xvar)

def buildBeta(yvar, paths):
    """ Beta is a square matrix where the rows and columns are equal to the
    endogenous genes. In other words, Beta is a 0,1 matrix of the relationships
    between endogenous (yvar) genes.
    """
    # flatten the yvar list to remove isoform grouping
    fyvar = semnet.utils.flatten_list(yvar, 1)

    # Initialize Beta as all 0's
    nyvar = len(fyvar)
    Beta = np.zeros((nyvar, nyvar))

    # Iterate through yvar and see if there are any relationships
    for colIndex, value in enumerate(yvar):
        # If the current gene has a downstream target continue
        if paths[value]:
            # get isoform count
            if hasattr(value, '__iter__'):
                isoCnt = len(value)
            else:
                isoCnt = 1

            # get flattened list of downstream targets
            targets = semnet.utils.flatten_list(paths[value],1)

            # Add ones to correct spot in Beta
            for target in targets:
                rowIndex = fyvar.index(target)
                Beta[rowIndex, colIndex:colIndex+isoCnt] = 1
    return Beta

def buildGamma(yvar, xvar, paths):
    """ Gamma has 'yvar' rows and 'xvar' columns.  Gamma is a 0,1 matrix of the
    relationships between exogenous (xvar) to endogenous (yvar) genes.
    """

    # flatten yvar and xvar to remove isoform grouping
    fyvar = semnet.utils.flatten_list(yvar, 1)
    fxvar = semnet.utils.flatten_list(xvar, 1)

    # Initialize Gamma as all 0's
    nyvar = len(fyvar)
    nxvar = len(fxvar)
    Gamma = np.zeros((nyvar, nxvar))

    # Iterate through xvar and see if there are any relationships
    for colIndex, value in enumerate(xvar):

        # If the current gene has a downstream target continue
        if paths[value]:
            # get isoform count
            if hasattr(value, '__iter__'):
                isoCnt = len(value)
            else:
                isoCnt = 1

            # get flattened list of downstream targets
            targets = semnet.utils.flatten_list(paths[value],1)

            # Add ones to correct spot in Gamma
            for target in targets:
                rowIndex = fyvar.index(target)
                Gamma[rowIndex, colIndex:colIndex+isoCnt] = 1
    return Gamma

def buildPhi(xvar, covs):
    """ Phi is the 0,1 variance covariance matrix of the exogenous variables
    (xvar). So it is a square matrix where the number or rows and columns is
    equal to the number of xvars.
    """

    # flatten yvar and xvar to remove isoform grouping
    fxvar = semnet.utils.flatten_list(xvar, 1)

    # Initialize Phi as all 0's with 1's along the diagonal.
    nxvar = len(fxvar)
    Phi = np.eye(nxvar)

    # Iterate through the covariance list. Find each genes in the pair is
    # located in phi and change covariances to 1. Remember that Phi is a
    # symmetric matrix so you need to change both locations
    for cov in covs:
        coord = []
        for gene in cov:
            coord.append(fxvar.index(gene))

        Phi[coord[0], coord[1]] = 1
        Phi[coord[1], coord[0]] = 1
    return Phi

def createPath(fname):
    paths, covs = parsePathFile(fname)
    yvar = identifyEndogenous(paths)
    xvar = identifyExogenous(yvar, paths)
    beta = buildBeta(yvar, paths)
    gamma = buildGamma(yvar, xvar, paths)
    phi = buildPhi(xvar, covs)
    graph = semnet.SemPath(xvar, yvar, beta, gamma, phi)
    return graph

if __name__ == '__main__':
    pass
