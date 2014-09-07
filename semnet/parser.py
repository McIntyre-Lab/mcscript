from collections import defaultdict
import numpy as np
import semnet

def rowSplitter(row, separator):
    """ Split the row of a path file based on a give separator. """
    rowSplit = [i.strip(' ') for i in row.strip().split(separator)]
    left = rowSplit[0].split(' ')
    right = rowSplit[1].split(' ')
    return(left, right)

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
    # Initialize Beta as all 0's
    nyvar = len(yvar)
    Beta = np.zeros((nyvar, nyvar))

    # Iterate through yvar and see if there are any relationships
    for colIndex, value in enumerate(yvar):
        if paths[value]:
            targets = paths[value]
            for target in targets:
                rowIndex = yvar.index(target)
                Beta[rowIndex, colIndex] = 1
    return Beta

def buildGamma(yvar, xvar, paths):
    """ Gamma has 'yvar' rows and 'xvar' columns.  Gamma is a 0,1 matrix of the
    relationships between exogenous (xvar) to endogenous (yvar) genes.
    """
    # Initialize Gamma as all 0's
    nyvar = len(yvar)
    nxvar = len(xvar)
    Gamma = np.zeros((nyvar, nxvar))
    # Iterate through xvar and see if there are any relationships
    for colIndex, value in enumerate(xvar):
        try:
            targets = paths[value]
            for target in targets:
                rowIndex = yvar.index(target)
                Gamma[rowIndex, colIndex] = 1
        except:
            pass
    return Gamma

def buildPhi(xvar, covs):
    """ Phi is the 0,1 variance covariance matrix of the exogenous variables
    (xvar). So it is a square matrix where the number or rows and columns is
    equal to the number of xvars.
    """
    # Initialize Phi as all 0's with 1's along the diagonal.
    nxvar = len(xvar)
    Phi = np.eye(nxvar)

    # Iterate through the covariance list. Find each genes in the pair is
    # located in phi and change covariances to 1. Remember that Phi is a
    # symmetric matrix so you need to change both locations
    for cov in covs:
        coord = []
        for gene in cov:
            coord.append(xvar.index(gene))

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
