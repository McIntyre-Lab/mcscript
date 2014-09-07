import numpy as np
import logging

class SemPath(object):
    """ A data storage class with everything you need to run a SEM. """
    def __init__(self, xvar=None, yvar=None, beta=None, gamma=None, phi=None):
        self.xvar = xvar
        self.yvar = yvar
        self.beta = beta
        self.gamma = gamma
        self.phi = phi

        # Make a private copy of the original yvar, xvar, beta, gamma and phi
        # matrices so that I can reset after adding a gene.
        self._yvar = list(yvar)
        self._xvar = list(xvar)
        self._beta = beta.copy()
        self._gamma = gamma.copy()
        self._phi = phi.copy()

        self.bRow, self.bCol = np.shape(beta)
        self.gRow, self.gCol = np.shape(gamma)
        self.pRow, self.pCol = np.shape(phi)

        # Initialize a model counter
        count = 1

    def __repr__(self):
        """ Simple method to print out the information stored in the __init__ """
        return "\n\nYvar: {0}\nXvar: {1}\n\nBeta:\n{2}\n\nGamma:\n{3}\n\nPhi:\n{4}\n\n".format(self.yvar, self.xvar, self.beta, self.gamma, self.phi)

    def reinit(self):
        """ Re-initialize the matrices. The adding genes process changes the
        original beta, gamma, phi matrices, this function re-initializes them
        to the original state.
        """
        self.yvar = list(self._yvar)
        self.xvar = list(self._xvar)
        self.beta = self._beta.copy()
        self.gamma = self._gamma.copy()
        self.phi = self._phi.copy()

    def count_increment(self):
        self.count += 1

    def expand_beta(self):
        """ Expands the self._beta matrix by adding a row of 0's and a column
        of 0's.
        """
        self.beta = np.insert(self.beta, self.bRow, values=0, axis=1)
        self.beta = np.insert(self.beta, self.bRow, values=0, axis=0)

    def expand_gamma(self, axis=0):
        """ Expands the self._gamma matrix by adding a row of 0's when axis=0.
        Or by adding a column of 0's when axis=1.
        """
        gPos = np.shape(self.gamma)[axis]
        if axis == 0:
            self.gamma = np.insert(self.gamma, self.gRow, values=0, axis=axis)
        elif axis == 1:
            self.gamma = np.insert(self.gamma, self.gCol, values=0, axis=axis)
        else:
            print "Axis can only be 0 or 1."
            raise ValueError

    def expand_phi(self):
        """ Expands the self._phi matrix by adding a row of 0's and a column
        of 0's. Then makes the diagonal 1.
        """
        self.phi = np.insert(self.phi, self.pRow, values=0, axis=0)
        self.phi = np.insert(self.phi, self.pRow, values=0, axis=1)

    def del_col_gamma(self, coords):
        """ Uses numpy to delete columns from the gamma matrix. coords is a
        tuple with the start and end coordinates for the columns you want
        delete.
        """
        self.gamma = np.delete(self.gamma, range(*coords), 1)

    def del_row_col_phi(self, coords):
        """ Uses numpy to delete columns and rows from the phi matrix. coords
        is a tuple with the start and end coordinates for the columns you want
        delete.
        """
        self.phi = np.delete(self.phi, range(*coords), 0)
        self.phi = np.delete(self.phi, range(*coords), 1)

class NewGene(object):
    """ When adding genes, users provide a list of genes/isoforms to add. It is
    important to expand the matrices using this isoform information.
    """
    def __init__(self, newGene):

        try:
            # Determine if the gene we are adding to the model has multiple isoforms
            # and create a counter
            self.name = newGene

            if hasattr(self.name,'__iter__'):
                self.count = len(self.name)
            else:
                self.count = 1
        except ValueError:
            logging.error("If you are adding genes you need to create a new gene object.")

    def __repr__(self):
        return "\n\nNew Gene: {0}\nNumber of isoforms: {1}\n\n".format(self.name, self.count)

if __name__ == '__main__':
    pass
