#!/usr/bin/env python
import numpy as np
from string import Template
from semnet import utils

class CalisOut(object):
    def __init__(self, path, model, libname, sasdat, gene, template='./templates/calis.sas'):
        self.model = model              # model number
        self.libname = libname          # sas libname to use in sas program
        self.sasdat = sasdat            # input sasdataset for the program
        self.gene = gene                # gene name
        self.template = template        # Get SAS code template from 

        # Flatten variable list into a string for proc calis
        self._yvar = utils.flatten_list(path.yvar)
        self._xvar = utils.flatten_list(path.xvar)

        # Create SAS proc calis matrix statements
        self._beta = self.build_matrix_output(path.beta,'beta')
        self._gamma = self.build_matrix_output(path.gamma,'gamma')
        il1 = np.triu_indices_from(path.phi)
        path.phi[il1] = 2
        self._phi = self.build_matrix_output(path.phi,find_zero=1)

        # Create SAS proc calis program
        self.build_calis()

    def build_matrix_output(self, myMatrix, prefix='', find_zero=0):
        """ Construct proc calis matrix statements from full matrices. Iterate over
        the full matrix and identify the important coordinates and construct a
        string with the correct syntax.
            
        myMatrix is a matrix of {1,0} either Beta, Gamma, or Phi from SEM.

        prefix {beta, gamma, ''} is the prefix for the parameter labels, leave
        blank if none.

        find_zero {0,1} if 0 then identify non-zero indices (for beta and gamma
        matrices), if 1 identify 0 indices (for covariance matrix phi).
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
        for index, (row, col) in enumerate(zip(rows, cols)):
            if prefix:
                myList.append("[{0},{1}] = {2}{3}".format(row+1,col+1,prefix,index+1))
            else:
                myList.append("[{0},{1}] = 0".format(row+1,col+1))
        return(",\n".join(myList))

    def build_calis(self):
        """ Create the proc calis statement for sas. Uses a template file to
        make the substitutions. 
        """
        FH = open(self.template, 'r')
        srcTemp = Template(FH.read())
        outDic = { 'libname': self.libname, 'data': self.sasdat, 
                   'yvar': self._yvar, 'xvar': self._xvar, 'beta': self._beta, 
                   'gamma': self._gamma, 'phi': self._phi, 'gene': self.gene, 
                   'model': self.model 
                }
        self.out = srcTemp.substitute(outDic)

if __name__ == '__main__':
    pass
