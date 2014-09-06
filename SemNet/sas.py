#!/usr/bin/env python
from os.path import splitext
import numpy as np

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

if __name__ == '__main__':
    pass
