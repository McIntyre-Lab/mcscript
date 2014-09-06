#!/usr/bin/env python
import logging

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
