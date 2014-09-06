#!/usr/bin/env python
import logging
from os.path import splitext
import numpy as np

# Functions to Add Gene to different locations in the network
def run_baseline(yvar, xvar, BETA, GAMMA, PHI, newGene, newCount, args):
    """ Iterate through beta matrix and place new genes downstream of each beta. """
    model_type = "Baseline"

    # Bind global instance of myCount variable.
    global myCount

    generate_sas(yvar, xvar, BETA, GAMMA, PHI, args, fbase=1)
    output_model_to_log(yvar, xvar, BETA, GAMMA, PHI, model_type)
