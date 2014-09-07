#!/usr/bin/env python
import logging
from os.path import splitext
import numpy as np
from semnet.writer import createOutput
from semnet import utils

################################################################################
# Construct baseline model
################################################################################
def baseline(path, args):
    """ Create out for running the baseline model. """
    model_type = "Baseline"
    createOutput(path, model_type, args)
