#!/usr/bin/env python
import numpy as np

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
