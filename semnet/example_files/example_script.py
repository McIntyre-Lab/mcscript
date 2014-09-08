#!/usr/bin/env python
# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys

# Add-on packages
import numpy as np

# McLab Packages
import mclib
import semnet

def getOptions():
    """ Function to pull in arguments """
    description = """"This script starts with a begining pathway in the format
    of a 'path' file. It then builds SAS code for SEMs after adding a new gene
    to all possible locations within a network.

        (a) Downstream of endogenous genes
        (b) Downstream of exogenous genes
        (c) Upstream of endogenous genes
        (d) Upstream of exogenous genes
        (e) In between exogenous genes and their endogenous targets
        (e) In betwee endogenous genes and their endogenous targets
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-p", dest="pname", action='store', required=True, help="Name of the 'path' file [Required]")
    parser.add_argument("-l", dest="lname", action='store', required=True, help="Path to sas library that has the SAS dataset that will be analyzed [Required]")
    parser.add_argument("-m", dest="mname", action='store', required=True, help="Name of SAS dataset that will be analyzed [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name with PATH of the output file ending with '.sas'; NOTE: model number will be appended to the filename [Required]")
    parser.add_argument("-n", dest="newGene", nargs='+', action='store', required=True, help="List of new gene isoforms separted by spaces [Required]")
    parser.add_argument("-g", dest="gname", action='store', required=True, help="The gene name of the new gene/isoforms being added [Required]")
    parser.add_argument("-t", dest="template", action='store', required=False, help="Name of the PROC CALIS template file [Optional]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 
    args = parser.parse_args()
    return(args)

if __name__ == '__main__':
    args = getOptions()

    # Turn on logging
    mclib.logger.set_logger(args.log, info)
    mclib.git.git_to_log()

    # Initialize base variable list
    path = semnet.createPath(args.pname)
    newGene = semnet.NewGene(args.newGene)

    # Run baseline model
    semnet.models_baseline.baseline(path, args)

    # Add genes to all locations in network
    ## Downstream
    semnet.models_add.add_genes_ds_beta(path, newGene, args)
    semnet.models_add.add_genes_ds_gamma(path, newGene, args)

    ## Upstream
    semnet.models_add.add_genes_above_beta(path, newGene, args)
    semnet.models_add.add_genes_above_gamma(path, newGene, args)

    ## In between
    semnet.models_add.add_genes_bt_beta(path, newGene, args)
    semnet.models_add.add_genes_bt_gamma_beta(path, newGene, args)

    ## Add new links for gene existing in the network
    #semnet.models_newlink.add_newlinks_beta_to_beta(path, args)
    #semnet.models_newlink.add_newlinks_gamma_to_beta(path, args)
    #semnet.models_newlink.add_newlinks_beta_to_gamma(path, args)
    #semnet.models_newlink.add_newlinks_gamma_to_gamma(path, args)
