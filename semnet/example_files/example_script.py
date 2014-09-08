#!/usr/bin/env python
import argparse
import logging
import sys
import numpy as np
import semnet

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="A script that adds genes through out a network and generates sas models.")
    parser.add_argument("-p", dest="pname", action='store', required=True, help="Name of the PATH file [Required]")
    parser.add_argument("-l", dest="lname", action='store', required=True, help="Path to sas library that has the original data [Required]")
    parser.add_argument("-m", dest="mname", action='store', required=True, help="Name of SAS dataset [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Output file ending with .sas, note model number will be appended to the filename. [Required]")
    parser.add_argument("-n", dest="newGene", nargs='+', action='store', required=True, help="List of new gene isoforms separted by spaces [Required]")
    parser.add_argument("-g", dest="gname", action='store', required=True, help="The gene name [Required]")
    parser.add_argument("-t", dest="template", action='store', required=False, help="Name of the calis template file [Optional]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    #args = parser.parse_args(['-l', 'bob', '-o','/home/jfear/tmp/CG10036.sas','-n','CG10036_PC CG10036_PE', '-g', 'CG10036', '--log', '/home/jfear/tmp/CG10036.log'])
    return(args)


def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, filemode='w', level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')


if __name__ == '__main__':
    args = getOptions()
    if args.log:                                         # Turn on Logging if option -g was given
        setLogger(args.log,logging.INFO)

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

    ## New Links to Betas
    semnet.models_newlink.add_newlinks_beta_to_beta(path, args)
    semnet.models_newlink.add_newlinks_gamma_to_beta(path, args)
    semnet.models_newlink.add_newlinks_beta_to_gamma(path, args)
    semnet.models_newlink.add_newlinks_gamma_to_gamma(path, args)
