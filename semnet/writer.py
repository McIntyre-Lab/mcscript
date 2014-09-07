#!/usr/bin/env python
import logging
from os.path import splitext
from semnet.sas import CalisOut

def output_model_to_log(path, model_type, count):
    """ Create output for the log """
    # Construct Output
    message = "\nModel number {0}\nModel type: {1}{2}".format(count, model_type, str(path)) 
    logging.info(message)

def createOutput(path, model_type, args):

    # If dealing with a baseline model, set the count to 0
    if model_type == 'Baseline':
        count = 0
        model = model_type
    else:
        count = path.count
        model = path.count

    # If the user supplied a template file then use it.
    if args.template:
        calis = CalisOut(path, model_type, args.lname, args.mname, args.gname, args.template)
    else:
        calis = CalisOut(path, model_type, args.lname, args.mname, args.gname)

    # Create a output sas file
    oname = "{0}_{1}.sas".format(splitext(args.oname)[0], model)
    with open(oname, 'w') as OUT:
        OUT.write(calis.out)

    # write log output
    output_model_to_log(path, model_type, count)

if __name__ == '__main__':
    pass
