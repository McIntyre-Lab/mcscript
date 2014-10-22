import logging
import sys

def set_logger(fname=None,logLevel='info'):
    """ Function to set up the handle error logging.
    fname = the name of the log file, if no logfile is given then logging
            information will be pritned to stdout.

    logLevel = level of information to print out, options are {info, debug} [Default: info]
    """

    # Determine log level
    if logLevel == 'info':
        _log = logging.INFO
    elif logLevel == 'debug':
        _log = logging.DEBUG

    # Output to file or STDOUT
    if fname:
        logging.basicConfig(filename=fname, filemode='w', level=_log, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(stream=sys.stderr, level=_log, format='%(asctime)s - %(levelname)s - %(message)s')
