import logging
import sys

def setLogger(logger, fname=None, logLevel='info'):
    """ Function to set up the handle error logging.
    logger (obj) = a logger object

    fname = the name of the log file, if no logfile is given then logging
            information will be pritned to stdout.

    logLevel = level of information to print out, options are {info, debug} [Default: info]
    """

    # Determine log level
    if logLevel == 'info':
        _level = logging.INFO
    elif logLevel == 'debug':
        _level = logging.DEBUG

    # Set the level in logger
    logger.setLevel(_level)

    # Set the log format
    logfmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Output to file or STDOUT
    if fname:
        fh = logging.FileHandler(filename=fname,mode='w')
        fh.setLevel(_level)
        fh.setFormatter(logfmt)
        logger.addHandler(fh)
    else:
        # Since no file was provide output to STDOUT and STDERR
        loghandler = logging.StreamHandler(stream=sys.stdout)
        errhandler = logging.StreamHandler(stream=sys.stderr)

        # Set logging level for the different output handlers. 
        # ERRORs to STDERR, and everything else to STDOUT
        loghandler.setLevel(_level)
        errhandler.setLevel(logging.ERROR)

        # Format the log handlers
        loghandler.setFormatter(logfmt)
        errhandler.setFormatter(logfmt)

        # Add handler to the main logger
        logger.addHandler(loghandler)
        logger.addHandler(errhandler)
