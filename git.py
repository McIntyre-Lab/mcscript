import logging
import subprocess
import os.path

def get_git():
    """ Parse the current Git commit version. This will
    allow recording of exactly which script was used to create a given
    output.
    """
    # get full path to script
    fullname = os.path.abspath(__file__)
    gitdir = os.path.dirname(fullname)
    label = subprocess.check_output(["git", "log", "-n1", "--pretty=%h", fullname])
    return(label.rstrip(), gitdir)

def git_to_log():
    """ Write current git commit information to the log. """
    try:
        git_status, gitdir = get_git()
        logging.info("Starting %s", __file__) 
        logging.info("Running script from  %s", gitdir) 
        logging.info("Git commit id: %s", git_status)
    except:
        pass
