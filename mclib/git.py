import logging
import subprocess
import os.path

def get_git(sname):
    """ Parse the current Git commit version. This will
    allow recording of exactly which script was used to create a given
    output.
    """
    # get full path to script
    fullname = os.path.abspath(sname)
    gitdir = os.path.dirname(fullname)
    label = subprocess.check_output(["git", "--git-dir="+gitdir+"/.git", "--work-tree="+gitdir,"log", "-n1", "--pretty=%h", fullname])

    return(label.rstrip(), gitdir)

def git_to_log(sname):
    """ Write current git commit information to the log. """
    try:
        git_status, gitdir = get_git(sname)
        logger.info("Starting %s", sname) 
        logger.info("Running script from  %s", gitdir) 
        logger.info("Git commit id: %s", git_status)
    except:
        print "You need access to the repository to output git information."
        pass
