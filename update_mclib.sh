#!/bin/bash - 
#===============================================================================
#
#   DESCRIPTION: mclib_Python is a subtree to a seperate repository. Run this
#   script to pull the most recent upstream changes from mclib_Python and merge
#   them into mcscript.
# 
#===============================================================================

set -o nounset                              # Treat unset variables as an error
set -e

git subtree pull --prefix=mclib_Python --squash mclib_Python master
